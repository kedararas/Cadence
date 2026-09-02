function mv_substrate_test
%MV_SUBSTRATE_TEST  Regression test for assess_arrhythmia_substrate.
%
%   mv_substrate_test
%
%   Runs the real path — synthetic stack -> analyzeAPAlternans (3-D) ->
%   assess_arrhythmia_substrate — and checks the two things changed on
%   2026-09-01:
%
%   1. THE APD GRADIENT IS NOT THE MASK BORDER.  The previous version zeroed
%      the background before the Sobel and then masked only the OUTSIDE, so
%      every tissue pixel next to background kept a ~100 ms -> 0 step.  On a
%      synthetic map with a true 0.32 ms/pixel gradient that reported
%      max_gradient = 467 ms/pixel against a true interior 0.64 — a 730x
%      inflation, which also inflated the 95th percentile normalising the
%      risk component.
%
%   2. RESTITUTION IS NOT IN THE RISK COMPOSITE.  Under fixed-CL pacing the
%      slope is 1 by algebra (DI = CL - APD, so APD_next is an exact linear
%      function of DI), so `restitution_map > 1` scored rounding noise while
%      carrying a full 1.0 of a 4.0 maximum.  It is still returned for
%      inspection; it must not be scored.
%
%   Not a study instrument — a pass/fail run, like mv_rise_test.

    fs = 1000; dt = 1000/fs;
    R  = 48; C = 48;
    cl = 200;                       % pacing cycle length, ms
    nBeats = 8;
    T  = cl * nBeats;

    % Round tissue so the mask has a border for the erosion to matter.
    [rr, cc] = ndgrid(1:R, 1:C);
    tissue   = hypot(rr - (R+1)/2, cc - (C+1)/2) < 18;

    % True APD: a gentle gradient across columns, plus 2:1 alternans.
    grad_per_px = 0.30;                                  % ms per pixel
    apd_base    = 100 + grad_per_px * (cc - 1);
    alt_amp     = 6;                                     % ms peak-to-peak

    t  = (0:T-1) * dt;
    X  = zeros(R, C, T);
    rng(11, 'twister');
    for b = 1:nBeats
        ta  = (b-1)*cl + 8;
        apd = apd_base + (alt_amp/2) * (-1)^b;
        for r = 1:R
            for c = 1:C
                X(r,c,:) = squeeze(X(r,c,:))' + ap_wave(t, ta, apd(r,c));
            end
        end
    end
    X = X + 0.02 * randn(R, C, T);

    mask = ones(R, C); mask(~tissue) = NaN;
    bf = zeros(nBeats, 2);
    for b = 1:nBeats
        bf(b,1) = (b-1)*cl/dt + 1;
        bf(b,2) = min(T, b*cl/dt);
    end

    fprintf('MV_SUBSTRATE_TEST  %dx%d, %d beats, true APD gradient %.2f ms/pixel\n', ...
            R, C, nBeats, grad_per_px);

    ap = analyzeAPAlternans(1:T, X, 'BeatFrames', bf, 'Mask', mask, 'PlotResults', false);
    sub = assess_arrhythmia_substrate(ap, 'BeatFrames', bf);

    % ---- 1. gradient ----------------------------------------------------
    % The absolute value cannot be asserted against the ramp, because the APD
    % map carries real measurement noise: at APD90 the crossing sits on a
    % shallow part of the decay, so sigma_t = sigma/|slope| is several ms per
    % pixel and the Sobel amplifies it.  That noise is a property of the data,
    % not of the gradient code.  What IS assertable is the border: reproduce
    % the old handling on the same map and confirm the ring no longer sets the
    % answer, while the interior is left alone.
    m = sub.apd_mean_map;
    vmask = isfinite(m);
    a_old = m; a_old(~vmask) = 0;
    g_old = imgradient(a_old, 'sobel'); g_old(~vmask) = nan;

    max_old = max(g_old(:), [], 'omitnan');
    p95_old = prctile(g_old(isfinite(g_old)), 95);
    p95_new = prctile(sub.gradient_mag(isfinite(sub.gradient_mag)), 95);
    med_old = median(g_old(isfinite(g_old)));
    med_new = median(sub.gradient_mag(isfinite(sub.gradient_mag)));

    fprintf('\n  max_gradient  old %8.2f -> new %8.2f ms/pixel\n', max_old, sub.max_gradient);
    fprintf('  p95 (risk norm) old %8.2f -> new %8.2f\n', p95_old, p95_new);
    fprintf('  median          old %8.2f -> new %8.2f  (interior should not move)\n', med_old, med_new);

    check('border no longer sets max_gradient', sub.max_gradient < max_old / 5, ...
          sprintf('%.1f vs %.1f, %.0fx lower', sub.max_gradient, max_old, max_old/sub.max_gradient));
    check('border no longer sets the risk p95', p95_new < p95_old / 5, ...
          sprintf('%.1f vs %.1f', p95_new, p95_old));
    check('interior gradient is preserved', abs(med_new - med_old) < 0.25 * med_old, ...
          sprintf('median %.2f vs %.2f', med_new, med_old));
    check('border ring is excluded', nnz(isfinite(sub.gradient_mag)) < nnz(isfinite(g_old)), ...
          sprintf('%d pixels dropped', nnz(isfinite(g_old)) - nnz(isfinite(sub.gradient_mag))));

    % ---- 2. restitution excluded from the composite ---------------------
    check('restitution_map still returned',  isfield(sub, 'restitution_map'), 'present');
    check('slope_gt1_mask still returned',   isfield(sub, 'slope_gt1_mask'), 'present');

    % The slope really is ~1 by construction here, so the old code would have
    % scored a large fraction of pixels. Show it, then prove it is not scored.
    sg = sub.slope_gt1_mask & isfinite(sub.restitution_map);
    fprintf('  restitution slope   : median %.4f over %d valid pixels\n', ...
            median(sub.restitution_map(isfinite(sub.restitution_map))), ...
            nnz(isfinite(sub.restitution_map)));
    fprintf('  slope_gt1 pixels    : %d (would have scored 1.0 each under the old composite)\n', nnz(sg));

    % Re-derive the composite from the returned components and confirm the
    % denominator is 3 (no Ca) and that no restitution term is present.
    valid = isfinite(sub.apd_mean_map);
    c1 = min(sub.alt_ratio_map / 0.20, 1); c1(isnan(c1) | ~valid) = 0;
    gp = prctile(sub.gradient_mag(~isnan(sub.gradient_mag)), 95);
    c2 = min(sub.gradient_mag / (gp + eps), 1); c2(isnan(c2)) = 0;
    c3 = double(sub.nodal_lines);
    expect_no_rest = (c1 + c2 + c3) / 3.0;
    expect_with_rest = (c1 + c2 + c3 + double(sub.slope_gt1_mask)) / 4.0;

    d_no   = max(abs(sub.risk_map(valid) - expect_no_rest(valid)));
    d_with = max(abs(sub.risk_map(valid) - expect_with_rest(valid)));
    fprintf('  |risk - (c1+c2+c3)/3|          : %.3g\n', d_no);
    fprintf('  |risk - (c1+c2+c3+c4)/4| (old) : %.3g\n', d_with);
    check('risk composite excludes restitution', d_no < 1e-9, sprintf('max diff %.3g', d_no));
    check('risk composite is not the old form',  d_with > 1e-6 || nnz(sg) == 0, ...
          sprintf('max diff %.3g with %d slope_gt1 pixels', d_with, nnz(sg)));
    check('risk_map stays in [0,1]', ...
          all(sub.risk_map(valid) >= -1e-12 & sub.risk_map(valid) <= 1 + 1e-12), ...
          sprintf('range %.3f-%.3f', min(sub.risk_map(valid)), max(sub.risk_map(valid))));

    fprintf('\nAll assertions passed.\n');
end


function y = ap_wave(t, ta, apd80)
    tu  = ta + 2;
    tau = (apd80 - 2) / log(5);
    y   = zeros(size(t));
    r = t >= ta & t < tu;  y(r) = (t(r) - ta) / 2;
    d = t >= tu;           y(d) = exp(-(t(d) - tu) / tau);
end


function check(name, cond, detail)
    if cond
        fprintf('   PASS  %-40s (%s)\n', name, detail);
    else
        error('mv_substrate_test:failed', 'FAIL  %s (%s)', name, detail);
    end
end
