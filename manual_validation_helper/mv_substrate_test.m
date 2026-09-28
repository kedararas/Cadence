function mv_substrate_test
%MV_SUBSTRATE_TEST  Regression test for assess_arrhythmia_substrate.
%
%   mv_substrate_test
%
%   Runs the real path — synthetic stack -> analyzeAPAlternans (3-D) ->
%   assess_arrhythmia_substrate — and checks the two things changed on
%   2026-09-01, plus the phase/mask level fix of 2026-09-28 (3):
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
%   3. PHASE AND SIGNIFICANCE USE THE SAME APD LEVEL.  The phase map used to
%      come from ap_result (pinned to APD80) while sig_mask used the
%      auto-selected level (APD50 at fast pacing).  Significant pixels that
%      were NaN at APD80 counted in n_sig but in neither sign, so
%      concordance_ratio dropped below its 0.5 floor (0.441 on the Fig 4
%      recording) and is_discordant was over-called.  Checked on a hand-built
%      ap_result where APD80 is missing on 40% of the tissue and opposite in
%      sign on another 15%, once with concordant and once with genuinely
%      discordant APD50 alternans.  concordance_ratio must lie in [0.5, 1].
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
    % /8 as the function now does (Sobel normalisation, 2026-09-04), so old and
    % new border handling are compared in the same ms/pixel units.
    g_old = imgradient(a_old, 'sobel') / 8; g_old(~vmask) = nan;

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

    % ---- 3. phase and significance at the same APD level ----------------
    check('concordance_ratio in [0.5, 1] (real path)', ...
          sub.concordance_ratio >= 0.5 && sub.concordance_ratio <= 1, ...
          sprintf('%.3f', sub.concordance_ratio));

    fprintf('\n  phase / mask level (APD80 missing on 40%%, reversed on 15%%):\n');
    % Concordant APD50: every pixel alternates the same way.
    ap = fast_rate_result(false);
    s3 = assess_arrhythmia_substrate(ap);
    old_ratio = stale_ratio(ap);
    fprintf('  concordant   ratio %.3f (APD80-phase ratio would be %.3f)\n', ...
            s3.concordance_ratio, old_ratio);
    check('APD80 phase would breach the floor', old_ratio < 0.5, sprintf('%.3f', old_ratio));
    check('concordance_ratio in [0.5, 1]', ...
          s3.concordance_ratio >= 0.5 && s3.concordance_ratio <= 1, sprintf('%.3f', s3.concordance_ratio));
    check('auto level falls back to APD50',   s3.apd_level == 50, sprintf('APD%d', s3.apd_level));
    check('uniform APD50 alternans is concordant', ...
          s3.concordance_ratio == 1 && ~s3.is_discordant, sprintf('%.3f', s3.concordance_ratio));
    check('phase_map is sign(APD50_alt)', ...
          isequal(s3.phase_map, sign(ap.APD50_alt)), 'all pixels');
    check('no nodal lines without a sign change', ~any(s3.nodal_lines(:)), ...
          sprintf('%d px', nnz(s3.nodal_lines)));
    check('phase angle defined where APD50 is', all(isfinite(s3.phase_angle_map(:))), ...
          sprintf('%d NaN', nnz(~isfinite(s3.phase_angle_map))));

    % Discordant APD50: the fix must not simply force the ratio to 1.
    ap = fast_rate_result(true);
    s3 = assess_arrhythmia_substrate(ap);
    fprintf('  discordant   ratio %.3f (APD80-phase ratio would be %.3f)\n', ...
            s3.concordance_ratio, stale_ratio(ap));
    check('concordance_ratio in [0.5, 1]', ...
          s3.concordance_ratio >= 0.5 && s3.concordance_ratio <= 1, sprintf('%.3f', s3.concordance_ratio));
    check('30/70 APD50 split is 0.70 and discordant', ...
          abs(s3.concordance_ratio - 0.70) < 1e-12 && s3.is_discordant, sprintf('%.3f', s3.concordance_ratio));
    check('nodal line on the APD50 reversal only', ...
          isequal(find(any(s3.nodal_lines, 1)), [12 13]), ...
          sprintf('columns %s', mat2str(find(any(s3.nodal_lines, 1)))));

    fprintf('\nAll assertions passed.\n');
end


function ap = fast_rate_result(discordant)
%FAST_RATE_RESULT  Minimal 3-D ap_result shaped like fast pacing: APD30 and
%   APD50 are finite everywhere, APD80 is NaN on columns 1-16 (coverage 60%,
%   so 'auto' picks APD50) and, where it exists, reversed in sign on columns
%   17-22 — the APD50/APD80 sign disagreement seen on real data.  ap.phase_map
%   is sign(APD80_alt), as analyzeAPAlternans stores it.
%   discordant = true reverses APD50 alternans on columns 1-12 (30%).
    R = 40; C = 40; nB = 8; lv = [30 50 80];
    s50 = ones(R, C);
    if discordant, s50(:, 1:12) = -1; end
    s80 = s50;  s80(:, 17:22) = -s80(:, 17:22);
    m80 = true(R, C);  m80(:, 1:16) = false;

    ap.mode = '3D'; ap.num_beats = nB; ap.dt = 1; ap.APD_levels = lv;
    ap.APD_beat = cell(nB, numel(lv));
    for b = 1:nB
        p = (-1)^(b+1);                              % +1 on odd beats
        ap.APD_beat{b,1} = 25 + 3 * p * s50;
        ap.APD_beat{b,2} = 40 + 4 * p * s50;         % 8 ms alternans, 20% of mean
        a80 = 60 + 4 * p * s80;  a80(~m80) = NaN;
        ap.APD_beat{b,3} = a80;
    end
    for k = 1:numel(lv)
        st = cat(3, ap.APD_beat{:, k});
        ap.(sprintf('APD%d_alt', lv(k))) = mean(st(:,:,1:2:end) - st(:,:,2:2:end), 3);
    end
    ap.amp_alt_map     = zeros(R, C);                % full tissue extent
    ap.phase_map       = sign(ap.APD80_alt);
    ap.phase_angle_map = zeros(R, C);
end


function r = stale_ratio(ap)
%STALE_RATIO  The pre-2026-09-28 ratio: APD50 significance, APD80 phase.
    st  = cat(3, ap.APD_beat{:, 2});
    alt = ap.APD50_alt;
    m   = mean(st, 3);
    sig = ~isnan(alt) & m > 0 & abs(alt) > 0.15 * m;
    ph  = ap.phase_map(sig);
    r   = max(nnz(ph > 0), nnz(ph < 0)) / nnz(sig);
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
