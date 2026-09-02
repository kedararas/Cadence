function mv_tau_test
%MV_TAU_TEST  Regression test for extract_tau, against analytic truth.
%
%   mv_tau_test
%
%   NOT a study instrument.  It takes no reviewer input and produces no
%   per-recording result: it is a single pass/fail run on synthetic calcium
%   transients whose decay time constant is known by construction.  Run it
%   once, and again after any edit to extract_tau.  It is what puts tau in
%   the "synthetic ground truth" column of the coverage table, the same way
%   mv_synth_cv does for conduction velocity and mv_rise_test does for rise
%   time.
%
%   MIRROR BELOW REPRODUCES THE PROPOSED extract_tau (extract_tau_FIX.m)
%   OPERATING ON A PLAIN [R x C x T] ARRAY.  CURRENT reproduces the version
%   presently inside Cadence_Feature_Extraction.mlapp.  If either is edited,
%   edit this file to match in the same commit, or the test quietly starts
%   validating code that no longer exists.
%
%   THE SIGNAL, and why its tau is analytic
%     y(t) = offset                                    t <  t0
%          = offset + A*(t-t0)/t_rise                  t0 <= t < t0+t_rise
%          = offset + A*exp(-(t-t_pk)/tau)             t >= t_pk = t0+t_rise
%     The peak is EXACTLY A above baseline and lands exactly on a sample, and
%     the decay is exactly mono-exponential from that sample.  So the 1/e
%     crossing sits exactly tau after the peak, and the decay slope over any
%     sub-interval is exactly -1/tau.  Both estimators share one truth, and
%     any departure is the estimator's, not the signal's.
%
%   SIX PROPERTIES ASSERTED
%     1. INDEPENDENT OF ANALYSIS WINDOW LENGTH.  The window is the beat, so a
%        window-dependent tau is a pacing-rate-dependent tau, i.e. a
%        cross-species confound.  The current version normalises with
%        normalize_data's default 1st-99th percentile band, whose bounds are
%        set by how much of the window is diastole, and then compares against
%        ABSOLUTE levels -- so its tau moves with the window even though the
%        transient is identical.
%     2. NO WORSE THAN THE CURRENT VERSION at any SNR, in BOTH slots.
%     3. A MISSING BASELINE TERM IS FATAL when diastole is not at zero.
%        MIRROR_NOBASE is the proposed engine with the baseline subtraction
%        removed and nothing else changed, so the offset sweep attributes the
%        error to that one term.  (normalize_data is affine, so CURRENT is
%        invariant to a pure DC offset; its version of this defect surfaces
%        as the window dependence in test 1, where the percentile floor does
%        not coincide with the true asymptote.)
%     4. A PURE-NOISE PIXEL RETURNS NaN in both slots.  The current version
%        returns a finite tau for a pixel that contains no transient at all.
%     5. INCOMPLETE RELAXATION (the fast-pacing case) STILL YIELDS BOTH
%        SLOTS.  The current version's `continue` fires before the model-free
%        estimator is computed, so a slot-1 condition silently voids slot 6.
%     6. THE SUB-FRAME CROSSING CONVENTION IS findCrossMs's, to 1e-9 ms.
%
%   Requires prctile (Statistics Toolbox), because the shipped code uses it.

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(fullfile(root, 'signal_conditioning_helper'));
    addpath(fullfile(root, 'feature_extraction_helper'));

    rng(20260901, 'twister');

    fs        = 1000;                  % Hz
    tau_true  = 62.5;                  % ms -- deliberately NOT a whole number,
                                       % so the current version's round() to
                                       % integer ms is visible on its own
    A         = 1.0;
    t0        = 20;                    % ms, upstroke start
    t_rise    = 10;                    % ms, linear upstroke (calcium-like)
    tau_start = 0.9;                   % the app's 90% / 20% dropdown setting
    tau_end   = 0.2;
    R         = 20;  C = 20;           % 400 pixels per condition

    fprintf('\nMV_TAU_TEST   fs = %g Hz   truth tau = %.4g ms   fit window %d%%-%d%%\n', ...
            fs, tau_true, round(100*tau_start), round(100*tau_end));
    fprintf('%s\n', repmat('=', 1, 78));

    % ================= 1. window-length sweep ================================
    % Same transient, same noise samples; only the amount of trailing diastole
    % changes.  A correct estimator cannot see the difference.
    wins   = [150 400 1000 2000];
    snr_w  = 20;
    Tmax   = round(max(wins) * fs / 1000);
    noise  = randn(R, C, Tmax) * (A / snr_w);     % ONE realisation, truncated

    fprintf('\n1. WINDOW-LENGTH SWEEP  (SNR %g, identical transient and noise)\n', snr_w);
    header();
    nf = nan(1, numel(wins)); ne = nf; cf = nf; ce = nf;
    for k = 1:numel(wins)
        T  = round(wins(k) * fs / 1000);
        Yk = synth(R, C, T, fs, A, 0, t0, t_rise, tau_true) + noise(:,:,1:T);
        [a, b] = mirror(Yk,  fs, tau_start, tau_end);
        [c, e] = current(Yk, fs, tau_start, tau_end);
        nf(k) = med(a); ne(k) = med(b); cf(k) = med(c); ce(k) = med(e);
        row(sprintf('%4d ms window', wins(k)), a, b, c, e, tau_true);
    end
    spread_new_fit = spr(nf);  spread_cur_fit = spr(cf);
    spread_new_1e  = spr(ne);  spread_cur_1e  = spr(ce);
    fprintf('\n   peak-to-peak across windows, as %% of truth:\n');
    fprintf('     slot 1 (fit)        new %6.2f%%     current %6.2f%%\n', ...
            100*spread_new_fit/tau_true, 100*spread_cur_fit/tau_true);
    fprintf('     slot 6 (1/e)        new %6.2f%%     current %6.2f%%\n', ...
            100*spread_new_1e/tau_true,  100*spread_cur_1e/tau_true);

    % ================= 2. SNR sweep ==========================================
    T_snr  = round(400 * fs / 1000);
    snrs   = [Inf 20 10 5];
    fprintf('\n2. SNR SWEEP  (400 ms window)\n');
    header();
    err_new_fit = nan(size(snrs)); err_cur_fit = err_new_fit;
    err_new_1e  = err_new_fit;     err_cur_1e  = err_new_fit;
    cov_new     = err_new_fit;     cov_cur     = err_new_fit;
    for k = 1:numel(snrs)
        Yk = synth(R, C, T_snr, fs, A, 0, t0, t_rise, tau_true);
        if isfinite(snrs(k))
            Yk = Yk + randn(R, C, T_snr) * (A / snrs(k));
        end
        [a, b] = mirror(Yk,  fs, tau_start, tau_end);
        [c, e] = current(Yk, fs, tau_start, tau_end);
        err_new_fit(k) = med(a) - tau_true;  err_cur_fit(k) = med(c) - tau_true;
        err_new_1e(k)  = med(b) - tau_true;  err_cur_1e(k)  = med(e) - tau_true;
        cov_new(k)     = mean(isfinite(a(:)));
        cov_cur(k)     = mean(isfinite(c(:)));
        if isfinite(snrs(k)), lbl = sprintf('SNR %5.0f', snrs(k));
        else,                 lbl = 'SNR   Inf'; end
        row(lbl, a, b, c, e, tau_true);
        fprintf('   %-16s %9s %9s %9s %9s   coverage new %5.1f%%  cur %5.1f%%\n', ...
                '', '', '', '', '', 100*cov_new(k), 100*cov_cur(k));
    end

    % ================= 3. diastolic offset ===================================
    % The transient is identical; only the level diastole sits at changes.
    offs = [0 0.25 1.0 4.0] * A;
    fprintf('\n3. DIASTOLIC-OFFSET SWEEP  (400 ms window, SNR 20)\n');
    fprintf('   %-16s %10s %10s %10s %10s\n', ...
            'offset (x amp)', 'new fit', 'NO-BASE', 'cur fit', 'new 1/e');
    fprintf('   %s\n', repmat('-', 1, 62));
    off_new = nan(size(offs)); off_nb = off_new; off_cur = off_new;
    noise_o = randn(R, C, T_snr) * (A / 20);      % same realisation throughout
    for k = 1:numel(offs)
        Yk = synth(R, C, T_snr, fs, A, offs(k), t0, t_rise, tau_true) + noise_o;
        [a, b]  = mirror(Yk, fs, tau_start, tau_end);
        nb      = mirror_nobase(Yk, fs, tau_start, tau_end);
        c       = current(Yk, fs, tau_start, tau_end);
        off_new(k) = med(a); off_nb(k) = med(nb); off_cur(k) = med(c);
        fprintf('   %-16.2f %10.2f %10.2f %10.2f %10.2f\n', ...
                offs(k)/A, off_new(k), off_nb(k), off_cur(k), med(b));
    end
    fprintf('   %s\n', repmat('-', 1, 62));
    fprintf('   peak-to-peak: new %.2f ms (%.2f%%),  no-baseline %.1f ms (%.0f%%)\n', ...
            spr(off_new), 100*spr(off_new)/tau_true, ...
            spr(off_nb),  100*spr(off_nb)/tau_true);

    % ================= 4. pure noise =========================================
    Yn = randn(R, C, T_snr) * (A / 10);           % no transient whatsoever
    [an, bn] = mirror(Yn,  fs, tau_start, tau_end);
    [cn, en] = current(Yn, fs, tau_start, tau_end);
    fprintf('\n4. PURE-NOISE PIXELS (no transient present)\n');
    fprintf('   finite tau returned:  new  slot1 %5.1f%%  slot6 %5.1f%%\n', ...
            100*mean(isfinite(an(:))), 100*mean(isfinite(bn(:))));
    fprintf('                         cur  slot1 %5.1f%%  slot6 %5.1f%%  (median %.1f, %.1f ms)\n', ...
            100*mean(isfinite(cn(:))), 100*mean(isfinite(en(:))), med(cn), med(en));

    % ================= 5. incomplete relaxation (fast pacing) ================
    % 100 ms window: the transient never falls to tau_end (0.2), so the
    % current version's `continue` fires -- and takes slot 6 with it, even
    % though the 1/e crossing is comfortably inside the window.
    % The window has to end BETWEEN tau_end (0.2) and the 1/e level (0.368)
    % for this to be the case under test at all, which leaves little room:
    % 115 ms stops at 0.257 of amplitude.  SNR 100 keeps both margins many
    % sigma wide, because what is being demonstrated is a structural defect
    % -- the order of two statements -- not a noise limit.  At SNR 20 the new
    % slot 6 still resolves 67-86% of these pixels; the current version
    % resolves none at any SNR.
    T_fast = 115;
    snr_f  = 100;
    Yf = synth(R, C, T_fast, fs, A, 0.3, t0, t_rise, tau_true) + ...
         randn(R, C, T_fast) * (A / snr_f);
    [af, bf] = mirror(Yf,  fs, tau_start, tau_end);
    [cf2, ef] = current(Yf, fs, tau_start, tau_end);
    fprintf('\n5. INCOMPLETE RELAXATION  (%d ms window, SNR %g, trace stops at %.2f of amp)\n', ...
            T_fast, snr_f, exp(-(T_fast - t0 - t_rise)/tau_true));
    fprintf('   %-16s %9s %9s %9s %9s\n', '', 'new fit', 'new 1/e', 'cur fit', 'cur 1/e');
    fprintf('   %-16s %9.2f %9.2f %9.2f %9.2f\n', 'median tau (ms)', ...
            med(af), med(bf), med(cf2), med(ef));
    fprintf('   %-16s %8.1f%% %8.1f%% %8.1f%% %8.1f%%\n', 'coverage', ...
            100*mean(isfinite(af(:))), 100*mean(isfinite(bf(:))), ...
            100*mean(isfinite(cf2(:))), 100*mean(isfinite(ef(:))));

    % ================= 6. crossing convention ================================
    Y1  = synth(1, 1, T_snr, fs, A, 0.3, t0, t_rise, tau_true);
    [~, got, dbg] = mirror(Y1, fs, tau_start, tau_end);
    ref = findCrossMs((dbg.Yc(dbg.pk:end) - dbg.base) / dbg.pk_amp, 1/exp(1), 1/fs);
    fprintf('\n6. findCrossMs AGREEMENT   helper %.6f ms   vectorised %.6f ms   diff %.2g\n', ...
            ref, got, abs(ref - got));

    % ================= assertions ============================================
    fprintf('\n%s\n', repmat('=', 1, 78));

    % (1) window independence -- the headline property
    assert(spread_new_fit < 0.01 * tau_true, ...
        'slot 1 is window-dependent: %.3f ms spread (%.2f%% of truth)', ...
        spread_new_fit, 100*spread_new_fit/tau_true);
    assert(spread_new_1e < 0.01 * tau_true, ...
        'slot 6 is window-dependent: %.3f ms spread (%.2f%% of truth)', ...
        spread_new_1e, 100*spread_new_1e/tau_true);
    % and the current version must NOT be flat, or test 1 proves nothing
    assert(spread_cur_fit > 5 * spread_new_fit && spread_cur_fit > 0.05 * tau_true, ...
        'current slot 1 unexpectedly window-independent (%.3f ms) -- test is toothless', ...
        spread_cur_fit);

    % (2) accuracy and coverage, at every SNR, in both slots
    for k = 1:numel(snrs)
        assert(abs(err_new_fit(k)) <= abs(err_cur_fit(k)) + 1e-9, ...
            'slot 1 worse than current at SNR %g: %+.3f vs %+.3f ms', ...
            snrs(k), err_new_fit(k), err_cur_fit(k));
        assert(abs(err_new_1e(k)) <= abs(err_cur_1e(k)) + 1e-9, ...
            'slot 6 worse than current at SNR %g: %+.3f vs %+.3f ms', ...
            snrs(k), err_new_1e(k), err_cur_1e(k));
        assert(abs(err_new_fit(k)) < 0.05 * tau_true, ...
            'slot 1 off by %+.3f ms at SNR %g (>5%% of truth)', err_new_fit(k), snrs(k));
    end

    % (3) the baseline term is what buys offset invariance
    assert(spr(off_new) < 0.01 * tau_true, ...
        'slot 1 depends on diastolic offset: %.3f ms spread', spr(off_new));
    assert(spr(off_nb) > 10 * spr(off_new), ...
        'dropping the baseline term did not hurt -- the offset sweep is toothless');

    % (4) no transient, no tau
    assert(~any(isfinite(an(:))), ...
        '%d pure-noise pixels returned a finite fitted tau', sum(isfinite(an(:))));
    assert(~any(isfinite(bn(:))), ...
        '%d pure-noise pixels returned a finite 1/e tau', sum(isfinite(bn(:))));

    % (5) a slot-1 condition must not be able to void slot 6
    assert(mean(isfinite(bf(:))) > 0.95, ...
        'slot 6 lost on incomplete relaxation: only %.1f%% finite', ...
        100*mean(isfinite(bf(:))));
    assert(abs(med(bf) - tau_true) < 0.05 * tau_true, ...
        'slot 6 inaccurate on a truncated window: %.2f vs %.2f ms', med(bf), tau_true);
    assert(mean(isfinite(cf2(:))) < 0.05 && mean(isfinite(ef(:))) < 0.05, ...
        'current version unexpectedly survives truncation -- test 5 is toothless');

    % (6) sub-frame convention preserved
    assert(abs(ref - got) < 1e-9, ...
        'vectorised 1/e crossing disagrees with findCrossMs by %.3g ms', abs(ref - got));

    fprintf('All assertions passed.\n\n');
end


% ======================================================================== %
% Synthetic signal                                                         %
% ======================================================================== %
function Y = synth(R, C, T, fs, A, offset, t0, t_rise, tau)
%SYNTH  [R x C x T] stack of identical analytic transients (see header).
    t  = (0:T-1) / fs * 1000;                 % ms
    tp = t0 + t_rise;                         % peak time, on a sample
    g  = zeros(1, T);
    up = t >= t0 & t < tp;
    dn = t >= tp;
    g(up) = (t(up) - t0) / t_rise;
    g(dn) = exp(-(t(dn) - tp) / tau);
    Y = repmat(reshape(offset + A * g, 1, 1, T), R, C, 1);
end


% ======================================================================== %
% MIRROR of the PROPOSED extract_tau (see extract_tau_FIX.m)               %
% ======================================================================== %
function [tau_fit, tau_1e, dbg] = mirror(win_data, acqFreq, tau_start, tau_end)
    [tau_fit, tau_1e, dbg] = tau_engine(win_data, acqFreq, tau_start, tau_end, true);
end

function tau_fit = mirror_nobase(win_data, acqFreq, tau_start, tau_end)
%MIRROR_NOBASE  Identical, except the decay model omits the baseline term
%   (defect 2 in isolation).  Everything else -- own-baseline landmarks,
%   gates, no rounding -- is unchanged.
    tau_fit = tau_engine(win_data, acqFreq, tau_start, tau_end, false);
end

function [tau_fit, tau_1e, dbg] = tau_engine(win_data, acqFreq, tau_start, tau_end, use_base)
    [num_rows, num_cols, T] = size(win_data);
    dt = 1 / acqFreq;

    % No normalize_data: every level below is referenced to THIS pixel's own
    % baseline and amplitude, so a global rescaling would cancel -- and the
    % default percentile band CLIPS, which does not cancel.
    Y  = double(reshape(win_data, num_rows * num_cols, T));
    N  = size(Y, 1);
    tt = 1:T;
    idx = (1:N)';

    % Running MEDIANS, used for landmarks and for the level crossing only.
    % A running median is exactly unbiased on any MONOTONE stretch -- the
    % median of a monotone window is its centre sample -- so it removes noise
    % from the decay without moving the crossing, which a moving average
    % would not do.  Both estimates themselves are taken from the raw samples.
    ws = max(3, round(0.005 * acqFreq));               % 5 ms, landmarks
    wc = max(5, 2*floor(0.015 * acqFreq / 2) + 1);     % 15 ms, odd, crossing
    Ys = movmedian(Y, ws, 2);
    Yc = movmedian(Y, wc, 2);

    % ---- landmarks: upstroke anchor and diastolic foot ---------------------
    b0  = prctile(Ys,  5, 2);
    p0  = prctile(Ys, 99, 2);
    rg  = p0 - b0;
    mid = b0 + 0.5 * rg;
    low = b0 + 0.1 * rg;

    hi_mid    = Ys >= mid;
    sustained = hi_mid(:, 1:end-1) & hi_mid(:, 2:end);
    [has_anchor, anchor] = max(sustained, [], 2);
    anchor(~has_anchor)  = T;

    % The FOOT, not the mid-amplitude anchor, bounds the baseline sample set.
    % Taking the median of everything before the anchor lets the lower half of
    % the upstroke into the estimate, and how much gets in depends on where
    % the percentile band put `mid` -- i.e. on the window length.  That alone
    % moved tau by 1% across a 150-2000 ms window.
    cand = (Ys <= low) & (tt < anchor);
    [has_foot, j_rev] = max(fliplr(cand), [], 2);
    foot = T - j_rev + 1;
    foot(~has_foot) = max(anchor(~has_foot) - 1, 1);

    pre     = Y;  pre(tt > foot) = NaN;
    base    = median(pre, 2, 'omitnan');
    n_pre   = sum(~isnan(pre), 2);
    min_pre = max(3, round(0.005 * acqFreq));
    base(n_pre < min_pre) = NaN;

    % ---- noise, from first differences over the WHOLE trace ---------------
    % A short diastolic segment gives a MAD that is itself noisy, and an
    % underestimate there opens the amplitude gate.  Successive differences
    % are dominated by noise wherever the signal is smooth, and the median is
    % robust to the upstroke.
    sigma = 1.4826 * median(abs(diff(Y, 1, 2)), 2, 'omitnan') / sqrt(2);

    % ---- peak on the filtered trace, height from the raw sample -----------
    post         = Ys;  post(tt < anchor) = NaN;
    [pk_val, pk] = max(post, [], 2);
    amp          = pk_val - base;                       % robust, for gating
    pk_amp       = Y(sub2ind([N T], idx, pk)) - base;   % unbiased height

    % ---- amplitude / SNR gate ---------------------------------------------
    % sqrt(2*log n) * sigma is the expected largest excursion of n zero-mean
    % samples; the running median leaves ~T/ws independent ones, each with
    % sd ~ sigma*sqrt(pi/(2*ws)).  Three further sigma of margin puts the gate
    % clear of the noise maximum's own spread, so a pixel with no transient
    % cannot clear it, while a transient at SNR 5 still can.
    sigma_f  = sigma * sqrt(pi / (2 * ws));
    amp_gate = (sqrt(2 * log(max(T / ws, 3))) + 3) * sigma_f;
    ok_amp   = isfinite(base) & isfinite(sigma) & (sigma > 0) & ...
               (amp > 0) & (amp > amp_gate) & (pk_amp > 0);

    amp_safe = amp;     amp_safe(~ok_amp) = NaN;
    pk_safe  = pk_amp;  pk_safe(~ok_amp)  = NaN;
    Yn       = (Y - base) ./ amp_safe;        % own-amplitude units

    % ======== slot 6: model-free 1/e, computed FIRST and independently =====
    % A slot-1 failure must not be able to discard this, and vice versa; at
    % fast pacing the transient may never relax to tau_end at all.
    lev  = 1 / exp(1);
    labs = base + lev .* pk_safe;
    [hit1, ix] = max(cumsum((Yc <= labs) & (tt > pk), 2) == 1, [], 2);
    ix   = max(ix, 2);
    y1   = Yc(sub2ind([N T], idx, ix - 1));
    y2   = Yc(sub2ind([N T], idx, ix));
    den1 = y1 - y2;
    frac = (y1 - labs) ./ den1;               % findCrossMs, exactly
    frac(~(den1 > 0)) = 0;                    % y1 == y2 -> frac 0
    frac = min(max(frac, 0), 1);
    t1e  = ((ix - 1 + frac) - pk) * dt * 1000;
    tau_1e = nan(N, 1);
    g1e    = hit1 & ok_amp & isfinite(t1e) & (t1e > 0);
    tau_1e(g1e) = t1e(g1e);
    tau_1e = reshape(tau_1e, num_rows, num_cols);   % ms, NOT rounded

    % ======== slot 1: decay fit with a baseline term =======================
    if tau_start >= 1
        rep_start = pk;
        has_start = ok_amp;
    else
        [has_start, rep_start] = max(cumsum((Yn <= tau_start) & (tt >= pk), 2) == 1, [], 2);
        rep_start(~has_start)  = T;
    end
    % If the trace never reaches tau_end (fast pacing, incomplete relaxation)
    % fit the decay that IS present rather than discarding the pixel -- but
    % only if it still covers a halving, so the estimate comes from real decay
    % and not from noise on a plateau.
    [has_end, rep_end] = max(cumsum((Yn <= tau_end) & (tt >= rep_start), 2) == 1, [], 2);
    rep_end(~has_end)  = T;
    rep_end = max(rep_end, rep_start);

    v_s = Yn(sub2ind([N T], idx, rep_start));
    v_e = Yn(sub2ind([N T], idx, rep_end));
    ok_span = (v_s > 0) & (v_e <= 0.5 * v_s);

    % Amplitude above the pixel's own diastolic level.  Fitting the raw trace
    % instead assumes the decay asymptotes to ZERO; a positive residual
    % asymptote flattens the decay and overestimates tau without limit.
    if use_base
        D = Y - base;
    else
        D = Y;
    end
    inwin = (tt >= rep_start) & (tt <= rep_end);
    Trel  = ((tt - 1) - (rep_start - 1)) * dt;      % s, per-pixel origin

    n_win = sum(inwin, 2);
    [hf, i_first] = max(inwin, [], 2);
    [~,  j2]      = max(fliplr(inwin), [], 2);
    span_ms = (T - j2 + 1 - i_first) * dt * 1000;
    min_pts = max(4, round(0.005 * acqFreq));
    min_span_ms = 5;

    % ---- seed: weighted log-linear regression, vectorised over pixels -----
    valid = inwin & (D > 0) & isfinite(D);
    Z = zeros(N, T);  Z(valid) = log(D(valid));
    w = zeros(N, T);  w(valid) = D(valid);
    S0  = sum(w, 2);
    S1  = sum(w .* Trel, 2);
    S2  = sum(w .* Trel.^2, 2);
    Sy  = sum(w .* Z, 2);
    Sty = sum(w .* Trel .* Z, 2);
    den = S0 .* S2 - S1.^2;
    den(~(den > 0)) = NaN;
    tau_s = -1 ./ ((S0 .* Sty - S1 .* Sy) ./ den);       % seconds

    % ---- refinement: least squares in the domain the noise lives in -------
    % The log transform inflates the variance of the small tail samples and,
    % worse, censors the ones that go negative -- both bias the slope, by
    % ~+15% at SNR 5.  Gauss-Newton (variable projection: the amplitude is
    % linear, so it is profiled out and the residual is orthogonal to it)
    % removes that; the log fit is only the starting point.  Samples at or
    % below baseline are kept here, which is what makes it unbiased.
    Dm = D;     Dm(~inwin) = 0;
    Tm = Trel;  Tm(~inwin) = 0;
    tau_n = tau_s;
    tau_n(~isfinite(tau_n) | tau_n <= 0) = 20 * dt;      % fallback seed
    tau_n = min(max(tau_n, dt), T * dt);
    for it = 1:8
        E    = exp(-Tm ./ tau_n) .* inwin;
        Aa   = sum(Dm .* E, 2) ./ max(sum(E.^2, 2), eps);
        Rr   = (Dm - Aa .* E) .* inwin;
        G    = Aa .* E .* Tm ./ tau_n.^2;
        step = sum(Rr .* G, 2) ./ max(sum(G.^2, 2), eps);
        step = min(max(step, -0.5 * tau_n), 0.5 * tau_n);   % damped
        tau_n = min(max(tau_n + step, dt), T * dt);
    end

    % A failed fit returns NaN.  No floor, no fallback, and in particular no
    % silent substitution of the model-free estimator into this slot.
    ok_fit = ok_amp & ok_span & has_start & hf & (n_win >= min_pts) & ...
             (span_ms >= min_span_ms) & isfinite(tau_n) & ...
             (tau_n > dt) & (tau_n < T * dt);
    tau_fit = nan(N, 1);
    tau_fit(ok_fit) = tau_n(ok_fit) * 1000;           % ms, NOT rounded
    tau_fit = reshape(tau_fit, num_rows, num_cols);

    dbg = struct('Yc', Yc, 'base', base, 'pk', pk, 'pk_amp', pk_amp, ...
                 'amp', amp, 'sigma', sigma, 'rep_start', rep_start, ...
                 'rep_end', rep_end, 'ok_fit', ok_fit);
end


% ======================================================================== %
% The CURRENT implementation, kept as the baseline the assertions beat.    %
% Transcribed from Cadence_Feature_Extraction.mlapp extract_tau(app, d).   %
% ======================================================================== %
function [tau_constant, tau_mode_constant] = current(win_data, acqFreq, tau_start, tau_end)
    dt = 1 / acqFreq;
    [num_rows, num_cols, ~] = size(win_data);
    tau_constant      = nan(num_rows, num_cols);
    tau_mode_constant = nan(num_rows, num_cols);

    for x = 1:num_rows
        for y = 1:num_cols
            win_oap = squeeze(normalize_data(win_data(x,y,:)));   % CLIPS to [0,1]
            [~, max_loc] = max(win_oap);
            if tau_start == 1.0
                rep_start = max_loc;
            else
                rep_start = max_loc + find(win_oap(max_loc:end) <= tau_start, 1) - 1;
            end
            if isempty(rep_start), continue; end

            rep_end = find(win_oap(rep_start:end) <= tau_end, 1);
            if isempty(rep_end)              % bad data -- also voids slot 6
                continue;
            else
                rep_end = rep_start + rep_end - 1;
            end

            decay_oap  = win_oap(rep_start:rep_end);
            decay_time = (0:numel(decay_oap)-1).' * dt;
            post       = win_oap(max_loc:end);
            decay_63   = findCrossMs(post, 1/exp(1), dt);

            Y = decay_oap(:);
            Y(Y <= 0) = NaN;
            m = isfinite(Y);
            tau_seed_s = NaN;
            if nnz(m) >= 3
                pLL = polyfit(decay_time(m), log(Y(m)), 1);   % no baseline term
                tau_seed_s = -1 / pLL(1);
            end
            if ~isfinite(tau_seed_s) || tau_seed_s <= 0
                tau_seed_s = max(2*dt, decay_63/1000);        % max(x,NaN) -> x
            end
            tau_seed_ms = tau_seed_s * 1000;

            if tau_seed_s > 0
                tau_constant(x,y) = round(tau_seed_ms);       % 1 ms quantisation
            end
            if isfinite(decay_63) && decay_63 > 0
                tau_mode_constant(x,y) = round(decay_63);
            end
        end
    end
end


% ======================================================================== %
% Reporting helpers                                                        %
% ======================================================================== %
function header()
    fprintf('   %-16s %9s %9s %9s %9s   %9s %9s\n', ...
            'condition', 'new fit', 'new 1/e', 'cur fit', 'cur 1/e', ...
            'new err%', 'cur err%');
    fprintf('   %s\n', repmat('-', 1, 78));
end

function row(lbl, a, b, c, e, truth)
    fprintf('   %-16s %9.2f %9.2f %9.2f %9.2f   %+8.2f%% %+8.2f%%\n', ...
            lbl, med(a), med(b), med(c), med(e), ...
            100*(med(a)-truth)/truth, 100*(med(c)-truth)/truth);
end

function v = med(x)
    v = median(x(:), 'omitnan');
    if isempty(v), v = NaN; end
end

function s = spr(v)
    v = v(isfinite(v));
    if numel(v) < 2, s = NaN; else, s = max(v) - min(v); end
end
