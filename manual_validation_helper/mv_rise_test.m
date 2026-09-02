function mv_rise_test
%MV_RISE_TEST  Regression test for extract_rise_time, against analytic truth.
%
%   mv_rise_test
%
%   NOT a study instrument.  It takes no reviewer input and produces no
%   per-recording result: it is a single pass/fail run on synthetic signals
%   whose true rise time is known by construction.  Run it once, and again
%   after any edit to extract_rise_time.  It is what puts rise time in the
%   "synthetic ground truth" column of the coverage table, the same way
%   mv_synth_cv does for conduction velocity.
%
%   SHIPPED_METHOD BELOW MIRRORS Cadence_Feature_Extraction.mlapp's
%   extract_rise_time LINE FOR LINE.  If that function is edited, edit this one
%   to match in the same commit, or the test quietly starts validating code
%   that no longer exists.
%
%   Signals, both with analytic truth:
%     linear ramp   time from rise_min to rise_max is exactly
%                   (rise_max - rise_min) * D, wherever the ramp starts
%     logistic      t(L) = t0 - k*log(1/L - 1), so 20-80% is 2*log(4)*k
%
%   Each pixel gets a different SUB-FRAME start offset.  That is what exposes
%   frame quantization: with the upstroke aligned to the sample grid the old
%   integer-frame code happened to be right, which is why the bug survived.
%
%   THREE PROPERTIES ASSERTED
%     1. on a linear ramp, accurate to the APEX-SAMPLING FLOOR.  Not exact:
%        amp is estimated as max(sample) - base, and unless the apex lands on
%        a sample that under-reads the true peak by up to ~1.2%, which pulls
%        the rise_min/rise_max thresholds inward and shortens the result by
%        about 1%.  A floor of the discrete sampling, not a defect, and it is
%        constant rather than window- or rate-dependent.
%     2. independent of ANALYSIS WINDOW LENGTH -- the window is the beat, so a
%        window-dependent result is a pacing-rate-dependent result, i.e. a
%        cross-species confound
%     3. no worse than the superseded version at any SNR, compared ON THE
%        PIXELS BOTH RESOLVE.  That pairing matters: the superseded code
%        silently drops the pixels it cannot handle (6 of 40 survive at SNR 5),
%        so its raw mean is taken over the easy ones and flatters itself.
%        Coverage is asserted separately, and is the more important of the two.
%
%   Requires prctile (Statistics Toolbox), because the shipped code uses it.
%   The marking harness proper avoids that dependency; this mirror cannot.

    fs       = 1000;
    D        = 5;                       % linear upstroke duration, ms
    rise_min = 0.2;
    rise_max = 0.8;
    truth    = (rise_max - rise_min) * D;
    offs     = linspace(0, 0.99, 40);   % sub-frame start offsets, ms
    R        = numel(offs);

    % ================= 1. linear ramp: exactness =============================
    fprintf('Linear ramp, 1 kHz, truth %.3f ms\n', truth);
    Y = ramp_stack(R, 60, offs, D, fs);
    old = original_method(Y, rise_min, rise_max, fs);
    new = shipped_method(Y,  rise_min, rise_max, fs);
    report('  superseded (integer frames)', old, truth);
    report('  shipped', new, truth);

    fprintf('\nAt 2 kHz (the superseded error is one FRAME PERIOD, so it only halves):\n');
    Y2 = ramp_stack(R, 120, offs, D, 2000);
    report('  superseded (integer frames)', original_method(Y2, rise_min, rise_max, 2000), truth);
    report('  shipped', shipped_method(Y2, rise_min, rise_max, 2000), truth);

    % Bound, not exactness: see the apex-sampling note in the header.  The
    % ~1% floor is expected; anything approaching a frame period is not.
    assert(max(abs(new(:) - truth)) < 0.10, ...
           'shipped method exceeds the apex-sampling floor on a linear ramp (max err %.4g ms)', ...
           max(abs(new(:) - truth)));

    % ================= 2. window independence ================================
    % The regression guard that matters most.  The analysis window is the beat,
    % so if this drifts, the metric has become a function of the pacing rate and
    % every cross-species comparison is confounded.
    fprintf('\nWindow-length sweep (same AP, longer window = slower pacing CL):\n');
    fprintf('  window(ms)    superseded        shipped\n');
    wins = [150 200 400 1000];
    ship_means = nan(size(wins));
    for w = 1:numel(wins)
        Yw = ramp_stack(R, wins(w), offs, D, fs);
        a  = original_method(Yw, rise_min, rise_max, fs);
        b  = shipped_method(Yw,  rise_min, rise_max, fs);
        ship_means(w) = mean(b(:), 'omitnan');
        fprintf('  %8d     %+7.1f%%         %+7.1f%%\n', wins(w), ...
                100*(mean(a(:),'omitnan')-truth)/truth, ...
                100*(ship_means(w)-truth)/truth);
    end
    spread = max(ship_means) - min(ship_means);
    assert(spread < 0.01, ...
           ['shipped rise time now depends on window length (spread %.4f ms across ' ...
            '%d-%d ms). That makes it depend on pacing rate.'], ...
           spread, min(wins), max(wins));
    fprintf('  shipped spread across windows: %.4f ms (window-independent)\n', spread);

    % ================= 3. behaviour under noise ==============================
    % Linear interpolation is exact on a ramp by construction, so that alone
    % flatters any candidate.  A logistic upstroke curves between samples and
    % leaves a real residual, and noise is where the landmark choice shows.
    k      = 1.2;
    truth2 = 2 * log(4) * k;
    fprintf('\nLogistic upstroke, 400 ms window, truth %.3f ms:\n', truth2);
    fprintf('  SNR      superseded            shipped\n');
    for snr = [Inf 20 10 5]
        Y3 = logistic_stack(R, 400, offs, k, fs, snr);
        a  = original_method(Y3, rise_min, rise_max, fs);
        b  = shipped_method(Y3,  rise_min, rise_max, fs);
        both = isfinite(a(:)) & isfinite(b(:));
        ea = mean(a(both) - truth2);
        eb = mean(b(both) - truth2);
        fprintf('  %-5s   %+6.2f ms (n=%2d)     %+6.2f ms (n=%2d)   | paired: %+6.2f vs %+6.2f (n=%2d)\n', ...
                num2str(snr), ...
                mean(a(:)-truth2,'omitnan'), sum(isfinite(a(:))), ...
                mean(b(:)-truth2,'omitnan'), sum(isfinite(b(:))), ...
                ea, eb, sum(both));

        % Coverage first: losing pixels is a worse failure than a small bias,
        % because dropped pixels leave holes in the map and bias the survivors.
        assert(sum(isfinite(b(:))) >= sum(isfinite(a(:))), ...
               'shipped method resolves FEWER pixels than the superseded one at SNR %s (%d vs %d)', ...
               num2str(snr), sum(isfinite(b(:))), sum(isfinite(a(:))));

        assert(abs(eb) <= abs(ea) + 1e-9, ...
               'shipped method is worse on the pixels both resolve at SNR %s (%+.3f vs %+.3f)', ...
               num2str(snr), eb, ea);
    end

    fprintf('\nAll assertions passed.\n');
    fprintf('(The residual ~1%% is the apex-sampling floor described in the header.)\n');
end


% ======================================================================
% Signals
% ======================================================================
function Y = ramp_stack(R, T, offs, D, fs)
    t = (0:T-1) / fs * 1000;
    Y = zeros(R, 1, T);
    for i = 1:R
        Y(i,1,:) = ramp_ap(t, 8 + offs(i), D);
    end
end

function Y = logistic_stack(R, T, offs, k, fs, snr)
    t = (0:T-1) / fs * 1000;
    rng(3, 'twister');
    Y = zeros(R, 1, T);
    for i = 1:R
        y = 1 ./ (1 + exp(-(t - (12 + offs(i))) / k)) .* (t < 30) + ...
            exp(-(t - 30) / 40) .* (t >= 30);
        if isfinite(snr); y = y + randn(size(y)) / snr; end
        Y(i,1,:) = y;
    end
end

function y = ramp_ap(t, ta, D)
    y = zeros(size(t));
    up = t >= ta & t < ta + D;
    y(up) = (t(up) - ta) / D;
    dn = t >= ta + D;
    y(dn) = exp(-(t(dn) - ta - D) / 40);
end

function report(name, v, truth)
    e = v(:) - truth;
    fprintf('%-30s  mean %+.3f ms   max |err| %.3f ms\n', ...
            name, mean(e, 'omitnan'), max(abs(e)));
end


% ======================================================================
% MIRROR of the SHIPPED extract_rise_time.  Keep in step with the .mlapp.
% ======================================================================
function rise_time = shipped_method(win_data, rise_min, rise_max, fs)
    [num_rows, num_cols, T] = size(win_data);

    Y  = double(reshape(win_data, num_rows * num_cols, T));
    N  = size(Y, 1);
    tt = 1:T;

    b0  = prctile(Y,  5, 2);
    p0  = prctile(Y, 99, 2);
    mid = b0 + 0.5 * (p0 - b0);

    hi_mid    = Y >= mid;
    sustained = hi_mid(:, 1:end-1) & hi_mid(:, 2:end);
    [has_anchor, anchor] = max(sustained, [], 2);
    anchor(~has_anchor)  = T;

    pre  = Y;  pre(tt >= anchor)  = NaN;
    base = median(pre, 2, 'omitnan');
    post = Y;  post(tt <  anchor) = NaN;
    pk   = max(post, [], 2);
    amp  = pk - base;

    lo_thr = base + rise_min * amp;
    hi_thr = base + rise_max * amp;

    [has_end,   i_end] = max((Y >= hi_thr) & (tt >= anchor), [], 2);
    [has_start, j_rev] = max(fliplr((Y <= lo_thr) & (tt <= anchor)), [], 2);
    i_start = T - j_rev + 1;

    idx = (1:N)';

    ya    = Y(sub2ind([N T], idx, max(i_end - 1, 1)));
    yb    = Y(sub2ind([N T], idx, i_end));
    den   = yb - ya;
    f_end = (hi_thr - ya) ./ den;
    f_end(den <= 0) = 0;
    t_end = (i_end - 1) + min(max(f_end, 0), 1);

    yc      = Y(sub2ind([N T], idx, i_start));
    yd      = Y(sub2ind([N T], idx, min(i_start + 1, T)));
    den2    = yd - yc;
    f_start = (lo_thr - yc) ./ den2;
    f_start(den2 <= 0) = 0;
    t_start = i_start + min(max(f_start, 0), 1);

    rt = t_end - t_start;
    ok = has_anchor & has_end & has_start & isfinite(base) & ...
         (amp > 0) & (i_end > 1) & (i_start < i_end) & (rt > 0);

    rise_frames     = nan(N, 1);
    rise_frames(ok) = rt(ok);
    rise_time       = reshape(rise_frames, num_rows, num_cols) * (1000 / fs);
end


% ======================================================================
% The SUPERSEDED implementation, kept as the baseline the assertions beat.
% Integer frames with outward-opening brackets, and max dV/dt as the landmark.
% ======================================================================
function rise_time = original_method(win_norm, rise_min, rise_max, fs)
    [~, ~, T] = size(win_norm);
    [~, act_loc] = max(diff(win_norm, 1, 3), [], 3);
    t_idx = reshape(1:T, 1, 1, T);

    [hit90, rise_end] = max(cumsum((win_norm >= rise_max) & (t_idx >= act_loc), 3) == 1, [], 3);

    before_act = t_idx <= act_loc;
    flip_valid = flip((win_norm <= rise_min) & before_act, 3);
    [hit20, rs_rev] = max(cumsum(flip_valid, 3) == 1, [], 3);
    rise_start = T - rs_rev + 1;

    rise_time = double(rise_end - rise_start);
    rise_time(~hit90 | ~hit20 | rise_time < 0) = NaN;
    rise_time = rise_time * (1000 / fs);
end
