function mv_df_interp_test()
%MV_DF_INTERP_TEST  Regression test for sub-bin DF refinement and the floor flag.
%
%   Synthetic paced traces at 1 kHz, 4 s (bin = 0.244 Hz, the geometry of the
%   paced-DF validation), at rates that fall between bins.  Asserts:
%     1. bin-quantised DF is within half a bin of the truth (the right bin wins)
%     2. refined DF is within 0.1 bin of the truth, at every rate tested
%     3. an on-bin rate is left exactly on its bin
%     4. the 3-D path gives the same DF as the 1-D path for a uniform stack
%     5. a trace whose sub-band power outranks the rhythm is flagged
%        peakAtEdge = -1 (the paced-DF harness reports it as floor-pinned)
%     6. RI/OI are unchanged by the refinement (they live on the bin grid)
%
%   Uses a train of AP-like pulses, not a sinusoid, so the spectrum has the
%   harmonic structure real recordings have.

    fs = 1000; T = 4; t = (0:T*fs-1)' / fs;
    binw = fs / 4096;
    rates = [1.333 2.0 3.333 5.0 6.667 10.0];        % Hz; 1.333/3.333/6.667 are off-bin
    fprintf('mv_df_interp_test: fs %d, %g s, bin %.4f Hz\n', fs, T, binw);

    for r = rates
        x = pulse_train(t, r, 0.02);
        [df, ri, oi, sp] = cardiacSpectralMetrics(x, fs);
        [df0, ri0, oi0]  = cardiacSpectralMetrics(x, fs, 'PeakInterp', false);
        eb  = abs(sp.dfBin - r) / binw;
        er  = abs(df - r) / binw;
        check(sprintf('%.3f Hz: bin peak within half a bin', r), eb <= 0.5 + 1e-9, sprintf('%.2f bins', eb));
        check(sprintf('%.3f Hz: refined DF within 0.1 bin', r), er <= 0.1, sprintf('%.3f bins (%.4f Hz)', er, df - r));
        check(sprintf('%.3f Hz: PeakInterp=false is the bin peak', r), abs(df0 - sp.dfBin) < 1e-12, sprintf('%.4f Hz', df0));
        check(sprintf('%.3f Hz: RI/OI unchanged by refinement', r), abs(ri - ri0) < 1e-12 && abs(oi - oi0) < 1e-12, ...
              sprintf('RI %.3f OI %.3f', ri, oi));
    end

    % On-bin rate: 10 bins = 2.4414 Hz exactly.
    r = 10 * binw;
    [df, ~, ~, sp] = cardiacSpectralMetrics(pulse_train(t, r, 0.02), fs);
    check('on-bin rate stays on its bin', abs(df - sp.dfBin) < 0.02 * binw, sprintf('shift %.4f bins', (df - sp.dfBin) / binw));

    % 3-D path: 4 x 4 stack of the same trace with a NaN corner.
    X = repmat(reshape(pulse_train(t, 3.333, 0.02), 1, 1, []), 4, 4, 1);
    X(1, 1, :) = NaN;
    [dfm, ~, ~, spm] = cardiacSpectralMetrics(X, fs);
    [df1, ~, ~, sp1] = cardiacSpectralMetrics(squeeze(X(2, 2, :)), fs);
    check('3-D path matches 1-D path', abs(dfm(2, 2) - df1) < 1e-9 && isnan(dfm(1, 1)) && ...
          abs(spm.dfBin(3, 3) - sp1.dfBin) < 1e-9, sprintf('3-D %.4f vs 1-D %.4f Hz', dfm(2, 2), df1));

    % Floor-pinned: a 0.7 Hz oscillation five times the pulse amplitude.
    x = pulse_train(t, 3.333, 0.02) + 5 * sin(2 * pi * 0.7 * t);
    ws = warning('off', 'cardiacSpectralMetrics:peakAtBandEdge'); c = onCleanup(@() warning(ws));
    [df, ~, ~, sp] = cardiacSpectralMetrics(x, fs);
    check('sub-band power pins the peak to the band floor', sp.autoBand && sp.peakAtEdge == -1 && df < 1, ...
          sprintf('DF %.3f Hz, band [%.2f %.2f]', df, sp.band(1), sp.band(2)));
    [df, ~, ~, sp] = cardiacSpectralMetrics(pulse_train(t, 3.333, 0.02), fs);
    check('clean trace is not flagged', sp.peakAtEdge == 0 && abs(df - 3.333) < 0.1 * binw, sprintf('DF %.4f Hz', df));

    fprintf('All checks passed.\n');
end


function x = pulse_train(t, rate, width_s)
%PULSE_TRAIN  AP-like pulses (fast rise, exponential fall) at a fixed rate.
    x = zeros(size(t));
    for k = 0:floor(t(end) * rate)
        t0 = k / rate + 0.05;
        m  = t >= t0;
        x(m) = x(m) + (1 - exp(-(t(m) - t0) / 0.002)) .* exp(-(t(m) - t0) / (3 * width_s));
    end
end


function check(name, cond, detail)
    if cond
        fprintf('   PASS  %-48s (%s)\n', name, detail);
    else
        error('mv_df_interp_test:failed', 'FAIL  %s (%s)', name, detail);
    end
end
