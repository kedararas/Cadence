function mv_alternans_test
%MV_ALTERNANS_TEST  Regression test for the 3-D alternans crossings.
%
%   mv_alternans_test
%
%   On 2026-09-01 both 3-D paths were changed to INTERPOLATE the
%   repolarization / decay crossing instead of taking the bracketing frame:
%
%     analyzeAPAlternans.m         APD at each level
%     analyzeCaTransientAlternans.m  CaTD (D50/D80) at each level
%
%   Before the change the crossing was a whole frame while the activation
%   anchor was already sub-frame, which quantised APD and made the 3-D maps
%   disagree with the 1-D readouts (which have always used findCrossMs).  For
%   alternans — a difference of two such values — the bias cancels but the
%   per-pixel scatter does not: 0.62 ms SD against 0.25 ms interpolated, which
%   is most of a 1 ms alternans signal.
%
%   This checks the 3-D APD against ANALYTIC truth, and checks that the
%   alternans magnitude is recovered.  Not a study instrument.

    fs = 1000; dt = 1000/fs;
    R  = 24; C = 24;
    % CL must allow FULL relaxation, or consecutive beats overlap and the
    % per-window min/max baseline shifts — an artefact of the stimulus design,
    % not of the crossing code.  tau ~= (APD80-2)/log(5) ~= 73 ms here, so 5*tau
    % plus the onset offset needs ~400 ms.
    cl = 600; nBeats = 8; T = cl*nBeats;

    % Sub-frame variation across pixels is what exposes quantisation: with the
    % upstroke on the sample grid the old code happened to be right.
    % (grid only needed for its size; offsets are built directly below)
    offs    = reshape(linspace(0, 0.99, R*C), R, C);

    apd_long = 120;                 % ms, APD80 of the long beat
    alt_true = 4;                   % ms, long-minus-short
    t = (0:T-1) * dt;
    X = zeros(R, C, T);
    for b = 1:nBeats
        this_apd = apd_long - (mod(b,2) == 0) * alt_true;
        for r = 1:R
            for c = 1:C
                ta = (b-1)*cl + 10 + offs(r,c);
                X(r,c,:) = squeeze(X(r,c,:))' + ap_wave(t, ta, this_apd);
            end
        end
    end

    mask = ones(R, C);
    bf = zeros(nBeats,2);
    for b = 1:nBeats
        bf(b,1) = (b-1)*cl/dt + 1;
        bf(b,2) = min(T, b*cl/dt);
    end

    fprintf('MV_ALTERNANS_TEST  truth: APD80 = %.1f ms long, alternans = %.1f ms\n', ...
            apd_long, alt_true);

    a = analyzeAPAlternans(1:T, X, 'BeatFrames', bf, 'Mask', mask, ...
                           'APD_levels', 80, 'PlotResults', false);

    % ---- per-beat APD against analytic truth ----------------------------
    apd1 = a.APD_beat{1,1};        % beat 1 = long
    apd2 = a.APD_beat{2,1};        % beat 2 = short
    e1 = apd1(:) - apd_long;
    e2 = apd2(:) - (apd_long - alt_true);
    fprintf('\n  beat 1 (long)  : mean %+.3f ms, SD %.3f ms\n', mean(e1,'omitnan'), std(e1,'omitnan'));
    fprintf('  beat 2 (short) : mean %+.3f ms, SD %.3f ms\n', mean(e2,'omitnan'), std(e2,'omitnan'));

    check('per-beat APD is not frame-quantised', std(e1,'omitnan') < 0.30, ...
          sprintf('SD %.3f ms (whole-frame would be ~0.29 minimum and biased)', std(e1,'omitnan')));
    % The residual bias is the APEX-SAMPLING floor, not quantisation: amp is
    % hi-lo with hi the max SAMPLE, which under-reads a between-sample apex by
    % up to ~1%, shifting the threshold and hence the crossing.  Same floor
    % found in extract_rise_time; ~1 ms on a 120 ms APD.
    check('per-beat APD is accurate',            abs(mean(e1,'omitnan')) < 1.5, ...
          sprintf('bias %+.3f ms (apex-sampling floor)', mean(e1,'omitnan')));
    check('short beat is equally accurate',      abs(mean(e2,'omitnan')) < 1.5, ...
          sprintf('bias %+.3f ms', mean(e2,'omitnan')));

    % A whole-frame crossing can only ever land on integer multiples of dt
    % relative to a fixed anchor; interpolated values must not.
    frac = mod(apd1(isfinite(apd1)), dt);
    check('APD takes sub-frame values', any(frac > 1e-6 & frac < dt - 1e-6), ...
          sprintf('%.0f%% of pixels off the frame grid', ...
                  100*mean(frac > 1e-6 & frac < dt - 1e-6)));

    % ---- alternans magnitude --------------------------------------------
    alt = a.APD80_alt;
    ea  = abs(alt(:)) - alt_true;
    fprintf('  alternans      : mean %+.3f ms, SD %.3f ms (truth %.1f)\n', ...
            mean(ea,'omitnan'), std(ea,'omitnan'), alt_true);
    check('alternans magnitude recovered', abs(mean(ea,'omitnan')) < 0.50, ...
          sprintf('bias %+.3f ms', mean(ea,'omitnan')));
    check('alternans scatter is small',    std(ea,'omitnan') < 0.40, ...
          sprintf('SD %.3f ms (was ~0.62 whole-frame)', std(ea,'omitnan')));

    % ---- units: the analysis must honour a real millisecond axis ---------
    % Both wrappers currently pass 1:N (frame indices), which pins dt at 1 and
    % is right only at 1 kHz.  alternans_time_units_FIX.m changes them to pass
    % real ms.  This proves the analysis handles that correctly, by running the
    % identical AP at 2 kHz: with a real ms axis the APD must be unchanged,
    % whereas a 1:N axis would double it.
    fs2 = 2000; dt2 = 1000/fs2;
    T2  = round(cl*nBeats/dt2);
    t2  = (0:T2-1) * dt2;
    X2  = zeros(4, 4, T2);
    for b = 1:nBeats
        this_apd = apd_long - (mod(b,2) == 0) * alt_true;
        for r = 1:4
            for c = 1:4
                ta = (b-1)*cl + 10 + offs(r,c);
                X2(r,c,:) = squeeze(X2(r,c,:))' + ap_wave(t2, ta, this_apd);
            end
        end
    end
    bf2 = zeros(nBeats,2);
    for b = 1:nBeats
        bf2(b,1) = round((b-1)*cl/dt2) + 1;
        bf2(b,2) = min(T2, round(b*cl/dt2));
    end

    a_ms  = analyzeAPAlternans(t2,     X2, 'BeatFrames', bf2, 'Mask', ones(4,4), ...
                               'APD_levels', 80, 'PlotResults', false);
    a_idx = analyzeAPAlternans(1:T2,   X2, 'BeatFrames', bf2, 'Mask', ones(4,4), ...
                               'APD_levels', 80, 'PlotResults', false);
    m_ms  = median(a_ms.APD_beat{1,1}(:),  'omitnan');
    m_idx = median(a_idx.APD_beat{1,1}(:), 'omitnan');
    fprintf('\n  2 kHz, real ms axis   : APD80 = %.2f ms (truth %.1f)\n', m_ms, apd_long);
    fprintf('  2 kHz, 1:N frame axis : APD80 = %.2f ms  <- what ships today\n', m_idx);
    check('real ms axis gives the right APD at 2 kHz', abs(m_ms - apd_long) < 1.5, ...
          sprintf('%.2f vs %.1f', m_ms, apd_long));
    check('frame axis is wrong by the rate ratio',      abs(m_idx - 2*apd_long) < 3.0, ...
          sprintf('%.2f, i.e. 2x too large', m_idx));

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
        fprintf('   PASS  %-38s (%s)\n', name, detail);
    else
        error('mv_alternans_test:failed', 'FAIL  %s (%s)', name, detail);
    end
end
