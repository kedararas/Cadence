function recs = validate_conditioned(d, recording)
%VALIDATE_CONDITIONED  QC gate for the signal-conditioning stage output.
%
%   recs = validate_conditioned(cmos_all_data, recording)
%
%   Run this in the batch loop right after conditioning produces cmos_all_data,
%   before it is saved.  Returns an array of qc_check records.  Each check
%   targets a failure mode that would otherwise pass silently in a batch run
%   (no human watching any single recording).
%
%   Checks
%     acqFreq            present, finite, in a physiological frame-rate band.
%     <CAM>_data         the conditioned stack has finite, non-flat content.
%     <CAM>_SNR          median in-mask SNR above a usable floor.
%     <CAM>_capture      tissue dominant frequency matches the pacing rate
%                        (1:1 capture); WARNs when the tissue runs at its own
%                        rate (spontaneous / arrhythmic / block), which makes
%                        the ensemble average and any paced-APD metric invalid.
%     <CAM>_average      ensemble average has dynamic range (catches the
%                        all-zeros / no-valid-beat ensemble-average failure)
%                        and relaxes within the cycle (decay-fraction).
%
%   See qc_check for the record shape and status semantics.

    recs = repmat(qc_check("","","" ,"PASS",NaN,""), 0, 1);   % empty typed array
    name = string(recording);

    % ---- nested append helper (shares recs/name) ----
    function add(check, status, value, message)
        recs(end+1, 1) = qc_check(name, "conditioning", check, status, value, message);
    end

    % ---- acquisition frame rate ----
    fs = NaN;   % kept in scope for the per-camera capture check below
    if isfield(d, 'acqFreq') && isscalar(d.acqFreq) && isfinite(d.acqFreq)
        fs = double(d.acqFreq);
        if fs >= 50 && fs <= 5000
            add("acqFreq", "PASS", fs, sprintf("%.0f Hz", fs));
        else
            add("acqFreq", "WARN", fs, sprintf("%.0f Hz outside [50,5000] Hz", fs));
        end
    else
        add("acqFreq", "FAIL", NaN, "acqFreq missing or non-finite");
    end

    % ---- per-camera checks ----
    cams = ["CAM1","CAM2","CAM3","CAM4"];
    for c = 1:numel(cams)
        cam = cams(c);
        if ~isfield(d, cam) || isempty(d.(cam))
            continue;   % camera not present in this recording
        end

        X      = d.(cam);
        finite = isfinite(X);
        if ~any(finite(:))
            add(cam + "_data", "FAIL", 0, "all non-finite");
            continue;
        end

        % Dynamic range of the conditioned signal — a flat stack means
        % conditioning collapsed the signal.
        dr = range(X(finite));
        if dr < 1e-6
            add(cam + "_data", "FAIL", dr, "conditioned stack is flat (no dynamic range)");
        else
            add(cam + "_data", "PASS", dr, "");
        end

        % SNR mask: median over tissue (SNR > 0) pixels.
        snrf = cam + "_SNR";
        if isfield(d, snrf) && ~isempty(d.(snrf))
            s    = d.(snrf);
            msnr = median(s(s > 0), 'omitnan');
            if isnan(msnr)
                add(cam + "_SNR", "WARN", NaN, "SNR mask empty / all zero");
            elseif msnr < 2
                add(cam + "_SNR", "WARN", msnr, sprintf("median SNR %.1f below floor 2.0", msnr));
            else
                add(cam + "_SNR", "PASS", msnr, sprintf("median SNR %.1f", msnr));
            end
        end

        % Capture: is the tissue following the pacing 1:1?  The ensemble average
        % and any paced-APD metric are only valid when each stimulus captures one
        % beat.  When the tissue runs at its own rate -- spontaneous activity, an
        % organized tachyarrhythmia, or 2:1 block -- the average keys to a
        % stimulus the tissue is not following and its APD is meaningless.
        % Nothing else here catches that: the stack is finite, has range, and an
        % organized arrhythmia even relaxes within its own (faster) cycle, so it
        % passes the decay-fraction check too.
        if isfield(d, 'analog1') && ~isempty(d.analog1) && isfinite(fs)
            [cr, fp, ftis] = capture_ratio(X, d.analog1, fs);
            if ~isnan(cr)
                if cr < 0.85 || cr > 1.15
                    add(cam + "_capture", "WARN", cr, sprintf( ...
                        "not 1:1 captured (tissue %.1f Hz vs pacing %.1f Hz, ratio %.2f) — do not use the ensemble average or paced APD; route to arrhythmia analysis", ...
                        ftis, fp, cr));
                else
                    add(cam + "_capture", "PASS", cr, sprintf("1:1 capture, tissue %.1f Hz", ftis));
                end
            end
        end

        % Ensemble average: must have dynamic range.  All-zeros means no beat
        % window fell inside the recording (e.g. too few beats / late pacing).
        avgf = cam + "_average";
        if isfield(d, avgf) && ~isempty(d.(avgf))
            % WARN (not FAIL): a degenerate ensemble average is not a conditioning
            % failure. CAM<n>_average is an optional convenience product that
            % alternans/arrhythmia analysis never use (those read the full CAM<n>
            % stack). A flat/zeros average just means no coherent beat window was
            % found -- expected and legitimate for unpaced or arrhythmic
            % (VF/AF/sustained-rotor) recordings -- so it should not fail the file.
            % Note: ensembleAverageFull normalizes its output to [0,1], so this
            % only triggers on its all-zeros "no valid beats" fallback.
            ar = range(d.(avgf)(:));
            if ar < 1e-6
                add(cam + "_average", "WARN", ar, ...
                    "ensemble average degenerate (no valid beat window) — expected for unpaced/arrhythmic recordings");
            else
                % Beat-relaxation quality: does the averaged beat return toward
                % baseline?  decay-fraction = (peak - final) / (peak - min) per
                % pixel, taken over the high-amplitude (tissue) pixels.  ~1 = a
                % clean transient that relaxes; ~0 = a non-decaying plateau.  The
                % latter happens when the signal kinetics exceed the pacing cycle
                % (e.g. calcium transients at fast rates do not relax within one
                % cycle), so the ensemble average is garbage even though it passes
                % the range check (ensembleAverageFull normalizes to [0,1]).
                df = ensemble_decay_fraction(d.(avgf));
                if isnan(df)
                    add(cam + "_average", "PASS", ar, "");
                elseif df < 0.5
                    add(cam + "_average", "WARN", df, sprintf( ...
                        "ensemble average does not relax (decay-fraction %.2f) — signal kinetics likely exceed the pacing cycle; do not use the average for decay/rise/APD/tau features (use beat-windowed data)", df));
                else
                    add(cam + "_average", "PASS", df, sprintf("decay-fraction %.2f", df));
                end
            end
        end
    end

    if isempty(recs)
        add("structure", "FAIL", NaN, "no recognised CAM fields or acqFreq");
    end
end


function df = ensemble_decay_fraction(A)
%ENSEMBLE_DECAY_FRACTION  How fully the averaged beat relaxes toward baseline.
%   df = ensemble_decay_fraction(A)
%
%   A is the ensemble-average array (rows x cols x frames): one representative
%   beat per pixel.  For each high-amplitude (tissue) pixel, computes
%       (peak - final) / (peak - min)
%   where 'final' is the mean of the last few frames (diastole).  Returns the
%   median over those pixels: ~1 = a clean transient that returns to baseline,
%   ~0 = a non-decaying plateau (kinetics exceed the pacing cycle).  Returns
%   NaN when there are too few frames or no clear tissue pixels to judge.

    df = NaN;
    sz = size(A);
    if numel(sz) < 3 || sz(3) < 5
        return;
    end
    T = sz(3);
    P = double(reshape(A, [], T));            % pixels x frames
    P = P(all(isfinite(P), 2), :);            % drop background / NaN pixels
    if isempty(P)
        return;
    end
    rng_px = max(P, [], 2) - min(P, [], 2);
    if max(rng_px) <= 0
        return;
    end
    P = P(rng_px > 0.5 * max(rng_px), :);     % high-amplitude (beat-bearing) pixels
    if isempty(P)
        return;
    end
    pk  = max(P, [], 2);
    mn  = min(P, [], 2);
    nf  = min(5, T);
    fin = mean(P(:, end-nf+1:end), 2);        % diastolic level (last few frames)
    df  = median((pk - fin) ./ max(pk - mn, eps));
end


function [ratio, f_pace, f_tissue] = capture_ratio(X, analog1, fs)
%CAPTURE_RATIO  Tissue dominant frequency / pacing frequency for a CAM stack.
%   [ratio, f_pace, f_tissue] = capture_ratio(X, analog1, fs)
%
%   ratio ~ 1  => 1:1 capture (each stimulus drives one beat)
%   ratio > 1  => tissue faster than the drive (spontaneous activity or an
%                 organized tachyarrhythmia not following the pacemaker)
%   ratio < 1  => loss of capture / conduction block (e.g. 2:1)
%
%   Frequency-domain, not beat-counting, on purpose: a peak detector whose
%   minimum spacing is set from the PACING period skips every other beat of a
%   faster rhythm and reports a false ~0.75 ratio.  The dominant frequency of
%   the tissue's own spatial-mean signal has no such trap.  Toolbox-free (no
%   findpeaks / hann).  Returns NaN when pacing cannot be established.

    ratio = NaN; f_pace = NaN; f_tissue = NaN;
    if isempty(analog1) || ~isfinite(fs) || fs <= 0
        return;
    end

    % --- pacing frequency from the stimulus trace (rising edges) ---
    a  = double(analog1(:));
    mx = max(a); mn = min(a);
    if numel(a) < 10 || mx <= mn
        return;
    end
    stim = find(diff(a > (mx + mn) / 2) > 0);
    if numel(stim) < 3            % too few stimuli to define a rate
        return;
    end
    f_pace = fs / median(diff(stim));

    % --- tissue dominant frequency (NaN-safe spatial mean, windowed FFT) ---
    T = size(X, ndims(X));
    s = mean(reshape(double(X), [], T), 1, 'omitnan')';
    if ~all(isfinite(s)) || std(s) == 0
        return;
    end
    w    = 0.5 * (1 - cos(2*pi*(0:T-1)'/(T-1)));   % Hann window, toolbox-free
    s    = (s - mean(s)) .* w;
    nfft = 2^nextpow2(4*T);                        % zero-pad for resolution
    Pw   = abs(fft(s, nfft)).^2;
    fr   = (0:nfft-1)' * (fs / nfft);
    band = fr >= 1 & fr <= 40;                     % physiological band
    Pb   = Pw(band);  fb = fr(band);
    [pmax, ix] = max(Pb);
    f_tissue   = fb(ix);

    % Harmonic guard: if the tallest peak is actually the 2nd harmonic of a
    % lower fundamental (strong power near f_tissue/2), use the fundamental.
    half = abs(fb - f_tissue/2) < 0.15 * f_tissue;
    if any(half) && max(Pb(half)) >= 0.5 * pmax
        f_tissue = f_tissue / 2;
    end

    ratio = f_tissue / f_pace;
end
