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
%     <CAM>_average      ensemble average has dynamic range (catches the
%                        all-zeros / no-valid-beat ensemble-average failure).
%
%   See qc_check for the record shape and status semantics.

    recs = repmat(qc_check("","","" ,"PASS",NaN,""), 0, 1);   % empty typed array
    name = string(recording);

    % ---- nested append helper (shares recs/name) ----
    function add(check, status, value, message)
        recs(end+1, 1) = qc_check(name, "conditioning", check, status, value, message);
    end

    % ---- acquisition frame rate ----
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
