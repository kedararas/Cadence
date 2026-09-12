function [needs_inversion, confidence, info] = check_polarity_paced(data, pacing, fs, mask, min_beats)
%CHECK_POLARITY_PACED  Signal polarity from the pacing stimulus, not from shape.
%
%   [needs_inversion, confidence, info] = check_polarity_paced(data, pacing, fs)
%   [needs_inversion, confidence, info] = check_polarity_paced(data, pacing, fs, mask, min_beats)
%
%   After a pacing stimulus captures, a correctly oriented optical signal must
%   deflect UP: voltage depolarizes, calcium rises.  Anchoring the decision to
%   the stimulus makes it independent of the signal's SHAPE, which is what the
%   shape heuristics in check_signal_inversion depend on and what fails at fast
%   pacing rates.
%
%   WHY THIS EXISTS.  check_signal_inversion combines a histogram-mode test
%   ("a correct signal spends most of its time near rest") with a slope-
%   asymmetry test ("the fast edge points up"), weighting slope double because
%   it is the duty-cycle-robust one.  Both assumptions break together as the
%   cycle shortens: the transient fills more of the cycle, so the signal spends
%   most of its time ELEVATED and the histogram test votes "inverted", while
%   the rise and decay speeds converge, so the slope test collapses toward zero
%   and stops outvoting it.  The duty-cycle-biased criterion is then left
%   deciding alone.  Calcium crosses that threshold first because its transient
%   is the longer one.
%
%   Measured on rat dual-camera recordings: within one heart and one camera,
%   consecutive recordings minutes apart returned calcium deflections of +0.86,
%   +0.83, -0.82, -0.80, +0.76, -0.70 as the cycle length stepped down.  The
%   magnitude is preserved while the sign alternates, which no physiological
%   change can produce -- it is the per-recording detector flip-flopping.
%
%   METHOD.  Over the tissue mask, take the spatial-mean trace.  For each
%   captured beat, compare the level just BEFORE the stimulus with the extremes
%   reached during the beat:
%       d_up = max(beat) - level_before        how far it rises
%       d_dn = level_before - min(beat)        how far it falls
%   Upright when the rise dominates.  This holds regardless of duty cycle: with
%   incomplete relaxation the pre-stimulus level is still the cycle's minimum
%   for an upright signal and its maximum for an inverted one, so the contrast
%   grows rather than shrinks at fast rates.
%
%   Inputs
%     data       [rows x cols x frames] optical stack (conditioned)
%     pacing     analog stimulus trace or binary onset vector; [] if none
%     fs         sampling rate (Hz)
%     mask       [rows x cols] tissue mask, or [] for all finite pixels
%     min_beats  minimum captured beats required (default 5)
%
%   Outputs
%     needs_inversion  true when the signal should be flipped
%     confidence       (d_up - d_dn) / (d_up + d_dn), in [-1, 1].  POSITIVE
%                      means upright, negative means inverted, and the
%                      magnitude is how decisive the evidence is.  NaN when
%                      there is no usable pacing channel -- the caller should
%                      then fall back to check_signal_inversion.
%     info             struct: .n_beats, .d_up, .d_dn, .method
%
%   Never throws: any internal failure returns NaN confidence so the caller
%   falls back rather than acting on a bad decision.
%
%   See also check_signal_inversion, extract_beat_frames.

    needs_inversion = false;
    confidence      = NaN;
    info = struct('n_beats', 0, 'd_up', NaN, 'd_dn', NaN, 'method', 'none');

    try
        if nargin < 4, mask = []; end
        if nargin < 5 || isempty(min_beats), min_beats = 5; end
        if isempty(data) || ndims(data) ~= 3 || isempty(pacing) || ~isfinite(fs) || fs <= 0
            return;
        end

        [nr, nc, T] = size(data);
        p = double(pacing(:));
        if numel(p) ~= T || ~any(isfinite(p)), return; end

        % ---- stimulus onsets (same convention as extract_beat_frames) --------
        if all(p == 0 | p == 1)
            onsets = find(p);
        else
            lvl    = 0.5 * max(p);
            onsets = find(p(2:end) >= lvl & p(1:end-1) < lvl) + 1;
        end
        if numel(onsets) < min_beats + 2, return; end
        ibi = round(median(diff(onsets)));
        if ~(ibi > 2), return; end

        % ---- mask-averaged trace ---------------------------------------------
        X = reshape(double(data), nr * nc, T);
        keep = all(isfinite(X), 2);
        if ~isempty(mask) && any(mask(:))
            keep = keep & logical(mask(:));
        end
        if nnz(keep) < 20, return; end
        y = mean(X(keep, :), 1, 'omitnan');
        clear X

        % ---- per-beat rise versus fall, referenced to the pre-stimulus level --
        pre = max(2, round(0.005 * fs));                 % 5 ms before the stimulus
        use = onsets(onsets > pre + 1 & onsets + ibi - 1 <= T);
        if numel(use) > 4, use = use(2:end-1); end       % drop the edge beats
        if numel(use) < min_beats, return; end

        up = nan(numel(use), 1);  dn = up;
        for k = 1:numel(use)
            o    = use(k);
            base = mean(y(o - pre : o));
            seg  = y(o : o + ibi - 1);
            up(k) = max(seg) - base;
            dn(k) = base - min(seg);
        end

        info.n_beats = numel(use);
        info.d_up    = median(up, 'omitnan');
        info.d_dn    = median(dn, 'omitnan');
        info.method  = 'post-stimulus deflection';

        den = info.d_up + info.d_dn;
        if ~(den > 0), return; end
        confidence      = (info.d_up - info.d_dn) / den;
        needs_inversion = confidence < 0;

    catch
        needs_inversion = false;
        confidence      = NaN;
    end
end
