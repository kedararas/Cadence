function [ok, msg, info] = check_analysis_window(win_data, mask, camera, label, pk_max, end_min, df_min)
%CHECK_ANALYSIS_WINDOW  Is an analysis window a complete, usable cardiac cycle?
%
%   [ok, msg, info] = check_analysis_window(win_data, mask, camera, label)
%   [ok, msg, info] = check_analysis_window(win_data, mask, camera, label, pk_max, end_min, df_min)
%
%   Guards every window-referenced feature (activation, APD / repolarization,
%   rise time, CV, calcium duration / decay / tau, V-Ca delay) against the
%   fast-pacing failure in which the transient does not relax within the cycle.
%   The averaged (or drawn) window then OPENS on the tail of the previous
%   transient, so the detected peak belongs to the previous beat and every
%   landmark downstream of it is wrong.
%
%   This has to be REFUSED rather than warned about, because the numbers it
%   produces look physiological: on rat calcium at 90 ms it gives a transient
%   duration of 13 ms, a decay constant of 5 ms and a NEGATIVE voltage-calcium
%   delay, none of which a reader would flag as impossible without checking.
%
%   A window is rejected when EITHER failure mode is present.
%
%   (1) WRAPPED -- both of these together:
%         peak position <= pk_max   the maximum sits at the very start of the
%                                   window (it belongs to the previous beat)
%         end level     >= end_min  the trace has not returned to its own
%                                   diastolic level by the window end
%       Either alone is normal: a slow-paced window has its peak early AND ends
%       at diastole; a window that captures the beat late has a late peak AND
%       may end high.  Only the combination means the window is not one cycle.
%       Both are measured on the mask-averaged trace normalized to its own
%       min-max range, so they are dimensionless and scale-free.
%
%   (2) NOT RELAXED -- decay fraction < df_min, per pixel, median over
%       beat-bearing pixels: (peak - final) / (peak - min).  ~1 = a clean
%       transient that returns to baseline, ~0 = a plateau.  When it is low,
%       every level referenced to the window's own range is misnamed: an
%       "80% decay" is 80% of a partial decay, so the duration returned is not
%       the duration it claims to be.  This is the same definition
%       validate_conditioned applies at the conditioning stage
%       (ensemble_decay_fraction), reused here so the two stages agree.
%
%   The two are complementary and both are needed.  Measured on the calibration
%   set below, the wrap test alone missed a window whose calcium decayed only
%   32% of the way back, and the decay test alone missed a window whose peak
%   had wrapped but which still relaxed 55% of the way.
%
%   Polarity is NOT this function's job.  Signal Conditioning owns it, and the
%   voltage camera is the one that normally needs it (voltage dyes fluoresce
%   DOWN on depolarization, calcium indicators UP), so by the time data reaches
%   feature extraction the polarity is already correct.  Do not rely on this
%   guard to catch an inversion that slipped through: tested on upright windows
%   and their exact negatives, an inverted VOLTAGE window still PASSES (the
%   maximum moves to the long diastolic stretch late in the window, so the
%   peak-position condition fails), while an inverted CALCIUM window happens to
%   be refused.  The behaviour is a side effect of AP versus CaT shape, not a
%   polarity test.  Fix polarity upstream.
%
%   Inputs
%     win_data  [rows x cols x frames] the window actually being analysed
%               (CAM<n>_average on the ensemble path, or the drawn window)
%     mask      [rows x cols] logical/numeric tissue mask, or [] for all pixels
%               (NaN-masked combo masks: pass ~isnan(combo_mask))
%     camera    camera index, for the message
%     label     what is being skipped, for the message (e.g. 'calcium duration')
%     pk_max    peak-position threshold   (default 0.35)
%     end_min   end-level threshold       (default 0.45)
%     df_min    decay-fraction threshold  (default 0.50, matching
%               validate_conditioned's WARN level)
%
%   Outputs
%     ok    true  = window usable; false = refuse and skip this feature
%     msg   '' when ok, otherwise a ready-to-print explanation
%     info  struct with .peak_position, .end_level, .decay_fraction,
%           .n_pixels and .thresholds
%
%   Calibration.  Measured on rat dual-camera recordings paced 250 -> 50 ms
%   (52 camera-windows, both cameras, baseline and drug).  With the defaults
%   the guard refuses 6 windows, every one of which has a demonstrable defect
%   (negative voltage-calcium delay, a 13 ms rat APD80, or a decay level that
%   is never reached), and passes the other 46.  Loosen end_min toward 0.6 and
%   df_min toward 0.35 to be more permissive.
%
%   The guard never throws: on any internal error it returns ok = true, so a
%   problem here can never block an extraction that would otherwise work.
%
%   See also extract_apd, extract_rise_time, extract_tau, ensembleAverageFull.

    ok   = true;
    msg  = '';
    info = struct('peak_position', NaN, 'end_level', NaN, 'decay_fraction', NaN, ...
                  'n_pixels', 0, 'thresholds', [NaN NaN NaN]);

    try
        if nargin < 2, mask   = [];  end
        if nargin < 3, camera = NaN; end
        if nargin < 4, label  = '';  end
        if nargin < 5 || isempty(pk_max),  pk_max  = 0.35; end
        if nargin < 6 || isempty(end_min), end_min = 0.45; end
        if nargin < 7 || isempty(df_min),  df_min  = 0.50; end
        info.thresholds = [pk_max end_min df_min];

        if isempty(win_data) || ndims(win_data) ~= 3, return; end
        [nr, nc, T] = size(win_data);
        if T < 20, return; end                  % too short to judge shape

        X = reshape(double(win_data), nr * nc, T);
        finite_px = all(isfinite(X), 2);
        if isempty(mask) || ~any(mask(:))
            keep = finite_px;
        else
            keep = logical(mask(:)) & finite_px;
        end
        info.n_pixels = nnz(keep);
        if info.n_pixels < 50, return; end       % not enough tissue to judge

        % Mask-averaged trace: a single high-SNR view of the window's shape.
        y = mean(X(keep, :), 1, 'omitnan');
        y = movmean(y, 5);
        y = y(3:end-2);                          % drop the smoothing edges
        span = max(y) - min(y);
        if ~(span > 0), return; end              % flat window: nothing to judge

        [~, pk] = max(y);
        info.peak_position = pk / numel(y);
        info.end_level     = (y(end) - min(y)) / span;

        % Criterion 2, per-pixel: how fully the averaged beat relaxes.  Same
        % definition validate_conditioned already applies at the conditioning
        % stage (ensemble_decay_fraction), reused here so the two stages agree.
        info.decay_fraction = decay_fraction(X(keep, :));

        wrapped     = info.peak_position <= pk_max && info.end_level >= end_min;
        not_relaxed = isfinite(info.decay_fraction) && info.decay_fraction < df_min;
        ok = ~(wrapped || not_relaxed);

        if wrapped
            msg = sprintf(['CAM%d %s SKIPPED: the analysis window is not a complete cycle ' ...
                '(peak at %.0f%% of the window, end still %.0f%% above diastole). The ' ...
                'transient has not relaxed within the cycle, so the window opens on the ' ...
                'previous beat and every window-referenced value would be measured against ' ...
                'the wrong landmarks. Draw a signal window over one full transient, or ' ...
                'analyse this rate beat-by-beat.'], ...
                camera, label, 100 * info.peak_position, 100 * info.end_level);
        elseif not_relaxed
            msg = sprintf(['CAM%d %s SKIPPED: the averaged beat does not relax within the ' ...
                'cycle (decay-fraction %.2f, i.e. it falls only %.0f%% of the way back to ' ...
                'diastole). Levels referenced to the window range are then not the levels ' ...
                'they are named after -- an "80%% decay" is 80%% of a partial decay. Use ' ...
                'beat-windowed data at this rate.'], ...
                camera, label, info.decay_fraction, 100 * info.decay_fraction);
        end

    catch
        ok  = true;      % a guard must never be the reason an extraction fails
        msg = '';
    end
end


% =========================================================================
function df = decay_fraction(P)
%DECAY_FRACTION  Median (peak - final) / (peak - min) over beat-bearing pixels.
%   Mirrors validate_conditioned/ensemble_decay_fraction so the conditioning
%   stage and this guard cannot disagree about the same recording.
    df = NaN;
    if isempty(P) || size(P, 2) < 5, return; end
    rng_px = max(P, [], 2) - min(P, [], 2);
    if max(rng_px) <= 0, return; end
    P = P(rng_px > 0.5 * max(rng_px), :);      % high-amplitude (beat-bearing) pixels
    if isempty(P), return; end
    pk  = max(P, [], 2);
    mn  = min(P, [], 2);
    nf  = min(5, size(P, 2));
    fin = mean(P(:, end-nf+1:end), 2);         % diastolic level (last few frames)
    df  = median((pk - fin) ./ max(pk - mn, eps));
end
