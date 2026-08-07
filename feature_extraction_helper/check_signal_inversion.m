function [needs_inversion, confidence, diagnostics] = check_signal_inversion(data, fs, varargin)
% CHECK_SIGNAL_INVERSION  Determine whether voltage/optical signals need inverting.
%
%   [needs_inversion, confidence, diagnostics] = check_signal_inversion(data, fs)
%
%   Works for both 1-D traces and 3-D optical mapping volumes (rows×cols×frames),
%   and for both PACED and ARRHYTHMIA (AF/VF) signals — the two criteria below
%   are orientation cues that hold regardless of regularity, so no signal-type
%   detection or branching is needed.
%
%   Two complementary criteria
%     1. Histogram mode   — a correctly-oriented signal spends most of its
%                           time near rest, so the amplitude histogram peaks
%                           in the LOWER half of the range.  An inverted
%                           signal pushes the mode into the upper half.
%     2. Slope asymmetry  — depolarisation (Na-channel kinetics) is faster
%                           than repolarisation, so dV/dt_max sits on the
%                           upstroke.  This is the most physically reliable
%                           cue and is robust to plateau / duty-cycle and to
%                           baseline offset, so it carries extra weight.
%
%   Each criterion reports a signed confidence in [-1, +1] (positive =
%   correctly oriented, negative = inverted).  The recommendation is a
%   magnitude-weighted sum of the two confidences, so a near-zero
%   "coin-flip" criterion cannot override a decisive one.
%
%   For a 3-D stack the function pools all in-ROI pixels into a SINGLE global
%   trace (SNR-weighted spatial mean) and tests that one high-SNR trace, then
%   applies the result to the whole field.  Polarity is a global property of
%   the acquisition, so pooling maximises SNR for the decision.
%
%   Inputs
%     data        [N×1] or [rows×cols×frames]  voltage or fluorescence
%                 (ideally baseline-corrected / filtered / ensemble-averaged)
%     fs          scalar   sampling rate (Hz)
%
%   Name/Value options
%     'Method'         'vote' (default) | 'histmode' | 'slope'
%                      'vote'     — magnitude-weighted combination of both
%                      'histmode' — histogram-mode criterion only
%                      'slope'    — slope-asymmetry criterion only
%     'Mask'           [rows×cols] logical/numeric ROI (3-D input only).
%                      Pixels outside the mask are excluded from pooling.
%                      Default: all finite, non-flat pixels.
%     'PixelWeighting' 'snr' (default) | 'uniform'   (3-D input only)
%                      'snr'     — weight each pixel by its SNR when pooling,
%                                  so noisy/background pixels contribute less.
%                      'uniform' — plain mean over the ROI pixels.
%     'Plot'           true/false (default false)
%
%   Outputs
%     needs_inversion   logical scalar
%     confidence        agreement of the weighted evidence, in [0, 1]
%                       (1 = both criteria agree, →0 = evidence is split)
%     diagnostics       struct with per-criterion scores, confidences and votes
%
%   Example — paced trace:
%     [inv, conf] = check_signal_inversion(trace, 1000);
%     if inv, trace = -trace; end
%
%   Example — AF optical mapping volume:
%     [inv, conf] = check_signal_inversion(cmos_data, frame_rate, 'Plot', true);
%     if inv, cmos_data = -cmos_data; end

% ---------- parse inputs ----------
p = inputParser;
p.addRequired('data', @(v) isnumeric(v) && (isvector(v) || ndims(v)==3));
p.addRequired('fs',   @(v) isscalar(v)  && v > 0);
p.addParameter('Method',         'vote', @(s) ismember(lower(s), {'vote','histmode','slope'}));
p.addParameter('Mask',           [],     @(v) isempty(v) || isnumeric(v) || islogical(v));
p.addParameter('PixelWeighting', 'snr',  @(s) ismember(lower(s), {'snr','uniform'}));
p.addParameter('Plot',           false,  @islogical);
p.parse(data, fs, varargin{:});
o = p.Results;

% ---------- reduce to the single trace tested for inversion ----------
% 1-D input  : test the trace as-is.
% 3-D stack  : pool all in-ROI pixels into ONE global trace.  Polarity is a
%              global property, so a high-SNR spatial mean is the best input.
roi_npix = NaN; roi_snr_mean = NaN;
if isvector(data)
    traces = double(data(:));
else
    [traces, roi_npix, roi_snr_mean] = build_roi_trace(double(data), fs, o.Mask, o.PixelWeighting);
end

% ---------- guard: flat / all-zero global trace ----------
% Occurs when the caller passes an uninitialised placeholder (e.g. an
% all-zero *_average field before ensemble averaging has been run), or
% when no valid pixels were found.  Both criteria collapse to conf = 0,
% the weighted sum = 0, and 0 < 0 = false — silently giving the wrong
% answer.  Catch this early and return with low confidence + a warning.
trace_range = max(traces) - min(traces);
if trace_range < 1e-6
    warning('check_signal_inversion:flatTrace', ...
        ['The global trace has no dynamic range (range = %.2e).\n' ...
         'This usually means an uninitialised placeholder (e.g. CAM_average\n' ...
         'before ensemble averaging) was passed instead of the conditioned\n' ...
         'data array.  Pass the main CAM1/CAM2/... field instead.\n' ...
         'Returning needs_inversion = false with confidence = 0.'], trace_range);
    needs_inversion = false;
    confidence      = 0;
    diagnostics     = struct( ...
        'histmode_score', 0, 'histmode_conf', 0, 'histmode_vote_inv', false, ...
        'slope_score',    0, 'slope_conf',    0, 'slope_vote_inv',    false, ...
        'weighted_total', 0, 'needs_inversion', false, 'confidence', 0, ...
        'roi_npix', roi_npix, 'roi_snr_mean', roi_snr_mean, ...
        'global_trace', traces(:, 1));
    return;
end

% ---------- compute the two criteria ----------
c_hist  = criterion_histogram_mode(traces);
c_slope = criterion_slope(traces, fs);
c_hist.name  = 'histogram_mode';
c_slope.name = 'slope_asym';

% ---------- combine ----------
% .conf is a signed, normalised confidence in [-1, +1]:
%   conf > 0  → evidence the signal is correctly oriented
%   conf < 0  → evidence the signal is inverted
% Magnitude-weighted sum, with the kinetic slope cue weighted more heavily
% because it is the most reliable orientation indicator.
SLOPE_WEIGHT = 2;
switch lower(o.Method)
    case 'histmode', total = c_hist.conf;   wabs = abs(c_hist.conf);
    case 'slope',    total = c_slope.conf;  wabs = abs(c_slope.conf);
    otherwise  % weighted evidence
        w     = [1, SLOPE_WEIGHT];
        cf    = [c_hist.conf, c_slope.conf];
        total = sum(w .* cf);
        wabs  = sum(w .* abs(cf));
end

needs_inversion = total < 0;
confidence      = abs(total) / (wabs + eps);   % [0,1]: 1 = full agreement

% ---------- diagnostics ----------
diagnostics = struct( ...
    'histmode_score',    c_hist.score, ...
    'histmode_conf',     c_hist.conf, ...
    'histmode_vote_inv', c_hist.vote_inv, ...
    'slope_score',       c_slope.score, ...
    'slope_conf',        c_slope.conf, ...
    'slope_vote_inv',    c_slope.vote_inv, ...
    'weighted_total',    total, ...
    'needs_inversion',   needs_inversion, ...
    'confidence',        confidence, ...
    'roi_npix',          roi_npix, ...
    'roi_snr_mean',      roi_snr_mean, ...
    'global_trace',      traces(:, 1));

% ---------- report ----------
if ~isnan(roi_npix)
    fprintf('Pooled %d ROI pixels into a global trace (mean SNR = %.1f, weighting = %s)\n', ...
            roi_npix, roi_snr_mean, lower(o.PixelWeighting));
end
fprintf('Inversion check:  needs_inversion = %d   confidence = %.0f%%\n', ...
        needs_inversion, confidence*100);
fprintf('  %-16s score = %+.3f  conf = %+.2f  vote_inv = %d\n', ...
        c_hist.name,  c_hist.score,  c_hist.conf,  c_hist.vote_inv);
fprintf('  %-16s score = %+.3f  conf = %+.2f  vote_inv = %d\n', ...
        c_slope.name, c_slope.score, c_slope.conf, c_slope.vote_inv);
if strcmpi(o.Method, 'vote')
    fprintf('  weighted evidence total = %+.3f  (slope weight = %d)\n', total, SLOPE_WEIGHT);
end

% ---------- plot ----------
if o.Plot
    ref_trace = traces(:, min(end, max(1, round(size(traces,2)/2))));  % median column
    plot_summary(ref_trace, fs, needs_inversion, diagnostics);
end
end


% =========================================================================
%  ROI pooling (3-D input)
% =========================================================================

function [global_trace, n_pix, snr_mean] = build_roi_trace(data, fs, mask, weighting)
% BUILD_ROI_TRACE  Pool a [rows x cols x frames] stack into one global trace.
%
%   The trace is a per-frame spatial mean over valid ROI pixels.  With
%   'snr' weighting each pixel is weighted by its SNR, so clean pixels
%   dominate and noisy/background pixels contribute little — giving the
%   highest-SNR single trace for the polarity decision.
%
%   Valid pixels = inside Mask (if given) AND finite over all frames AND
%   non-flat (dynamic range > eps).

    [nr, nc, nf] = size(data);
    X = reshape(data, nr*nc, nf);          % pixels x frames

    % ---- per-pixel SNR: (peak-to-peak of smoothed) / std(residual) ----
    win    = max(3, round(0.01 * nf));
    Xs     = movmean(X, win, 2);           % smooth along time (per pixel)
    amp    = max(Xs, [], 2) - min(Xs, [], 2);
    noise  = std(X - Xs, 0, 2);
    snr    = amp ./ (noise + eps);         % pixels x 1

    % ---- valid pixel mask ----
    valid = all(isfinite(X), 2) & amp > eps & noise > 0;
    if ~isempty(mask)
        valid = valid & logical(mask(:));
    end
    if ~any(valid)                          % fall back to anything finite
        valid = all(isfinite(X), 2);
    end

    % ---- pooling weights (valid pixels only) ----
    % Operate on the valid subset only: those pixels are guaranteed finite,
    % so the weighted sum cannot be poisoned by 0*NaN from dead pixels.
    Xv = X(valid, :);
    if strcmpi(weighting, 'uniform')
        wv = ones(size(Xv, 1), 1);
    else                                    % 'snr'
        wv = max(snr(valid), 0);
        if ~any(wv > 0), wv = ones(size(Xv, 1), 1); end   % guard: all-zero SNR
    end

    % ---- weighted spatial mean per frame ----
    global_trace = (wv' * Xv) ./ (sum(wv) + eps);   % 1 x frames
    global_trace = global_trace(:);                 % column

    n_pix    = sum(valid);
    snr_mean = mean(snr(valid), 'omitnan');
end


% =========================================================================
%  Criteria
% =========================================================================

function r = criterion_histogram_mode(traces)
% The mode of the amplitude histogram indicates where the signal spends most
% of its time:
%
%   correct orientation:  mode in lower half of amplitude range  →  score > 0
%   inverted:             mode in upper half of amplitude range  →  score < 0
%
% mode_position = (mode_amplitude - min) / range  ∈ [0, 1]
% score = 0.5 - mode_position  (positive means mode in lower half)
    num_traces = size(traces, 2);
    scores     = zeros(1, num_traces);
    for i = 1:num_traces
        col = traces(:, i);
        % Robust amplitude bounds: use the 1st/99th percentiles rather than
        % raw min/max so a single outlier or a startup/edge transient cannot
        % define the range.  With raw min/max, one artifact sample at the top
        % stretches the range and drags the computed mode position into the
        % lower half — flipping this criterion's vote to "not inverted".
        col_lo    = prctile(col,  1);
        col_hi    = prctile(col, 99);
        amp_range = col_hi - col_lo;
        if amp_range < eps
            scores(i) = 0;
            continue;
        end
        col = min(max(col, col_lo), col_hi);   % clip outliers into the range
        [counts, edges] = histcounts(col, 20);
        [~, idx]        = max(counts);
        mode_amp        = (edges(idx) + edges(idx+1)) / 2;
        mode_pos        = (mode_amp - col_lo) / amp_range;  % 0..1
        scores(i)       = 0.5 - mode_pos;  % >0 → mode in lower half → correct
    end
    r.score    = median(scores, 'omitnan');
    r.conf     = 2 * r.score;              % scores in [-0.5,0.5] → conf in [-1,1]
    r.vote_inv = r.score < 0;
end

function r = criterion_slope(traces, fs)
% dV/dt_max should be on the upstroke (positive slope).
% Works for paced signals (large, obvious upstroke) and arrhythmia
% (smaller individual events but the same kinetic asymmetry).
%
% The signal is lightly smoothed (~6 ms moving average) before
% differentiating.  The raw single-sample max derivative is dominated by
% high-frequency noise spikes — which are roughly symmetric and wash out
% the true up/down kinetic asymmetry — so on noisy traces the raw slope
% collapses toward zero (and can even flip sign).  Smoothing first recovers
% the genuine AP up/down slopes while leaving clean signals essentially
% unchanged.
    dt       = 1 / fs;
    win      = max(3, round(0.006 * fs));    % ~6 ms moving average
    smoothed = movmean(traces, win, 1);      % along time (dim 1)
    dv       = diff(smoothed, 1, 1) / dt;
    max_pos = median(max( dv, [], 1), 'omitnan');
    max_neg = median(max(-dv, [], 1), 'omitnan');
    r.score    = max_pos - max_neg;
    r.conf     = r.score / (max_pos + max_neg + eps);  % slope asymmetry [-1,1]
    r.vote_inv = r.score < 0;
end


% =========================================================================
%  Plot
% =========================================================================

function plot_summary(trace, fs, needs_inversion, d)
    t = (0:numel(trace)-1) / fs;
    figure('Color','w','Position',[100 100 1100 420]);

    % Raw vs corrected
    subplot(1,3,1);
    plot(t, trace, 'k', 'LineWidth', 1); hold on;
    if needs_inversion
        plot(t, -trace, 'r--', 'LineWidth', 1);
        legend('Raw','Inverted','Location','best');
    end
    xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('invert = %d', needs_inversion));
    grid on;

    % Amplitude histogram
    subplot(1,3,2);
    histogram(trace, 30, 'FaceColor', [.6 .6 .9], 'EdgeColor','none');
    xline(prctile(trace,10), 'b--', 'P10');
    xline(prctile(trace,50), 'k-',  'P50');
    xline(prctile(trace,90), 'r--', 'P90');
    xlabel('Amplitude'); ylabel('Count');
    title(sprintf('Histogram  histmode score = %.3f', d.histmode_score));
    grid on;

    % dV/dt
    subplot(1,3,3);
    dt = 1 / fs;
    dv = diff(trace) / dt;
    plot((0:numel(dv)-1)/fs, dv, 'Color',[.3 .5 .3], 'LineWidth', 0.8);
    yline(0, 'k--');
    xlabel('Time (s)'); ylabel('dV/dt');
    title(sprintf('Slope  score = %.3f', d.slope_score));
    grid on;

    sgtitle(sprintf('Inversion check — invert=%d  confidence=%.0f%%', ...
            needs_inversion, d.confidence*100));
end
