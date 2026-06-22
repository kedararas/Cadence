function [df, ri, oi, spec] = cardiacSpectralMetrics(x, fs, varargin)
%CARDIACSPECTRALMETRICS  Dominant frequency, regularity index and
% organization index of a cardiac action-potential (or electrogram) signal.
%
%   [DF, RI, OI, SPEC] = cardiacSpectralMetrics(X, FS)
%
%   Accepts either:
%     X   [N x 1]            1-D time series  → scalar DF, RI, OI outputs
%     X   [rows x cols x N]  optical-mapping volume → rows×cols map outputs
%
%   For 3-D input every pixel is processed independently using a single
%   batched pwelch call (one column per pixel). Pixels that are all-NaN
%   (background mask) are skipped and returned as NaN in the output maps.
%   spec.pxx is returned as (nfft/2+1) × rows × cols; spec.xProc is omitted
%   (storing a preprocessed copy of the full volume is impractical).
%
%   Computes three standard spectral-organization metrics:
%
%     DF  Dominant Frequency (Hz)
%         The frequency that carries the peak power-spectral-density
%         within a configurable search band.
%
%     RI  Regularity Index (dimensionless, 0..1)
%         RI = P(DF +/- BW) / P(search band)
%         How tightly the signal's power is concentrated at DF. A perfect
%         sinusoid at DF gives RI -> 1, broadband noise gives RI -> 0.
%
%     OI  Organization Index (dimensionless, 0..1), Everett et al. 2001:
%         OI = sum_{k=1..K} P(k*DF +/- BW) / P(search band)
%         Shares the fundamental with RI but also counts the harmonic
%         content of a periodic non-sinusoidal signal -- e.g. an action
%         potential train, whose square-ish shape puts substantial power
%         in 2*DF, 3*DF, etc. A highly organized AP train therefore has
%         OI >> RI; a noisy/fibrillatory signal has OI ~= RI ~= low.
%
%   Inputs
%     X   [N x 1] or [rows x cols x N] real-valued signal(s)
%     FS  scalar sampling rate in Hz
%
%   Name / Value options
%     'Band'            [fLo fHi] DF search band        default [3 20]
%                       Use [3 12] for AF electrograms, [4 10] for rabbit
%                       paced APs, [5 20] for VF.
%     'PeakBW'          half-width (Hz) of peak window  default 0.75
%                       Width BW around DF and each k*DF used for RI/OI.
%     'Harmonics'       number of harmonics for OI      default 4
%     'NFFT'            Welch FFT length (samples)      default 4*fs
%     'WindowSec'       Welch window length (s)         default 2
%     'OverlapFrac'     Welch overlap fraction          default 0.5
%     'Preproc'         'none' | 'detrend' | 'botteron' default 'detrend'
%                       'botteron' = bandpass 40-250 Hz, rectify,
%                       low-pass 20 Hz (Botteron & Smith 1995). Designed
%                       for fractionated atrial electrograms. For clean
%                       optical AP traces use 'detrend' (just remove DC
%                       and linear trend).
%     'TotalPowerBand'  [fLo fHi] for RI/OI denominator default = Band
%                       Set to e.g. [0.5 30] to normalize to the full
%                       physiological band rather than just the search.
%     'Plot'            true/false                      default false
%                       For 3-D input plots the spatially-averaged PSD.
%
%   Outputs — 1-D input
%     DF    dominant frequency (Hz)
%     RI    regularity index (0..1)
%     OI    organization index (0..1)
%     SPEC  struct with fields:
%             .f         frequency vector (Hz)
%             .pxx       PSD estimate (nfft/2+1 × 1)
%             .df        same as DF
%             .ri, .oi   same as RI, OI
%             .peakBW, .harmonics, .band, .totalBand
%             .xProc     preprocessed signal
%
%   Outputs — 3-D input
%     DF    rows × cols dominant-frequency map (Hz)
%     RI    rows × cols regularity-index map
%     OI    rows × cols organization-index map
%     SPEC  struct with fields:
%             .f         frequency vector (Hz)
%             .pxx       (nfft/2+1) × rows × cols PSD volume
%             .df, .ri, .oi  same as map outputs
%             .peakBW, .harmonics, .band, .totalBand
%             (xProc omitted for 3-D — too large to store per pixel)
%
%   Example — 1-D
%     load('ap_trace.mat');               % x, fs
%     [df, ri, oi] = cardiacSpectralMetrics(x, fs, ...
%                        'Band', [3 12], 'Plot', true);
%
%   Example — 3-D optical mapping
%     % data: rows × cols × frames double array
%     [df_map, ri_map, oi_map, spec] = cardiacSpectralMetrics(data, fs, ...
%                        'Band', [3 12], 'Preproc', 'detrend');
%     imagesc(df_map); colorbar; title('DF map (Hz)');
%
%   References
%     Botteron GW, Smith JM. A technique for measurement of the extent of
%       spatial organization of atrial activation during atrial
%       fibrillation in the intact human heart. IEEE Trans Biomed Eng 42,
%       1995.
%     Everett TH et al. Assessment of global atrial fibrillation
%       organization to optimize timing of atrial defibrillation.
%       Circulation 103, 2001.
%     Sanders P et al. Spectral analysis identifies sites of high-
%       frequency activity maintaining atrial fibrillation in humans.
%       Circulation 112, 2005.
%
% ------------------------------------------------------------------------

    % ---------- parse inputs ----------
    p = inputParser;
    p.addRequired('x',  @(v) isnumeric(v) && (isvector(v) || ndims(v) == 3));
    p.addRequired('fs', @(v) isnumeric(v) && isscalar(v) && v > 0);
    p.addParameter('Band',           [3 20],    @(v) numel(v)==2 && v(2)>v(1));
    p.addParameter('PeakBW',         0.75,      @(v) v > 0);
    p.addParameter('Harmonics',      4,         @(v) v >= 1 && v == round(v));
    p.addParameter('NFFT',           [],        @(v) isempty(v) || (v > 0 && v == round(v)));
    p.addParameter('WindowSec',      2,         @(v) v > 0);
    p.addParameter('OverlapFrac',    0.5,       @(v) v >= 0 && v < 1);
    p.addParameter('Preproc',        'detrend', ...
        @(s) any(strcmpi(s, {'none','detrend','botteron'})));
    p.addParameter('TotalPowerBand', [],        @(v) isempty(v) || numel(v) == 2);
    p.addParameter('Plot',           false,     @islogical);
    p.parse(x, fs, varargin{:});
    o = p.Results;

    if isempty(o.TotalPowerBand), o.TotalPowerBand = o.Band; end
    if isempty(o.NFFT),           o.NFFT = max(256, 2^nextpow2(4*fs)); end

    % ---------- dispatch ----------
    if ~isvector(x) && ndims(double(x)) == 3
        [df, ri, oi, spec] = compute3D(double(x), fs, o);
    else
        [df, ri, oi, spec] = compute1D(double(x(:)), fs, o);
    end
end


% ========================================================================
function [df, ri, oi, spec] = compute1D(x, fs, o)
% Original single-channel path — unchanged.

    % --- preprocessing ---
    xProc = applyPreproc(x, fs, o.Preproc);

    % --- Welch PSD ---
    winN  = min(numel(xProc), max(64, round(o.WindowSec * fs)));
    nOvlp = round(o.OverlapFrac * winN);
    [pxx, f] = pwelch(xProc, hann(winN), nOvlp, o.NFFT, fs);

    % --- dominant frequency ---
    inBand = f >= o.Band(1) & f <= o.Band(2);
    assert(any(inBand), 'Search band [%g %g] Hz contains no FFT bins.', ...
           o.Band(1), o.Band(2));
    fBand   = f(inBand);
    pxxBand = pxx(inBand);
    [~, iPk] = max(pxxBand);
    df = fBand(iPk);

    % --- RI and OI ---
    inTot  = f >= o.TotalPowerBand(1) & f <= o.TotalPowerBand(2);
    totalP = trapz(f(inTot), pxx(inTot));

    ri = bandPower(f, pxx, df, o.PeakBW) / totalP;

    oi = 0;
    for k = 1:o.Harmonics
        fk = k * df;
        if fk + o.PeakBW > f(end), break; end
        oi = oi + bandPower(f, pxx, fk, o.PeakBW);
    end
    oi = oi / totalP;

    % --- pack struct ---
    spec = struct( ...
        'f',         f, ...
        'pxx',       pxx, ...
        'df',        df, ...
        'ri',        ri, ...
        'oi',        oi, ...
        'peakBW',    o.PeakBW, ...
        'harmonics', o.Harmonics, ...
        'band',      o.Band, ...
        'totalBand', o.TotalPowerBand, ...
        'xProc',     xProc);

    if o.Plot
        plotSpectrum(spec);
    end
end


% ========================================================================
function [df_map, ri_map, oi_map, spec] = compute3D(x, fs, o)
% Batched 3-D path: single pwelch call across all valid pixels.
%
% x is rows × cols × frames.
% Valid pixels are those with no NaN samples (background mask is all-NaN).
% All preprocessing functions (detrend, botteronSmith) accept matrices
% with one signal per column, so they apply to every pixel in one call.

    [num_rows, num_cols, num_frames] = size(x);
    num_pixels = num_rows * num_cols;

    % Reshape to frames × pixels — pwelch expects samples along rows
    x2d = reshape(permute(x, [3 1 2]), num_frames, num_pixels);

    % Pixels with any NaN frame are excluded (background mask)
    valid_mask = ~any(isnan(x2d), 1);   % 1 × num_pixels logical

    % --- preprocess valid pixels only ---
    x2d_valid = applyPreproc(x2d(:, valid_mask), fs, o.Preproc);

    % --- single pwelch call for all valid pixels ---
    % pwelch accepts a matrix: each column is an independent channel.
    % Output pxx_valid: (nfft/2+1) × num_valid_pixels
    winN  = min(num_frames, max(64, round(o.WindowSec * fs)));
    nOvlp = round(o.OverlapFrac * winN);
    [pxx_valid, f] = pwelch(x2d_valid, hann(winN), nOvlp, o.NFFT, fs);

    nF  = numel(f);
    pxx = NaN(nF, num_pixels);
    pxx(:, valid_mask) = pxx_valid;

    % --- dominant frequency — vectorised max across all pixels at once ---
    inBand = f >= o.Band(1) & f <= o.Band(2);
    assert(any(inBand), 'Search band [%g %g] Hz contains no FFT bins.', ...
           o.Band(1), o.Band(2));
    fBand    = f(inBand);
    [~, iPk] = max(pxx(inBand, :), [], 1);     % 1 × num_pixels
    df_vec   = fBand(iPk);                      % 1 × num_pixels
    df_vec(~valid_mask) = NaN;

    % --- total power per pixel (denominator for RI and OI) ---
    inTot  = f >= o.TotalPowerBand(1) & f <= o.TotalPowerBand(2);
    totalP = trapz(f(inTot), pxx(inTot, :), 1);  % 1 × num_pixels

    % --- RI and OI — per-pixel loop (lightweight after vectorised pwelch) ---
    % RI and OI require per-pixel df to define the integration window,
    % so full vectorisation would need per-column interpolation. The loop
    % cost is negligible compared to the pwelch call above.
    ri_vec = NaN(1, num_pixels);
    oi_vec = NaN(1, num_pixels);

    valid_idx = find(valid_mask);
    for vi = 1:numel(valid_idx)
        px  = valid_idx(vi);
        tp  = totalP(px);
        if tp <= 0, continue; end
        df_px = df_vec(px);
        ri_vec(px) = bandPower(f, pxx(:, px), df_px, o.PeakBW) / tp;
        oi_px = 0;
        for k = 1:o.Harmonics
            fk = k * df_px;
            if fk + o.PeakBW > f(end), break; end
            oi_px = oi_px + bandPower(f, pxx(:, px), fk, o.PeakBW);
        end
        oi_vec(px) = oi_px / tp;
    end

    % --- reshape vectors back to spatial maps ---
    df_map = reshape(df_vec, num_rows, num_cols);
    ri_map = reshape(ri_vec, num_rows, num_cols);
    oi_map = reshape(oi_vec, num_rows, num_cols);

    % --- pack struct (xProc omitted — storing all-pixel preprocessed data
    %     would replicate the entire input volume) ---
    spec = struct( ...
        'f',         f, ...
        'pxx',       reshape(pxx, nF, num_rows, num_cols), ...
        'df',        df_map, ...
        'ri',        ri_map, ...
        'oi',        oi_map, ...
        'peakBW',    o.PeakBW, ...
        'harmonics', o.Harmonics, ...
        'band',      o.Band, ...
        'totalBand', o.TotalPowerBand);

    % --- optional plot: spatially-averaged PSD ---
    if o.Plot
        mean_pxx = mean(pxx_valid, 2);   % average over valid pixels only
        s_avg = struct( ...
            'f',         f, ...
            'pxx',       mean_pxx, ...
            'df',        median(df_vec(valid_mask), 'omitnan'), ...
            'ri',        median(ri_vec(valid_mask), 'omitnan'), ...
            'oi',        median(oi_vec(valid_mask), 'omitnan'), ...
            'peakBW',    o.PeakBW, ...
            'harmonics', o.Harmonics, ...
            'band',      o.Band, ...
            'totalBand', o.TotalPowerBand);
        plotSpectrum(s_avg);
        title(sprintf( ...
            'Spatially-averaged PSD  (%d pixels)  median DF=%.2f Hz  RI=%.2f  OI=%.2f', ...
            numel(valid_idx), s_avg.df, s_avg.ri, s_avg.oi));
    end
end


% ========================================================================
function xProc = applyPreproc(x, fs, method)
% Apply preprocessing to x. Works for both column vectors (1-D) and
% matrices (one signal per column, used by the 3-D path).
    switch lower(method)
        case 'none'
            xProc = x;
        case 'detrend'
            xProc = detrend(x, 'linear');   % detrends each column independently
        case 'botteron'
            xProc = botteronSmith(x, fs);   % filtfilt operates column-wise
    end
end


% ========================================================================
function P = bandPower(f, pxx, fc, bw)
%BANDPOWER  Trapezoidal integration of pxx over [fc-bw, fc+bw].
    idx = f >= (fc - bw) & f <= (fc + bw);
    if nnz(idx) < 2
        P = 0;
    else
        P = trapz(f(idx), pxx(idx));
    end
end


% ========================================================================
function y = botteronSmith(x, fs)
%BOTTERONSMITH  Preprocessing chain designed for fractionated atrial
% electrograms: bandpass 40-250 Hz, full-wave rectify, low-pass 20 Hz.
% Accepts a matrix — filtfilt operates on each column independently.

    % Bandpass 40 - min(250, 0.45*fs)
    hi = min(250, 0.45 * fs);
    lo = min(40,  hi - 5);
    [b1, a1] = butter(4, [lo hi] / (fs/2), 'bandpass');
    xb = filtfilt(b1, a1, x);

    xr = abs(xb);                               % full-wave rectify

    % Low-pass 20 Hz
    [b2, a2] = butter(4, 20 / (fs/2), 'low');
    y = filtfilt(b2, a2, xr);
    y = detrend(y, 'linear');
end


% ========================================================================
function plotSpectrum(s)
    figure('Color', 'w', 'Position', [100 100 900 500]);
    plot(s.f, s.pxx, 'k', 'LineWidth', 1); hold on;

    % shade search band
    yl = ylim;
    fill([s.band(1) s.band(2) s.band(2) s.band(1)], ...
         [yl(1) yl(1) yl(2) yl(2)], [.85 .9 .95], ...
         'EdgeColor', 'none', 'FaceAlpha', .4);

    % DF marker + harmonic windows
    for k = 1:s.harmonics
        fk = k * s.df;
        if fk > s.f(end), break; end
        xline(fk, '--', sprintf('%d*DF', k), ...
              'Color', [.8 .2 .2], 'LabelOrientation', 'horizontal');
        lo = fk - s.peakBW;  hi = fk + s.peakBW;
        fill([lo hi hi lo], [yl(1) yl(1) yl(2) yl(2)], ...
             [.95 .8 .8], 'EdgeColor', 'none', 'FaceAlpha', .35);
    end
    plot(s.f, s.pxx, 'k', 'LineWidth', 1);   % redraw on top of fills

    xlim([0, min(s.f(end), 3 * s.band(2))]);
    xlabel('Frequency (Hz)');
    ylabel('PSD');
    title(sprintf(['Welch PSD   DF = %.2f Hz   RI = %.2f   OI = %.2f' ...
                   '   (BW = \\pm%.2f Hz, K = %d)'], ...
                  s.df, s.ri, s.oi, s.peakBW, s.harmonics));
    grid on;
end
