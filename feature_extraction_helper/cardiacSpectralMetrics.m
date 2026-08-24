function [df, ri, oi, spec] = cardiacSpectralMetrics(x, fs, varargin)
%CARDIACSPECTRALMETRICS  Dominant frequency, regularity index and
% organization index of a cardiac action-potential (or electrogram) signal.
%
%   [DF, RI, OI, SPEC] = cardiacSpectralMetrics(X, FS)
%
%   Accepts either:
%     X   [N x 1]            1-D time series  -> scalar DF, RI, OI outputs
%     X   [rows x cols x N]  optical-mapping volume -> rows x cols map outputs
%
%   For 3-D input every pixel is processed independently using a single
%   batched pwelch call (one column per pixel). Pixels that are all-NaN
%   (background mask) are skipped and returned as NaN in the output maps.
%   spec.pxx is returned as (nfft/2+1) x rows x cols; spec.xProc is omitted
%   (storing a preprocessed copy of the full volume is impractical).
%
%   Computes three standard spectral-organization metrics:
%
%     DF  Dominant Frequency (Hz)
%         The frequency that carries the peak power-spectral-density
%         within the DF search band.
%
%     RI  Regularity Index (dimensionless, 0..1)
%         RI = P(DF +/- BW) / P(total-power band)
%         How tightly the signal's power is concentrated at DF. A perfect
%         sinusoid at DF gives RI -> 1, broadband noise gives RI -> 0.
%
%     OI  Organization Index (dimensionless, 0..1), Everett et al. 2001:
%         OI = sum_{k=1..K} P(k*DF +/- BW) / P(total-power band)
%         Shares the fundamental with RI but also counts the harmonic
%         content of a periodic non-sinusoidal signal -- e.g. an action
%         potential train, whose square-ish shape puts substantial power
%         in 2*DF, 3*DF, etc. A highly organized AP train therefore has
%         OI >> RI; a noisy/fibrillatory signal has OI ~= RI ~= low.
%
%   MULTI-SPECIES / VF SUPPORT (adaptive defaults)
%   VF dominant frequency ranges roughly an order of magnitude across
%   species (human/swine ~4-10 Hz ... rat ~15-40 Hz ... mouse ~25-50 Hz),
%   so a single fixed band cannot serve all of them and a fixed denominator
%   makes OI ill-posed (out-of-band harmonics inflate it > 1). By default
%   the three band/width parameters therefore SCALE WITH THE RHYTHM:
%     * 'Band' []            -> auto-detected from the spectrum: the peak of
%                               the (pooled) PSD over a wide 0.5..min(60,0.9*Nyq)
%                               Hz scan defines f0, and the search band is
%                               [max(0.5,0.5*f0), min(0.9*Nyq, 2*f0)].
%                               The 0.5 Hz floor (30 bpm) is below any cardiac
%                               rhythm worth studying; it exists to keep DF off
%                               residual drift, NOT to bound the rhythm. An
%                               earlier 2-3 Hz floor made paced human and
%                               large-animal recordings (1 Hz / 60 bpm)
%                               structurally unmeasurable -- DF locked onto a
%                               harmonic and returned a plausible wrong answer.
%     * 'TotalPowerBand' []  -> DF-relative denominator [0.5, (K+0.5)*DF]
%                               (per pixel), so the counted harmonics always
%                               lie inside it -> OI is bounded 0..1 and
%                               comparable across species.
%     * 'PeakBW' []          -> DF-scaled half-width max(0.5, 0.05*DF), so the
%                               peak window means the same thing at 5 Hz and
%                               45 Hz.
%   Passing any of these explicitly restores the classic fixed behaviour for
%   that parameter (useful to reproduce a specific published band).
%
%   Inputs
%     X   [N x 1] or [rows x cols x N] real-valued signal(s)
%     FS  scalar sampling rate in Hz
%
%   Name / Value options
%     'Band'            [fLo fHi] DF search band        default [] = auto
%                       Pass e.g. [3 12] (AF), [5 20] (large-animal VF) to
%                       force a fixed search band instead of auto-detection.
%     'PeakBW'          half-width (Hz) of peak window  default [] = adaptive
%                       Adaptive = max(0.5, 0.05*DF). Pass a scalar to fix it.
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
%     'TotalPowerBand'  [fLo fHi] for RI/OI denominator default [] = adaptive
%                       Adaptive = [0.5, (Harmonics+0.5)*DF] per pixel. Pass
%                       a fixed [fLo fHi] to normalize to a fixed band.
%     'Plot'            true/false                      default false
%                       For 3-D input plots the spatially-averaged PSD.
%
%   Outputs -- 1-D input
%     DF    dominant frequency (Hz)
%     RI    regularity index (0..1)
%     OI    organization index (0..1)
%     SPEC  struct with fields:
%             .f         frequency vector (Hz)
%             .pxx       PSD estimate (nfft/2+1 x 1)
%             .df        same as DF
%             .ri, .oi   same as RI, OI
%             .peakBW    the (possibly adaptive) half-width used
%             .harmonics .band .totalBand  the resolved settings used
%             .autoBand  true if the search band was auto-detected
%             .adaptiveNorm true if RI/OI used the DF-relative denominator
%             .xProc     preprocessed signal
%
%   Outputs -- 3-D input
%     DF    rows x cols dominant-frequency map (Hz)
%     RI    rows x cols regularity-index map
%     OI    rows x cols organization-index map
%     SPEC  struct with fields:
%             .f         frequency vector (Hz)
%             .pxx       (nfft/2+1) x rows x cols PSD volume
%             .df, .ri, .oi  same as map outputs
%             .harmonics .band  the resolved settings used
%             .peakBW .totalBand  [] when adaptive (they vary per pixel)
%             .autoBand .adaptiveNorm  as above
%             (xProc omitted for 3-D -- too large to store per pixel)
%
%   Example -- 1-D, fully adaptive (works for any species)
%     load('ap_trace.mat');               % x, fs
%     [df, ri, oi] = cardiacSpectralMetrics(x, fs, 'Plot', true);
%
%   Example -- 3-D optical mapping, force a fixed AF band
%     [df_map, ri_map, oi_map, spec] = cardiacSpectralMetrics(data, fs, ...
%                        'Band', [3 12], 'TotalPowerBand', [1 40]);
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
    p.addParameter('Band',           [],   @(v) isempty(v) || (numel(v)==2 && v(2)>v(1)));
    p.addParameter('PeakBW',         [],   @(v) isempty(v) || (isscalar(v) && v > 0));
    p.addParameter('Harmonics',      4,    @(v) v >= 1 && v == round(v));
    p.addParameter('NFFT',           [],   @(v) isempty(v) || (v > 0 && v == round(v)));
    p.addParameter('WindowSec',      2,    @(v) v > 0);
    p.addParameter('OverlapFrac',    0.5,  @(v) v >= 0 && v < 1);
    p.addParameter('Preproc',        'detrend', ...
        @(s) any(strcmpi(s, {'none','detrend','botteron'})));
    p.addParameter('TotalPowerBand', [],   @(v) isempty(v) || (numel(v)==2 && v(2)>v(1)));
    p.addParameter('Plot',           false, @islogical);
    p.parse(x, fs, varargin{:});
    o = p.Results;

    if isempty(o.NFFT), o.NFFT = max(256, 2^nextpow2(4*fs)); end

    % Which behaviours are adaptive (parameter left empty) vs. fixed (given).
    flags.autoBand   = isempty(o.Band);            % auto-detect DF search band
    flags.adaptiveTP = isempty(o.TotalPowerBand);  % DF-relative RI/OI denominator
    flags.adaptiveBW = isempty(o.PeakBW);          % DF-scaled peak half-width

    % ---------- dispatch ----------
    if ~isvector(x) && ndims(double(x)) == 3
        [df, ri, oi, spec] = compute3D(double(x), fs, o, flags);
    else
        [df, ri, oi, spec] = compute1D(double(x(:)), fs, o, flags);
    end
end


% ========================================================================
function band = resolveBand(f, pxxMean, o)
%RESOLVEBAND  Return an explicit search band, or auto-detect one from the PSD.
%   When o.Band is empty the dominant peak of pxxMean over a wide physiological
%   scan sets f0, and the band is [max(0.5,0.5*f0), min(0.9*Nyq, 2*f0)] -- this
%   brackets the rhythm for any species (mouse..human) without clipping.
%
%   The 0.5 Hz floor is deliberate and load-bearing. It is low enough that a
%   1 Hz paced human heart is reachable, and high enough to keep the peak off
%   residual baseline drift. Validated against the pacing rate across paced
%   recordings: at floors of 0.5 and 0.75 Hz a 1 Hz rhythm resolves correctly,
%   at 1.0 Hz and above it does not, and rat paced (5-8 Hz) and rat VF
%   (25.15 Hz) are bit-identical at every floor tested. See
%   manual_validation_helper/mv_paced_df.
    if ~isempty(o.Band)
        band = o.Band;
        return;
    end
    fNyq   = f(end);
    scanLo = 0.5;
    scanHi = min(60, 0.9 * fNyq);
    fallback = [0.5, min(20, 0.9 * fNyq)];
    scan = f >= scanLo & f <= scanHi;
    if nnz(scan) < 2
        band = fallback;
        return;
    end
    fS = f(scan);  pS = pxxMean(scan);
    [~, ip] = max(pS);
    f0 = fS(ip);
    lo = max(0.5, 0.5 * f0);
    hi = min(0.9 * fNyq, 2 * f0);
    if ~(hi > lo)
        band = fallback;
        return;
    end
    band = [lo hi];
end


% ========================================================================
function [bw, totBand] = resolveWidths(df, f, o, flags)
%RESOLVEWIDTHS  Peak half-width and RI/OI denominator band for a given DF.
%   Adaptive (default): bw scales with DF and the total-power band is
%   DF-relative so the counted harmonics stay inside it (OI bounded 0..1).
    if flags.adaptiveBW
        bw = max(0.5, 0.05 * df);
    else
        bw = o.PeakBW;
    end
    if flags.adaptiveTP
        % The denominator must contain the whole peak window, or RI counts power
        % the total does not and comes out > 1. With a fixed 0.5 Hz lower edge
        % that happens whenever bw is comparable to df -- i.e. at slow rhythms,
        % where bw floors at 0.5 Hz while df is only ~1 Hz, putting the peak
        % window's lower edge (df - bw) below 0.5. Track the peak window down.
        totBand = [min(0.5, max(0, df - bw)), min(f(end), (o.Harmonics + 0.5) * df)];
    else
        totBand = o.TotalPowerBand;
    end
end


% ========================================================================
function [ri, oi] = regOrgIndices(f, pxx, df, bw, totBand, K)
%REGORGINDICES  RI and OI from a PSD given DF, peak half-width and denominator.
%   The K harmonic windows [k*df - bw, k*df + bw] are combined as a UNION of
%   frequency bins, not a sum of per-window powers. With bw floored at 0.5 Hz
%   the windows overlap whenever df < 1 Hz (2*bw > df); summing them then
%   counts the shared bins once per window and OI climbs past 1 (-> ~2 at the
%   0.5 Hz scan floor). Integrating each bin at most once keeps OI <= 1 by
%   construction for any df.
    inTot  = f >= totBand(1) & f <= totBand(2);
    totalP = maskPower(f, pxx, inTot);
    if ~(totalP > 0)
        ri = NaN;  oi = NaN;  return;
    end
    % RI numerator clipped to the denominator band, exactly like OI below:
    % unclipped, a bw-floored peak window that extends past a fixed
    % TotalPowerBand edge counts power the denominator excludes, so RI can
    % exceed OI (and 1). Sharing maskPower also keeps RI and OI on one
    % integrator by construction.
    inPeak = f >= df - bw & f <= df + bw;
    ri = maskPower(f, pxx, inPeak & inTot) / totalP;
    % No early break for windows crossing totBand(2): the inTot clip below
    % keeps the integration inside the denominator band, so a partially
    % in-band harmonic contributes its in-band part instead of being dropped
    % whole (dropping it stepped OI discontinuously at df = 2*bw, i.e. the
    % df = 1 Hz contour with the 0.5 Hz bw floor).
    inHarm = false(size(f));
    for k = 1:K
        fk = k * df;
        inHarm = inHarm | (f >= fk - bw & f <= fk + bw);
    end
    oi = maskPower(f, pxx, inHarm & inTot) / totalP;
end


% ========================================================================
function P = maskPower(f, pxx, idx)
%MASKPOWER  Trapezoidal power over a logical bin mask, integrated run by run.
%   trapz over a non-contiguous index set would bridge the gaps between runs
%   and count power that lies outside the mask, so each contiguous run of
%   true bins is integrated separately and the results summed.
    P = 0;
    d = diff([false; idx(:); false]);
    runStart = find(d == 1);
    runEnd   = find(d == -1) - 1;
    for r = 1:numel(runStart)
        ii = runStart(r):runEnd(r);
        if numel(ii) >= 2
            P = P + trapz(f(ii), pxx(ii));
        end
    end
end


% ========================================================================
function [df, ri, oi, spec] = compute1D(x, fs, o, flags)
% Single-channel path.

    % --- preprocessing ---
    xProc = applyPreproc(x, fs, o.Preproc);

    % --- Welch PSD ---
    winN  = min(numel(xProc), max(64, round(o.WindowSec * fs)));
    nOvlp = round(o.OverlapFrac * winN);
    [pxx, f] = pwelch(xProc, hann(winN), nOvlp, o.NFFT, fs);

    % --- resolve search band (explicit or auto-detected) ---
    band   = resolveBand(f, pxx, o);
    inBand = f >= band(1) & f <= band(2);
    assert(any(inBand), 'Search band [%g %g] Hz contains no FFT bins.', ...
           band(1), band(2));
    fBand   = f(inBand);
    pxxBand = pxx(inBand);
    [~, iPk] = max(pxxBand);
    df = fBand(iPk);

    % --- sanity guards -------------------------------------------------------
    % These are about TRUSTWORTHINESS, deliberately kept separate from the band
    % floor, which is about REACHABILITY. Folding the two together (e.g. raising
    % the floor on short records) reintroduces the very failure the 0.5 Hz floor
    % fixes. Warn instead, and let the caller decide.
    %
    % Only compute1D warns. compute3D runs the same maths per pixel and would
    % emit thousands of identical messages.
    if flags.autoBand && (iPk == 1 || iPk == numel(pxxBand))
        warning('cardiacSpectralMetrics:peakAtBandEdge', ...
            ['DF %.3g Hz sits at the edge of the auto search band [%.3g %.3g] Hz — ' ...
             'the true peak may lie outside it. Pass an explicit ''Band'' to check.'], ...
            df, band(1), band(2));
    end
    nCycles = df * numel(x) / fs;
    if nCycles < 3
        warning('cardiacSpectralMetrics:fewCycles', ...
            ['DF %.3g Hz over %.3g s is only %.1f cycles — too few for a stable ' ...
             'estimate. Use a longer recording.'], df, numel(x) / fs, nCycles);
    end

    % --- widths + RI/OI ---
    [bw, totBand] = resolveWidths(df, f, o, flags);
    [ri, oi]      = regOrgIndices(f, pxx, df, bw, totBand, o.Harmonics);

    % --- pack struct ---
    spec = struct( ...
        'f',            f, ...
        'pxx',          pxx, ...
        'df',           df, ...
        'ri',           ri, ...
        'oi',           oi, ...
        'peakBW',       bw, ...
        'harmonics',    o.Harmonics, ...
        'band',         band, ...
        'totalBand',    totBand, ...
        'autoBand',     flags.autoBand, ...
        'adaptiveNorm', flags.adaptiveTP, ...
        'xProc',        xProc);

    if o.Plot
        plotSpectrum(spec);
    end
end


% ========================================================================
function [df_map, ri_map, oi_map, spec] = compute3D(x, fs, o, flags)
% Batched 3-D path: single pwelch call across all valid pixels.
%
% x is rows x cols x frames.
% Valid pixels are those with no NaN samples (background mask is all-NaN).

    [num_rows, num_cols, num_frames] = size(x);
    num_pixels = num_rows * num_cols;

    % Reshape to frames x pixels -- pwelch expects samples along rows
    x2d = reshape(permute(x, [3 1 2]), num_frames, num_pixels);

    % Pixels with any NaN frame are excluded (background mask)
    valid_mask = ~any(isnan(x2d), 1);   % 1 x num_pixels logical

    % --- preprocess valid pixels only ---
    x2d_valid = applyPreproc(x2d(:, valid_mask), fs, o.Preproc);

    % --- single pwelch call for all valid pixels ---
    winN  = min(num_frames, max(64, round(o.WindowSec * fs)));
    nOvlp = round(o.OverlapFrac * winN);
    [pxx_valid, f] = pwelch(x2d_valid, hann(winN), nOvlp, o.NFFT, fs);

    nF  = numel(f);
    pxx = NaN(nF, num_pixels);
    pxx(:, valid_mask) = pxx_valid;

    % --- resolve search band from the pooled (mean valid) PSD ---
    mean_pxx = mean(pxx_valid, 2, 'omitnan');
    band     = resolveBand(f, mean_pxx, o);
    inBand   = f >= band(1) & f <= band(2);
    assert(any(inBand), 'Search band [%g %g] Hz contains no FFT bins.', ...
           band(1), band(2));
    fBand    = f(inBand);

    % --- dominant frequency -- vectorised max across all pixels at once ---
    [~, iPk] = max(pxx(inBand, :), [], 1);     % 1 x num_pixels
    df_vec   = fBand(iPk);                      % 1 x num_pixels
    df_vec(~valid_mask) = NaN;

    % --- RI and OI -- per-pixel (df sets the windows; loop is cheap vs pwelch) ---
    ri_vec = NaN(1, num_pixels);
    oi_vec = NaN(1, num_pixels);
    valid_idx = find(valid_mask);
    for vi = 1:numel(valid_idx)
        px = valid_idx(vi);
        df_px = df_vec(px);
        [bw, totBand] = resolveWidths(df_px, f, o, flags);
        [ri_vec(px), oi_vec(px)] = ...
            regOrgIndices(f, pxx(:, px), df_px, bw, totBand, o.Harmonics);
    end

    % --- reshape vectors back to spatial maps ---
    df_map = reshape(df_vec, num_rows, num_cols);
    ri_map = reshape(ri_vec, num_rows, num_cols);
    oi_map = reshape(oi_vec, num_rows, num_cols);

    % --- pack struct (peakBW/totalBand vary per pixel when adaptive -> []) ---
    if flags.adaptiveBW, peakBW_out = [];        else, peakBW_out = o.PeakBW;        end
    if flags.adaptiveTP, totBand_out = [];       else, totBand_out = o.TotalPowerBand; end
    spec = struct( ...
        'f',            f, ...
        'pxx',          reshape(pxx, nF, num_rows, num_cols), ...
        'df',           df_map, ...
        'ri',           ri_map, ...
        'oi',           oi_map, ...
        'peakBW',       peakBW_out, ...
        'harmonics',    o.Harmonics, ...
        'band',         band, ...
        'totalBand',    totBand_out, ...
        'autoBand',     flags.autoBand, ...
        'adaptiveNorm', flags.adaptiveTP);

    % --- optional plot: spatially-averaged PSD ---
    if o.Plot
        med_df = median(df_vec(valid_mask), 'omitnan');
        [rep_bw, rep_tot] = resolveWidths(med_df, f, o, flags);
        s_avg = struct( ...
            'f',         f, ...
            'pxx',       mean_pxx, ...
            'df',        med_df, ...
            'ri',        median(ri_vec(valid_mask), 'omitnan'), ...
            'oi',        median(oi_vec(valid_mask), 'omitnan'), ...
            'peakBW',    rep_bw, ...
            'harmonics', o.Harmonics, ...
            'band',      band, ...
            'totalBand', rep_tot);
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
function y = botteronSmith(x, fs)
%BOTTERONSMITH  Preprocessing chain designed for fractionated atrial
% electrograms: bandpass 40-250 Hz, full-wave rectify, low-pass 20 Hz.
% Accepts a matrix -- filtfilt operates on each column independently.

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
