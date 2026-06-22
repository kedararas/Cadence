function [alternans] = analyzeCaTransientAlternans(time, cSignal, varargin)
% ANALYZECATRANSIENTALTERNANS  Calcium transient alternans — 1-D or 3-D.
%
%   ── 1-D mode (single trace) ─────────────────────────────────────────────
%   alternans = analyzeCaTransientAlternans(time, cSignal)
%   alternans = analyzeCaTransientAlternans(time, cSignal, Name, Value, ...)
%
%   Inputs
%     time    - Time vector (ms), length T
%     cSignal - Calcium signal (F/F0, ratio, or raw fluorescence), length T
%
%   ── 3-D mode (full optical-mapping stack) ────────────────────────────────
%   alternans = analyzeCaTransientAlternans(time, cData3D, 'BeatFrames', bf, ...)
%
%   Inputs
%     time    - Time vector (ms), length T (number of frames)
%     cData3D - [rows × cols × T] fluorescence array (any numeric class)
%     BeatFrames - [num_beats × 2] [start_frame, end_frame] per beat (REQUIRED)
%
%   Optional Name-Value Pairs (both modes)
%     'BaselineWindow' - Duration (ms) for sliding-minimum baseline (default: 50)
%     'PlotResults'    - true/false (default: true for 1-D, false for 3-D)
%
%   Optional Name-Value Pairs (1-D mode only)
%     'PacingRate'     - Stimulation rate Hz (default: auto-detect)
%     'Threshold'      - Peak detection threshold, fraction of max (default: 0.3)
%
%   Optional Name-Value Pairs (3-D mode only)
%     'BeatFrames'     - [num_beats × 2] start/end frame index of each beat
%     'Mask'           - [rows × cols] logical/numeric; pixels outside = NaN
%     'DecayLevels'    - Decay % levels for CaTD maps (default: [50 80])
%                        D50 = time from peak to 50% decay; D80 = 80% decay
%
%   3-D output struct fields
%     .mode           '3D'
%     .amp_alt_map    Amplitude alternans map
%     .D50_alt_map    CaTD50 alternans map (ms)
%     .D80_alt_map    CaTD80 alternans map (ms)
%     .TTP_alt_map    Time-to-peak alternans map (ms)
%     .pval_map       Pixel-wise paired t-test p-value (amplitude)
%     .tstat_map      Pixel-wise t-statistic
%     .spectral_map   Pixel-wise spectral alternans index (power at 0.5 cyc/beat)
%     .amp_beat       {num_beats × 1} cell of per-beat amplitude maps
%     .CaTD_beat      {num_beats × num_levels} cell of per-beat CaTD maps (ms)
%     .TTP_beat       {num_beats × 1} cell of per-beat TTP maps (ms)

    %% ---- Parse inputs ----
    p = inputParser;
    addRequired(p,  'time',          @isvector);
    addRequired(p,  'cSignal',       @isnumeric);   % vector OR 3-D array
    addParameter(p, 'PacingRate',    [],    @isnumeric);
    addParameter(p, 'Threshold',     0.3,   @isnumeric);
    addParameter(p, 'BaselineWindow',50,    @isnumeric);
    addParameter(p, 'PlotResults',   [],    @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'BeatFrames',    [],    @isnumeric);
    addParameter(p, 'Mask',          [],    @(x) isnumeric(x) || islogical(x));
    addParameter(p, 'DecayLevels',   [50 80], @isvector);
    addParameter(p, 'SignalType',    'CaTransient', @(x) ischar(x) || isstring(x));
    % 'CaTransient' : intracellular Ca2+, signal rises during systole (default)
    % 'SR'          : intra-SR Ca2+ (e.g. Fluo-5N), signal falls during release
    parse(p, time, cSignal, varargin{:});
    opts = p.Results;

    opts.SignalType = char(opts.SignalType);

    %% ---- Route: 3-D spatial mode ----------------------------------------
    if ndims(cSignal) == 3
        if isempty(opts.BeatFrames)
            error('analyzeCaTransientAlternans:noBeatFrames', ...
                ['BeatFrames is required for 3-D input.\n' ...
                 'Supply a [num_beats x 2] array of [start_frame, end_frame].']);
        end
        if isempty(opts.PlotResults), opts.PlotResults = false; end
        alternans = caAlternans3D(time, cSignal, opts);
        return;
    end

    %% ---- 1-D single-trace mode ------------------------------------------
    if isempty(opts.PlotResults), opts.PlotResults = true; end

    time    = time(:);
    cSignal = cSignal(:);
    dt      = mean(diff(time));

    fprintf('\n======== Ca2+ Transient Alternans Analysis ========\n');

    %% ---- Step 1: Baseline correction ----
    baselineSamples = round(opts.BaselineWindow / dt);
    baseline        = movmin(cSignal, baselineSamples);
    cNorm           = cSignal - baseline;
    cNorm           = cNorm / max(cNorm);

    %% ---- Step 2: Detect peaks ----
    minPeakHeight = opts.Threshold * max(cNorm);
    if ~isempty(opts.PacingRate)
        minPeakDist = round((1000 / opts.PacingRate) * 0.6 / dt);
    else
        minPeakDist = round(0.1 * length(cNorm));
    end

    [peakAmps, peakLocs] = findpeaks(cNorm, ...
        'MinPeakHeight',     minPeakHeight, ...
        'MinPeakDistance',   minPeakDist,   ...
        'MinPeakProminence', 0.05);

    nBeats = length(peakLocs);
    fprintf('Detected beats     : %d\n', nBeats);

    if nBeats < 4
        warning('Too few beats (%d) for reliable alternans analysis.', nBeats);
        alternans = struct(); return;
    end

    %% ---- Step 3: Onset detection ----
    onsetLocs = zeros(nBeats, 1);
    for b = 1:nBeats
        searchStart = max(1, peakLocs(b) - minPeakDist);
        searchSeg   = cNorm(searchStart:peakLocs(b));
        [~, relMin] = min(searchSeg);
        onsetLocs(b) = searchStart + relMin - 1;
    end

    %% ---- Step 4: Per-beat metrics ----
    amplitude      = zeros(nBeats, 1);
    TTP            = zeros(nBeats, 1);
    D50            = zeros(nBeats, 1);
    D80            = zeros(nBeats, 1);
    tau            = zeros(nBeats, 1);
    diastolicLevel = zeros(nBeats, 1);   % pre-release baseline per beat

    for b = 1:nBeats
        baseVal        = cNorm(onsetLocs(b));
        diastolicLevel(b) = baseVal;     % diastolic [Ca2+] before each transient
        amplitude(b)   = peakAmps(b) - baseVal;
        TTP(b)       = (peakLocs(b) - onsetLocs(b)) * dt;

        if b < nBeats
            decayEnd = onsetLocs(b+1);
        else
            decayEnd = min(peakLocs(b) + round(500/dt), length(cNorm));
        end
        decaySeg  = cNorm(peakLocs(b):decayEnd);
        decayTime = (0:length(decaySeg)-1)' * dt;

        thresh50 = peakAmps(b) - 0.50 * amplitude(b);
        thresh80 = peakAmps(b) - 0.80 * amplitude(b);
        idx50    = find(decaySeg <= thresh50, 1, 'first');
        idx80    = find(decaySeg <= thresh80, 1, 'first');
        D50(b)   = idx50 * dt;
        D80(b)   = idx80 * dt;

        try
            fitObj  = fit(decayTime, decaySeg, 'a*exp(-x/b)+c', ...
                          'StartPoint', [amplitude(b), 100, baseVal], ...
                          'Lower', [0, 1, -0.5], 'Upper', [2, 2000, 0.5]);
            tau(b)  = fitObj.b;
        catch
            tau(b) = NaN;
        end
    end

    %% ---- Step 5: Alternans metrics ----
    DI = zeros(nBeats-1, 1);
    for b = 1:nBeats-1
        DI(b) = (amplitude(b+1) - amplitude(b)) / ...
                (amplitude(b+1) + amplitude(b) + eps);
    end

    evenBeats = amplitude(2:2:end);
    oddBeats  = amplitude(1:2:end);
    nPairs    = min(length(evenBeats), length(oddBeats));
    ALT       = mean(abs(evenBeats(1:nPairs) - oddBeats(1:nPairs)));
    ALT_ratio = ALT / mean(amplitude);

    if nPairs >= 3
        [~, pValue] = ttest(evenBeats(1:nPairs), oddBeats(1:nPairs));
    else
        pValue = NaN;
    end

    [spectralALT, f_alt] = spectralAlternansIndex1D(amplitude);
    isAlternans = (ALT_ratio > 0.10) && (~isnan(pValue) && pValue < 0.05);

    D50_even = D50(2:2:end); D50_odd = D50(1:2:end);
    D80_even = D80(2:2:end); D80_odd = D80(1:2:end);
    nD = min(length(D50_even), length(D50_odd));
    ALT_D50 = mean(abs(D50_even(1:nD) - D50_odd(1:nD)));
    ALT_D80 = mean(abs(D80_even(1:nD) - D80_odd(1:nD)));
    ALT_ratio_D50 = ALT_D50 / (mean(D50, 'omitnan') + eps);
    ALT_ratio_D80 = ALT_D80 / (mean(D80, 'omitnan') + eps);

    % Phase: +1 if odd beats have larger amplitude than even, -1 otherwise
    alt_phase = sign(mean(oddBeats(1:nPairs) - evenBeats(1:nPairs)));

    %% ---- Step 6: Wang load and release alternans (Wang et al. 2014) ----
    % Release alternans = 1 - S/L
    %   L = mean amplitude of large beats, S = mean amplitude of small beats
    % Load alternans = D/L
    %   D = mean |diastolicLevel_L - diastolicLevel_S|, normalised by L
    [releaseALT_Wang, loadALT_Wang, wang_L_mean, wang_S_mean, wang_L_idx, wang_S_idx] = ...
        wangLoadRelease(amplitude, diastolicLevel);

    %% ---- Assemble output ----
    alternans.mode             = '1D';
    alternans.time             = time;
    alternans.signal           = cNorm;
    alternans.peakLocs         = peakLocs;
    alternans.onsetLocs        = onsetLocs;
    alternans.amplitude        = amplitude;
    alternans.diastolicLevel   = diastolicLevel;
    alternans.TTP              = TTP;
    alternans.D50              = D50;
    alternans.D80              = D80;
    alternans.tau              = tau;
    alternans.DI               = DI;
    alternans.ALT              = ALT;
    alternans.ALT_ratio        = ALT_ratio;
    alternans.ALT_D50          = ALT_D50;
    alternans.ALT_D80          = ALT_D80;
    alternans.ALT_ratio_D50    = ALT_ratio_D50;   % |ALT_D50| / mean_D50
    alternans.ALT_ratio_D80    = ALT_ratio_D80;   % |ALT_D80| / mean_D80
    alternans.alt_phase        = alt_phase;        % +1 odd>even, -1 odd<even
    alternans.spectralALT      = spectralALT;
    alternans.f_alt            = f_alt;
    alternans.pValue           = pValue;
    alternans.isAlternans      = isAlternans;
    alternans.nBeats           = nBeats;
    % Wang et al. 2014 load and release alternans metrics
    alternans.releaseALT_Wang  = releaseALT_Wang;  % 1 - S/L
    alternans.loadALT_Wang     = loadALT_Wang;     % D/L
    alternans.wang_L_mean      = wang_L_mean;      % mean large-beat amplitude
    alternans.wang_S_mean      = wang_S_mean;      % mean small-beat amplitude
    alternans.wang_L_idx       = wang_L_idx;       % indices of large beats
    alternans.wang_S_idx       = wang_S_idx;       % indices of small beats

    fprintf('Mean amplitude     : %.4f +/- %.4f\n', mean(amplitude), std(amplitude));
    fprintf('ALT magnitude      : %.4f\n', ALT);
    fprintf('ALT ratio          : %.2f%%\n', ALT_ratio*100);
    fprintf('Spectral ALT index : %.4f  (at %.4f cycles/beat)\n', spectralALT, f_alt);
    fprintf('p-value (t-test)   : %.4f\n', pValue);
    fprintf('\n--- Wang et al. 2014 SR Ca2+ Alternans Metrics ---\n');
    fprintf('Release ALT (1-S/L): %.4f  (L=%.4f, S=%.4f)\n', ...
            releaseALT_Wang, wang_L_mean, wang_S_mean);
    fprintf('Load ALT    (D/L)  : %.4f\n', loadALT_Wang);
    if isAlternans
        fprintf('>> ALTERNANS DETECTED (ALT=%.2f%%, p=%.3f)\n', ALT_ratio*100, pValue);
    else
        fprintf('>> No significant alternans detected.\n');
    end
    fprintf('====================================================\n\n');

    if opts.PlotResults
        plotCaAlternans(alternans);
    end
end


% =========================================================================
function alternans = caAlternans3D(time, cData3D, opts)
% CAALTERNANS3D  Vectorized Ca2+ transient alternans across [R x C x T].
%
%   Pipeline per beat window
%     1. Sliding-minimum baseline correction (movmin, 3-D native)
%     2. Per-pixel [0,1] normalisation
%     3. Peak location via max along time
%     4. Onset: minimum in the first half of the beat window
%     5. Amplitude = peak - onset value (on baseline-corrected signal)
%     6. TTP = (peak_frame - onset_frame) * dt
%     7. CaTD at each decay level via vectorised cumsum crossing
%
%   Alternans magnitude = mean(odd - even) over all beat pairs.
%   Statistics: vectorised paired t-test on amplitude (requires ≥ 3 pairs).
%   Spectral:   FFT of per-beat amplitude series along beat dimension.

    dt          = mean(diff(time(:)));
    [R, C, ~]   = size(cData3D);
    bf          = opts.BeatFrames;          % [num_beats x 2]
    nBeats      = size(bf, 1);
    decayLevels = opts.DecayLevels(:)';     % row vector, e.g. [50 80]
    nLev        = numel(decayLevels);
    blWin       = max(3, round(opts.BaselineWindow / dt));

    % ── Mask ───────────────────────────────────────────────────────────────
    if ~isempty(opts.Mask)
        %nan_mask = double(logical(opts.Mask));
        %nan_mask(nan_mask == 0) = nan;
        nan_mask = opts.Mask;
    else
        nan_mask = ones(R, C);
    end

    % ── Global baseline correction (movmin along time, 3-D native) ─────────
    % movmin(A, k, 3) slides a k-frame minimum window along dim 3 of the
    % entire [R x C x T] array in one call — no pixel loop needed.
    cDouble   = double(cData3D);
    baseline3D = movmin(cDouble, blWin, 3);        % R x C x T
    cCorr      = cDouble - baseline3D;             % baseline-subtracted

    % ── Per-beat storage ───────────────────────────────────────────────────
    amp_beat      = cell(nBeats, 1);       % amplitude (baseline-subtracted, norm.)
    TTP_beat      = cell(nBeats, 1);       % time to peak (ms)
    CaTD_beat     = cell(nBeats, nLev);    % CaTD at each decay level (ms)
    diastolic_beat = cell(nBeats, 1);      % diastolic level per beat (pre-release baseline)

    for j = 1:nBeats
        sf    = bf(j, 1);
        ef    = bf(j, 2);
        T_win = ef - sf + 1;

        win = cCorr(:,:, sf:ef);           % R x C x T_win, baseline-subtracted

        % Per-pixel [0,1] normalisation
        lo  = min(win, [], 3);
        hi  = max(win, [], 3);
        rng = max(hi - lo, eps('single'));
        wn  = (win - lo) ./ rng;           % R x C x T_win, each pixel in [0,1]

        % Amplitude (raw, before normalisation — in baseline-subtracted units)
        amp_beat{j} = (hi - lo) .* nan_mask;

        % Diastolic level: the pre-release baseline (lo in the beat window,
        % i.e. the minimum of the baseline-corrected signal = onset level).
        % For SR Ca2+ (Fluo-5N, signal falls on release), use hi instead;
        % the SignalType field is checked below after the loop.
        diastolic_beat{j} = lo .* nan_mask;

        % Onset: minimum in the first half of the window (per pixel)
        T_half = max(1, floor(T_win / 2));
        [~, onset_f] = min(wn(:,:, 1:T_half), [], 3);   % R x C, frame in [1,T_half]

        % Peak: maximum across full window
        [~, pk] = max(wn, [], 3);                        % R x C

        % TTP (ms): from onset to peak
        TTP_beat{j} = max(double(pk) - double(onset_f), 0) * dt .* nan_mask;

        % CaTD: first frame AFTER peak where signal <= (1 - decay_pct/100)
        t_idx    = reshape(1:T_win, 1, 1, T_win);
        after_pk = t_idx > reshape(pk, R, C, 1);        % R x C x T_win gate

        for lv = 1:nLev
            thresh        = 1 - decayLevels(lv) / 100;  % e.g. 0.50 for D50, 0.20 for D80
            below         = wn <= thresh;
            [hit, rf]     = max(cumsum(below & after_pk, 3) == 1, [], 3);
            catd          = (double(rf) - double(pk)) * dt;  % ms from peak to crossing
            invalid       = ~logical(hit) | catd <= 0 | isnan(TTP_beat{j});
            catd(invalid) = nan;
            CaTD_beat{j, lv} = catd;
        end
    end

    % ── SR Ca2+ mode: diastolic = pre-release peak (hi), not onset (lo) ───
    % For Fluo-5N the signal is high at diastole and falls on release, so the
    % diastolic level within each beat window is the maximum (hi), not lo.
    if strcmpi(opts.SignalType, 'SR')
        for j = 1:nBeats
            sf  = bf(j, 1);
            ef  = bf(j, 2);
            win = cCorr(:,:, sf:ef);
            hi_sr = max(win, [], 3);
            diastolic_beat{j} = hi_sr .* nan_mask;
        end
    end

    % ── Per-level reachability diagnostic ──────────────────────────────────
    % Fraction of (in-mask pixel × beat) entries that yield a finite CaTD at
    % each decay level.  A low value at the deepest level means the transient
    % does not decay that far within the beat window (common at fast pacing)
    % -> the corresponding alternans map will be largely NaN.
    if islogical(nan_mask) || all(ismember(nan_mask(isfinite(nan_mask)), [0 1]))
        tissue = isfinite(nan_mask) & nan_mask ~= 0;
    else
        tissue = isfinite(nan_mask);
    end
    nT = max(nnz(tissue), 1);
    level_coverage = zeros(1, nLev);
    for lv = 1:nLev
        cnt = 0;
        for j = 1:nBeats
            m   = CaTD_beat{j, lv};
            cnt = cnt + nnz(isfinite(m) & tissue);
        end
        level_coverage(lv) = cnt / (nT * nBeats);
    end
    fprintf('Ca alternans reachability (tissue-beats with measurable CaTD):\n');
    for lv = 1:nLev
        fprintf('  CaTD%-2d : %5.1f%%\n', decayLevels(lv), 100*level_coverage(lv));
    end
    if level_coverage(end) < 0.5
        warning('analyzeCaTransientAlternans:lowCoverage', ...
            ['CaTD%d reached in only %.0f%% of tissue-beats; its alternans map ' ...
             'will be largely NaN. The transient likely does not decay to %d%% ' ...
             'within the beat window at this pacing rate -- consider a lower ' ...
             'decay level.'], decayLevels(end), 100*level_coverage(end), decayLevels(end));
    end

    % ── Odd / even split ──────────────────────────────────────────────────
    odd_idx  = 1:2:nBeats;
    even_idx = 2:2:nBeats;
    n_pairs  = min(numel(odd_idx), numel(even_idx));

    if n_pairs < 1
        warning('analyzeCaTransientAlternans:toofewbeats', ...
            'Need at least 2 beats. Returning empty struct.');
        alternans = struct('mode','3D'); return;
    end

    % ── Amplitude alternans ────────────────────────────────────────────────
    odd_amp  = cat(3, amp_beat{odd_idx(1:n_pairs)});   % R x C x n_pairs
    even_amp = cat(3, amp_beat{even_idx(1:n_pairs)});
    amp_alt  = mean(odd_amp - even_amp, 3, 'omitnan');

    % ── CaTD alternans at each decay level ────────────────────────────────
    CaTD_alt_maps = cell(nLev, 1);
    for lv = 1:nLev
        odd_s  = cat(3, CaTD_beat{odd_idx(1:n_pairs),  lv});
        even_s = cat(3, CaTD_beat{even_idx(1:n_pairs), lv});
        CaTD_alt_maps{lv} = mean(odd_s - even_s, 3, 'omitnan');
    end

    % ── TTP alternans ─────────────────────────────────────────────────────
    odd_ttp  = cat(3, TTP_beat{odd_idx(1:n_pairs)});
    even_ttp = cat(3, TTP_beat{even_idx(1:n_pairs)});
    TTP_alt  = mean(odd_ttp - even_ttp, 3, 'omitnan');

    % ── Wang et al. 2014: release and load alternans maps (per-pixel L/S) ──
    % Per-pixel classification: at each pixel independently, the group with
    % the larger mean amplitude is L and the smaller is S.  This is correct
    % for spatially discordant alternans where the alternans phase inverts
    % across the tissue — a global classification would swap L/S in the
    % minority-phase region and produce spurious negative release_alt values.
    % With per-pixel classification, L >= S everywhere so release_alt_map is
    % guaranteed to lie in [0, 1] and load_alt_map >= 0.
    %
    % odd_larger: [R x C] logical — true where odd beats are locally larger.
    % This is the per-pixel amplitude phase, analogous to phase_map for APD.
    odd_mean_map  = mean(odd_amp,  3, 'omitnan');   % R x C
    even_mean_map = mean(even_amp, 3, 'omitnan');   % R x C
    odd_larger    = odd_mean_map >= even_mean_map;  % R x C logical

    odd_diast_stack  = cat(3, diastolic_beat{odd_idx(1:n_pairs)});
    even_diast_stack = cat(3, diastolic_beat{even_idx(1:n_pairs)});

    % Initialise L stacks as odd, then swap pixels where even is locally larger
    L_amp_stack   = odd_amp;
    S_amp_stack   = even_amp;
    L_diast_stack = odd_diast_stack;
    S_diast_stack = even_diast_stack;

    swap = repmat(~odd_larger, 1, 1, n_pairs);   % R x C x n_pairs swap mask
    L_amp_stack(swap)   = even_amp(swap);
    S_amp_stack(swap)   = odd_amp(swap);
    L_diast_stack(swap) = even_diast_stack(swap);
    S_diast_stack(swap) = odd_diast_stack(swap);

    L_mean_map = mean(L_amp_stack, 3, 'omitnan');   % R x C, always >= S_mean_map
    S_mean_map = mean(S_amp_stack, 3, 'omitnan');   % R x C

    % Release alternans map: 1 - S/L  (Wang et al. 2014) — result in [0, 1]
    release_alt_map = 1 - S_mean_map ./ (L_mean_map + eps);
    % Load alternans map: D/L  (Wang et al. 2014) — result >= 0
    D_map        = mean(abs(L_diast_stack - S_diast_stack), 3, 'omitnan');
    load_alt_map = D_map ./ (L_mean_map + eps);

    % ── Paired t-test map on amplitude (vectorised, requires >= 3 pairs) ──
    pval_map  = nan(R, C);
    tstat_map = nan(R, C);
    if n_pairs >= 3
        d       = odd_amp - even_amp;              % R x C x n_pairs
        d_mean  = mean(d, 3, 'omitnan');
        d_std   = std(d,  0, 3, 'omitnan');
        tstat   = d_mean ./ (d_std / sqrt(n_pairs));
        tstat_map = tstat;
        try
            pval_map = 2 * (1 - tcdf(abs(tstat), n_pairs - 1));
        catch
            pval_map = 2 * (1 - normcdf(abs(tstat)));
        end
    end

    % ── Spectral alternans index map ──────────────────────────────────────
    % Power at 0.5 cycles/beat in the per-beat amplitude series.
    spectral_map = nan(R, C);
    if nBeats >= 4
        all_amp = cat(3, amp_beat{:});              % R x C x nBeats
        all_amp(isnan(all_amp)) = 0;
        all_amp = all_amp - mean(all_amp, 3);       % zero-mean detrend
        Y        = fft(all_amp, [], 3);
        P        = abs(Y / nBeats).^2;
        alt_bin  = floor(nBeats/2) + 1;
        if alt_bin <= size(P, 3)
            spectral_map = P(:,:, alt_bin);
        end
    end

    % ── Alternans ratio maps ───────────────────────────────────────────────
    amp_mean_map  = mean(cat(3, amp_beat{:}), 3, 'omitnan');
    amp_ratio_map = abs(amp_alt) ./ (amp_mean_map + eps);

    CaTD_ratio_maps = cell(nLev, 1);
    for lv = 1:nLev
        catd_mean          = mean(cat(3, CaTD_beat{:, lv}), 3, 'omitnan');
        CaTD_ratio_maps{lv} = abs(CaTD_alt_maps{lv}) ./ (catd_mean + eps);
    end

    TTP_mean_map  = mean(cat(3, TTP_beat{:}), 3, 'omitnan');
    TTP_ratio_map = abs(TTP_alt) ./ (TTP_mean_map + eps);

    % ── Phase maps ─────────────────────────────────────────────────────────
    % Binary phase map: +1 where odd-beat Ca amplitude > even-beat,
    %                   -1 where odd < even, 0 where no data.
    phase_map = sign(amp_alt);

    % Continuous phase angle (radians, −π to +π) from FFT at 0.5 cyc/beat.
    % Gradual spatial phase gradients identify travelling Ca alternans waves.
    % Anti-phase pixels (angle ≈ ±π apart) are discordant.
    phase_angle_map = nan(R, C);
    if nBeats >= 4
        all_amp = cat(3, amp_beat{:});       % R x C x nBeats
        all_amp(isnan(all_amp)) = 0;
        all_amp = all_amp - mean(all_amp, 3);
        Y_amp   = fft(all_amp, [], 3);
        alt_bin = floor(nBeats/2) + 1;
        if alt_bin <= size(Y_amp, 3)
            phase_angle_map = angle(Y_amp(:,:, alt_bin));  % [-pi, pi]
        end
    end

    % ── Package output ─────────────────────────────────────────────────────
    alternans.mode         = '3D';
    alternans.num_beats    = nBeats;
    alternans.n_pairs      = n_pairs;
    alternans.DecayLevels  = decayLevels;
    alternans.dt           = dt;
    alternans.level_coverage = level_coverage;   % fraction of tissue-beats reachable per decay level

    % Primary alternans maps
    alternans.amp_alt_map  = amp_alt;
    alternans.TTP_alt_map  = TTP_alt;

    % Named CaTD alternans maps (e.g. alternans.D50_alt_map, alternans.D80_alt_map)
    for lv = 1:nLev
        alternans.(sprintf('D%d_alt_map', decayLevels(lv))) = CaTD_alt_maps{lv};
    end
    alternans.CaTD_alt_maps = CaTD_alt_maps;  % cell for programmatic access

    % Statistical and spectral maps
    alternans.pval_map      = pval_map;
    alternans.tstat_map     = tstat_map;
    alternans.spectral_map  = spectral_map;

    % Alternans ratio maps  (|ALT| / mean_metric per pixel)
    alternans.amp_ratio_map  = amp_ratio_map;
    alternans.TTP_ratio_map  = TTP_ratio_map;
    for lv = 1:nLev
        alternans.(sprintf('D%d_ratio_map', decayLevels(lv))) = CaTD_ratio_maps{lv};
    end
    alternans.CaTD_ratio_maps = CaTD_ratio_maps;  % cell for programmatic access

    % Phase maps
    alternans.phase_map       = phase_map;        % Ca amplitude phase: +1 odd>even, -1 even>odd, 0 no data
    %   Direct analogue of substrate.phase_map (APD-based, from assess_arrhythmia_substrate).
    %   Product  phase_map .* substrate.phase_map  = +1 (Ca/AP in-phase) or -1 (out-of-phase).
    alternans.phase_angle_map = phase_angle_map;  % continuous  [-pi, pi] rad

    % Wang et al. 2014 load and release alternans maps (per-pixel L/S)
    alternans.release_alt_map = release_alt_map;  % 1 - S/L per pixel, always in [0,1]
    alternans.load_alt_map    = load_alt_map;     % D/L   per pixel, always >= 0
    alternans.L_mean_map      = L_mean_map;       % mean large-beat amplitude map
    alternans.S_mean_map      = S_mean_map;       % mean small-beat amplitude map
    %   phase_map (line 528) is the Ca-amplitude phase — identical to what
    %   ca_phase_map would be.  Compare alternans.phase_map against
    %   substrate.phase_map (APD-based) to reveal spatial concordance/discordance
    %   between Ca and AP alternans phase.

    % Per-beat raw maps (for downstream trace extraction)
    alternans.amp_beat        = amp_beat;          % {nBeats x 1}    cell of [R x C] maps
    alternans.diastolic_beat  = diastolic_beat;    % {nBeats x 1}    cell of [R x C] maps
    alternans.CaTD_beat       = CaTD_beat;         % {nBeats x nLev} cell of [R x C] maps
    alternans.TTP_beat        = TTP_beat;          % {nBeats x 1}    cell of [R x C] maps
end


% =========================================================================
function [releaseALT, loadALT, Lmean, Smean, Lidx, Sidx] = ...
        wangLoadRelease(amplitude, diastolicLevel)
% WANGLOADRELEASE  Wang et al. 2014 load and release alternans (1-D).
%
%   Release alternans = 1 - S/L
%     L = mean amplitude of large beats, S = mean amplitude of small beats.
%     Beats are classified as L or S by comparing the mean of odd vs even
%     groups; the larger group becomes L.
%
%   Load alternans = D/L
%     D = mean absolute difference in diastolic level between L and S beats.
%     Normalised by L (large-beat amplitude) per Wang et al. Fig 3B.

    amplitude      = amplitude(:);
    diastolicLevel = diastolicLevel(:);
    nB = length(amplitude);

    odd_idx  = (1:2:nB)';
    even_idx = (2:2:nB)';
    nP       = min(numel(odd_idx), numel(even_idx));
    odd_amp  = amplitude(odd_idx(1:nP));
    even_amp = amplitude(even_idx(1:nP));

    if mean(odd_amp) >= mean(even_amp)
        Lidx = odd_idx(1:nP);
        Sidx = even_idx(1:nP);
    else
        Lidx = even_idx(1:nP);
        Sidx = odd_idx(1:nP);
    end

    Lmean = mean(amplitude(Lidx));
    Smean = mean(amplitude(Sidx));

    % Release alternans (Wang Eq.: 1 − S/L)
    releaseALT = 1 - Smean / (Lmean + eps);

    % Load alternans (Wang Eq.: D/L)
    nD        = min(numel(Lidx), numel(Sidx));
    D         = mean(abs(diastolicLevel(Lidx(1:nD)) - diastolicLevel(Sidx(1:nD))));
    loadALT   = D / (Lmean + eps);
end


% =========================================================================
function [SAI, f_peak] = spectralAlternansIndex1D(amplitude)
    N   = length(amplitude);
    Y   = fft(amplitude - mean(amplitude));
    P   = abs(Y/N).^2;
    f   = (0:N-1)/N;
    altBand = f >= 0.4 & f <= 0.5;
    [SAI, idx] = max(P(altBand));
    fBand   = f(altBand);
    f_peak  = fBand(idx);
end


% =========================================================================
function plotCaAlternans(a)
    nBeats  = a.nBeats;
    oddIdx  = 1:2:nBeats;
    evenIdx = 2:2:nBeats;
    t       = a.time;

    figure('Name','Ca2+ Transient Alternans','Color','w','Position',[40 40 1400 950]);

    subplot(3,3,[1 2 3]);
    plot(t, a.signal,'Color',[0.15 0.15 0.15],'LineWidth',1.2); hold on;
    plot(t(a.peakLocs(oddIdx)),  a.signal(a.peakLocs(oddIdx)), ...
         'bo','MarkerFaceColor','b','MarkerSize',8);
    plot(t(a.peakLocs(evenIdx)), a.signal(a.peakLocs(evenIdx)), ...
         'rs','MarkerFaceColor','r','MarkerSize',8);
    plot(t(a.onsetLocs), a.signal(a.onsetLocs), ...
         'k^','MarkerSize',5,'MarkerFaceColor','k');
    xlabel('Time (ms)'); ylabel('dF/F0 (norm.)');
    title('Ca2+ Transients  (blue=odd, red=even beats)','FontWeight','bold');
    ylim([-0.05 1.15]); grid on;

    subplot(3,3,4);
    bar(oddIdx,  a.amplitude(oddIdx),  0.5,'FaceColor',[0.2 0.4 0.8],'EdgeColor','none'); hold on;
    bar(evenIdx, a.amplitude(evenIdx), 0.5,'FaceColor',[0.8 0.2 0.2],'EdgeColor','none');
    yline(mean(a.amplitude),'k--','LineWidth',1.5);
    xlabel('Beat #'); ylabel('Amplitude (norm.)');
    title(sprintf('Amplitude ALT=%.4f  (%.1f%%)', a.ALT, a.ALT_ratio*100),'FontWeight','bold');
    grid on;

    subplot(3,3,5);
    bar(2:nBeats, a.DI,'FaceColor',[0.3 0.7 0.4],'EdgeColor','none','FaceAlpha',0.85);
    yline(0,'k-','LineWidth',1.5);
    yline( 0.1,'r--','LineWidth',1); yline(-0.1,'r--','LineWidth',1);
    xlabel('Beat #'); ylabel('DI (alternans ratio)');
    title('Beat-to-Beat Alternans Index','FontWeight','bold');
    ylim([-1 1]); grid on;

    subplot(3,3,6);
    plot(1:nBeats, a.D50,'b-o','LineWidth',1.5,'MarkerFaceColor','b','MarkerSize',7); hold on;
    plot(1:nBeats, a.D80,'r-s','LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',7);
    xlabel('Beat #'); ylabel('Duration (ms)');
    title(sprintf('CaTD  D50 ALT=%.1f ms | D80 ALT=%.1f ms', a.ALT_D50, a.ALT_D80),'FontWeight','bold');
    legend('D50','D80'); grid on;

    subplot(3,3,7);
    bar(oddIdx,  a.tau(oddIdx),  0.5,'FaceColor',[0.2 0.4 0.8],'EdgeColor','none'); hold on;
    bar(evenIdx, a.tau(evenIdx), 0.5,'FaceColor',[0.8 0.2 0.2],'EdgeColor','none');
    xlabel('Beat #'); ylabel('tau decay (ms)');
    title('Decay Constant tau','FontWeight','bold'); grid on;

    subplot(3,3,8);
    N  = length(a.amplitude);
    Y  = fft(a.amplitude - mean(a.amplitude));
    P  = abs(Y/N).^2; f = (0:N-1)/N;
    stem(f(1:floor(N/2)), P(1:floor(N/2)),'filled', ...
         'Color',[0.5 0.2 0.7],'LineWidth',1.5,'MarkerSize',6);
    xline(0.5,'r--','0.5 cyc/beat','LineWidth',2,'LabelOrientation','horizontal');
    xlabel('Freq (cyc/beat)'); ylabel('Power');
    title(sprintf('Amplitude Spectrum  SAI=%.4f', a.spectralALT),'FontWeight','bold'); grid on;

    subplot(3,3,9);
    scatter(a.amplitude(1:end-1), a.amplitude(2:end), 60, ...
            1:nBeats-1,'filled','MarkerEdgeColor','k','LineWidth',0.5);
    colormap(gca, parula); colorbar;
    refLine = linspace(min(a.amplitude)*0.9, max(a.amplitude)*1.1, 100);
    hold on; plot(refLine, refLine,'k--','LineWidth',1.5);
    xlabel('A_n'); ylabel('A_{n+1}');
    title('Poincare Plot','FontWeight','bold'); axis equal; grid on;

    sgtitle('Ca2+ Transient Alternans Analysis','FontSize',14,'FontWeight','bold');
end


% =========================================================================
function demoCaAlternans()
% DEMOCAALTERNANS  Simulate Ca2+ transient alternans and run full analysis.
    dt   = 0.5;
    tSim = 0:dt:8000;
    BCL  = 500;
    nB   = floor(tSim(end)/BCL);
    cSim = zeros(size(tSim));
    rng(7);
    for b = 1:nB
        t0    = (b-1)*BCL + 20;
        Amp   = 1.0 - 0.4*mod(b,2) + 0.02*randn;
        tau_r = 15 + randn;
        tau_d = 120 + 20*mod(b,2) + 2*randn;
        tB    = tSim - t0;
        pulse = Amp*(1-exp(-max(tB,0)/tau_r)).*exp(-max(tB,0)/tau_d);
        pulse(tB < 0) = 0;
        cSim  = cSim + pulse;
    end
    cSim = cSim + 0.01*randn(size(cSim));
    analyzeCaTransientAlternans(tSim, cSim, 'PacingRate', 1000/BCL, 'Threshold', 0.25);
end
