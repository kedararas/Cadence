function [alternans] = analyzeAPAlternans(time, voltage, varargin)
% ANALYZEAPALTERNANS  Action potential alternans analysis — 1-D or 3-D.
%
%   ── 1-D mode (single trace) ─────────────────────────────────────────────
%   alternans = analyzeAPAlternans(time, voltage)
%   alternans = analyzeAPAlternans(time, voltage, Name, Value, ...)
%
%   Inputs
%     time    - Time vector (ms), length T
%     voltage - Membrane voltage vector (mV), length T
%
%   ── 3-D mode (full optical-mapping stack) ────────────────────────────────
%   alternans = analyzeAPAlternans(time, voltage3D, 'BeatFrames', bf, ...)
%
%   Inputs
%     time       - Time vector (ms), length T (number of frames)
%     voltage3D  - [rows × cols × T] fluorescence array (any numeric class)
%     BeatFrames - [num_beats × 2] [start_frame, end_frame] per beat (REQUIRED)
%
%   Optional Name-Value Pairs (both modes)
%     'APD_levels'  - Repolarization % levels (default: [30 50 80 90])
%     'PlotResults' - true/false (default: true for 1-D, false for 3-D)
%
%   Optional Name-Value Pairs (1-D mode only)
%     'PacingRate'  - Stimulation rate Hz (default: auto-detect)
%     'Threshold'   - AP detection threshold mV (default: -20 mV)
%     'RestingVm'   - Resting membrane potential mV (default: auto)
%
%   Optional Name-Value Pairs (3-D mode only)
%     'BeatFrames'  - [num_beats × 2] start/end frame index of each beat
%     'Mask'        - [rows × cols] logical/numeric; pixels outside = NaN
%
%   3-D output struct fields
%     .mode          '3D'
%     .APD##_alt     APD alternans map at each requested level (e.g. APD80_alt)
%     .amp_alt_map   Amplitude alternans map
%     .dvdt_alt_map  Max-upstroke-velocity alternans map (normalised units/ms)
%     .tri_alt_map   Triangulation alternans map (APD90 - APD30, ms)
%     .pval_map      Pixel-wise paired t-test p-value (APD80)
%     .tstat_map     Pixel-wise t-statistic
%     .spectral_map  Pixel-wise spectral alternans index (power at 0.5 cyc/beat)
%     .APD_beat      {num_beats × num_levels} cell of per-beat APD maps (ms)
%     .Vamp_beat     {num_beats × 1}  cell of per-beat raw amplitude maps

    %% ---- Parse inputs ----
    p = inputParser;
    addRequired(p,  'time',       @isvector);
    addRequired(p,  'voltage',    @isnumeric);   % vector OR 3-D array
    addParameter(p, 'PacingRate',  [],            @isnumeric);
    addParameter(p, 'Threshold',   -20,           @isnumeric);
    addParameter(p, 'RestingVm',   [],            @isnumeric);
    addParameter(p, 'APD_levels',  [30 50 80], @isvector);
    addParameter(p, 'PlotResults', [],            @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'BeatFrames',  [],            @isnumeric);
    addParameter(p, 'Mask',        [],            @(x) isnumeric(x) || islogical(x));
    parse(p, time, voltage, varargin{:});
    opts = p.Results;

    %% ---- Route: 3-D spatial mode ----------------------------------------
    if ndims(voltage) == 3
        if isempty(opts.BeatFrames)
            error('analyzeAPAlternans:noBeatFrames', ...
                ['BeatFrames is required for 3-D input.\n' ...
                 'Supply a [num_beats x 2] array of [start_frame, end_frame].']);
        end
        if isempty(opts.PlotResults), opts.PlotResults = false; end
        alternans = apAlternans3D(time, voltage, opts);
        return;
    end

    %% ---- 1-D single-trace mode ------------------------------------------
    if isempty(opts.PlotResults), opts.PlotResults = true; end

    time    = time(:);
    voltage = voltage(:);
    dt      = mean(diff(time));     % ms/sample

    fprintf('\n======== Action Potential Alternans Analysis ========\n');
    fprintf('Signal duration    : %.1f ms  (%.0f samples)\n', time(end)-time(1), length(time));
    fprintf('Sampling interval  : %.4f ms\n', dt);

    %% ---- Step 1: Resting membrane potential ----
    if isempty(opts.RestingVm)
        [counts, edges] = histcounts(voltage, 200);
        [~, idx]        = max(counts);
        restVm          = mean(edges(idx:idx+1));
    else
        restVm = opts.RestingVm;
    end
    fprintf('Resting Vm         : %.1f mV\n', restVm);

    %% ---- Step 2: Detect AP upstrokes ----
    aboveThresh = voltage > opts.Threshold;
    rising      = [false; diff(aboveThresh) > 0];

    if ~isempty(opts.PacingRate)
        minBeatDist = round((1000 / opts.PacingRate) * 0.5 / dt);
    else
        minBeatDist = round(0.05 * length(voltage));
    end

    crossIdx = find(rising);
    [peakVoltage, peakLocs, upstrokeIdx] = ...
        refineAPUpstrokes(voltage, crossIdx, minBeatDist, dt);

    nBeats = length(peakLocs);
    fprintf('Detected beats     : %d\n', nBeats);

    if nBeats < 4
        warning('Too few beats (%d) detected. Check threshold or signal.', nBeats);
        alternans = struct(); return;
    end

    %% ---- Step 3: Per-beat AP metrics ----
    APD     = zeros(nBeats, length(opts.APD_levels));
    Vmax    = zeros(nBeats, 1);
    Vrest   = zeros(nBeats, 1);
    Vamp    = zeros(nBeats, 1);
    dVdtMax = zeros(nBeats, 1);
    RMP     = zeros(nBeats, 1);
    TTP     = zeros(nBeats, 1);
    tri     = zeros(nBeats, 1);
    notch   = zeros(nBeats, 1);

    for b = 1:nBeats
        segStart = upstrokeIdx(b);
        if b < nBeats
            segEnd = upstrokeIdx(b+1) - 1;
        else
            segEnd = min(upstrokeIdx(b) + round(600/dt), length(voltage));
        end
        segV = voltage(segStart:segEnd);

        preWin   = max(1, segStart - round(30/dt));
        Vrest(b) = mean(voltage(preWin:segStart));
        RMP(b)   = Vrest(b);
        Vmax(b)  = peakVoltage(b);
        Vamp(b)  = Vmax(b) - Vrest(b);
        TTP(b)   = (peakLocs(b) - upstrokeIdx(b)) * dt;

        upWin      = upstrokeIdx(b) : min(upstrokeIdx(b)+round(5/dt), length(voltage)-1);
        dVdt       = diff(voltage(upWin)) / dt;
        dVdtMax(b) = max(dVdt);

        for lv = 1:length(opts.APD_levels)
            pct     = opts.APD_levels(lv) / 100;
            vThresh = Vmax(b) - pct * Vamp(b);
            repIdx  = find(segV(peakLocs(b)-segStart+1:end) <= vThresh, 1, 'first');
            if ~isempty(repIdx)
                APD(b,lv) = repIdx * dt;
            else
                APD(b,lv) = NaN;
            end
        end

        if ~any(isnan(APD(b,[1 4])))
            tri(b) = APD(b,4) - APD(b,1);
        else
            tri(b) = NaN;
        end

        notchWin = peakLocs(b) : min(peakLocs(b)+round(20/dt), segStart+length(segV)-1);
        notchWin = notchWin(notchWin <= length(voltage));
        notchSeg = voltage(notchWin);
        if length(notchSeg) > 3
            notch(b) = Vmax(b) - min(notchSeg);
        else
            notch(b) = 0;
        end
    end

    apd80col = find(opts.APD_levels == 80, 1);
    apd90col = find(opts.APD_levels == 90, 1);
    if isempty(apd80col), apd80col = size(APD,2); end
    if isempty(apd90col), apd90col = size(APD,2); end
    APD80 = APD(:, apd80col);
    APD90 = APD(:, apd90col);

    %% ---- Step 4: Diastolic interval ----
    DI = zeros(nBeats-1, 1);
    for b = 1:nBeats-1
        repEnd = (upstrokeIdx(b) + round(APD90(b)/dt));
        DI(b)  = max((upstrokeIdx(b+1) - repEnd) * dt, 0);
    end

    %% ---- Step 5: Alternans quantification ----
    [ALT_amp,   DI_amp,   pVal_amp,   evenOdd_amp]   = calcAlternans(Vamp,    'Amplitude');
    [ALT_apd80, DI_apd80, pVal_apd80, evenOdd_apd80] = calcAlternans(APD80,   'APD80');
    [ALT_apd90, DI_apd90, pVal_apd90, evenOdd_apd90] = calcAlternans(APD90,   'APD90');
    [ALT_dvdt,  DI_dvdt,  pVal_dvdt,  evenOdd_dvdt]  = calcAlternans(dVdtMax, 'dVdt');
    [ALT_tri,   ~,        pVal_tri,   evenOdd_tri]    = calcAlternans(tri,     'Triangulation');

    [SAI_amp,  f_amp]  = spectralAlternansIndex(Vamp);
    [SAI_apd,  f_apd]  = spectralAlternansIndex(APD80);
    [SAI_dvdt, ~]      = spectralAlternansIndex(dVdtMax);

    [shapeALT, shapePCs] = morphologicalAlternans(voltage, upstrokeIdx, peakLocs, nBeats, dt);

    %% ---- Step 5b: Wang load and release alternans (Wang et al. 2014) ----
    % Applied to AP signals following Wang et al. 2014 framework:
    %   Release alternans = 1 - S/L  (AP amplitude: large vs small beats)
    %   Load alternans    = D/L      (diastolic Vm alternans / large-beat amplitude)
    % Vrest is the pre-AP diastolic membrane potential — the AP analogue of
    % diastolic SR Ca2+ load. Vamp is the AP-amplitude analogue of Ca2+ release.
    [AP_releaseALT_Wang, AP_loadALT_Wang, AP_wang_L_mean, AP_wang_S_mean, ...
     AP_wang_L_idx, AP_wang_S_idx] = wangLoadReleaseAP(Vamp, Vrest);

    isAlternans = (ALT_apd80/mean(APD80,'omitnan') > 0.05 && pVal_apd80 < 0.05) || ...
                  (ALT_amp/mean(Vamp)               > 0.05 && pVal_amp   < 0.05);

    %% ---- Step 6: Restitution ----
    [restitSlope, restitFit] = calcRestitution(APD90(1:end-1), DI);

    %% ---- Assemble output struct ----
    alternans.mode          = '1D';
    alternans.time          = time;
    alternans.voltage       = voltage;
    alternans.upstrokeIdx   = upstrokeIdx;
    alternans.peakLocs      = peakLocs;
    alternans.nBeats        = nBeats;
    alternans.dt            = dt;
    alternans.restVm        = restVm;
    alternans.Vmax          = Vmax;
    alternans.Vrest         = Vrest;
    alternans.Vamp          = Vamp;
    alternans.dVdtMax       = dVdtMax;
    alternans.APD           = APD;
    alternans.APD80         = APD80;
    alternans.APD90         = APD90;
    alternans.APD_levels    = opts.APD_levels;
    alternans.TTP           = TTP;
    alternans.tri           = tri;
    alternans.notch         = notch;
    alternans.DI            = DI;
    alternans.ALT_amp       = ALT_amp;
    alternans.ALT_apd80     = ALT_apd80;
    alternans.ALT_apd90     = ALT_apd90;
    alternans.ALT_dvdt      = ALT_dvdt;
    alternans.ALT_tri       = ALT_tri;
    alternans.DI_amp        = DI_amp;
    alternans.DI_apd80      = DI_apd80;
    alternans.DI_dvdt       = DI_dvdt;
    alternans.pVal_amp      = pVal_amp;
    alternans.pVal_apd80    = pVal_apd80;
    alternans.pVal_apd90    = pVal_apd90;
    alternans.pVal_dvdt     = pVal_dvdt;
    alternans.pVal_tri      = pVal_tri;
    alternans.SAI_amp       = SAI_amp;
    alternans.SAI_apd       = SAI_apd;
    alternans.SAI_dvdt      = SAI_dvdt;
    alternans.f_alt_amp     = f_amp;
    alternans.f_alt_apd     = f_apd;
    alternans.shapeALT          = shapeALT;
    alternans.shapePCs          = shapePCs;
    alternans.evenOdd_amp       = evenOdd_amp;

    % Alternans ratio: normalised ALT magnitude per metric (0..1 range)
    % ratio = ALT / mean_metric  (same as evenOdd.ratio but stored top-level)
    alternans.ALT_ratio_amp   = evenOdd_amp.ratio;
    alternans.ALT_ratio_apd80 = evenOdd_apd80.ratio;
    alternans.ALT_ratio_apd90 = evenOdd_apd90.ratio;
    alternans.ALT_ratio_dvdt  = evenOdd_dvdt.ratio;
    alternans.ALT_ratio_tri   = evenOdd_tri.ratio;

    % Phase: sign of the mean odd−even APD80 difference (+1 or −1).
    % Tells you which beat class (odd/even) carries the longer AP in this
    % recording — useful for aligning phase across multi-site recordings.
    nP80 = evenOdd_apd80.nPairs;
    alternans.alt_phase = sign(mean( ...
        evenOdd_apd80.odd(1:nP80) - evenOdd_apd80.even(1:nP80)));

    alternans.restitSlope       = restitSlope;
    alternans.restitFit         = restitFit;
    alternans.isAlternans       = isAlternans;
    % Wang et al. 2014 load and release alternans for AP
    alternans.AP_releaseALT_Wang = AP_releaseALT_Wang;  % 1 - S/L (AP amplitude)
    alternans.AP_loadALT_Wang    = AP_loadALT_Wang;     % D/L (diastolic Vm)
    alternans.AP_wang_L_mean     = AP_wang_L_mean;      % mean large-beat Vamp
    alternans.AP_wang_S_mean     = AP_wang_S_mean;      % mean small-beat Vamp
    alternans.AP_wang_L_idx      = AP_wang_L_idx;       % large-beat indices
    alternans.AP_wang_S_idx      = AP_wang_S_idx;       % small-beat indices

    printSummary(alternans);
    if opts.PlotResults
        plotAPAlternans(alternans);
    end
end


% =========================================================================
function alternans = apAlternans3D(time, voltage3D, opts)
% APALTERNANS3D  Vectorized AP alternans analysis across a full [R x C x T]
%               optical-mapping stack.
%
%   For each beat window the signal is per-pixel normalised to [0,1] so
%   the APD threshold formula (peak - pct*amplitude) collapses to the
%   simple scalar  thresh = 1 - APD_level/100.  This makes it possible
%   to process all pixels simultaneously using cumsum-based crossing
%   detection instead of a pixel loop.
%
%   Alternans magnitude per pixel: mean(APD_odd - APD_even) over all pairs.
%   Statistical map: vectorised paired t-test on APD80 (odd vs even).
%   Spectral map: FFT of the per-beat APD80 series along the beat dimension.

    dt         = mean(diff(time(:)));    % ms per frame
    [R, C, ~]  = size(voltage3D);
    APD_levels = opts.APD_levels(:)';    % row vector, e.g. [30 50 80 90]
    nLev       = numel(APD_levels);
    bf         = opts.BeatFrames;        % [num_beats x 2]
    nBeats     = size(bf, 1);

    % ── Mask ───────────────────────────────────────────────────────────────
    if ~isempty(opts.Mask)
        % nan_mask = double(logical(opts.Mask));
        % nan_mask(nan_mask == 0) = nan;
        nan_mask = opts.Mask;
    else
        nan_mask = ones(R, C);
    end

    % ── Per-beat storage ───────────────────────────────────────────────────
    APD_beat      = cell(nBeats, nLev);  % APD maps  [R x C] per beat per level
    Vamp_beat     = cell(nBeats, 1);     % raw amplitude (pre-normalisation)
    dVdt_beat     = cell(nBeats, 1);     % max upstroke velocity (norm. units/ms)
    act_beat      = cell(nBeats, 1);     % activation time from window start (ms)
    diastolic_beat = cell(nBeats, 1);    % pre-AP diastolic Vm level per beat [R x C]

    for j = 1:nBeats
        sf    = bf(j, 1);
        ef    = bf(j, 2);
        T_win = ef - sf + 1;

        raw = double(voltage3D(:,:, sf:ef));   % R x C x T_win

        % Per-pixel [0,1] normalisation
        lo  = min(raw, [], 3);
        hi  = max(raw, [], 3);
        rng = max(hi - lo, eps('single'));
        wn  = (raw - lo) ./ rng;              % R x C x T_win

        % Store raw amplitude for amplitude-alternans map
        Vamp_beat{j} = (hi - lo) .* nan_mask;

        % Diastolic Vm: minimum of the beat window (pre-AP resting level)
        % This is the AP analogue of diastolic SR Ca2+ load (Wang et al. 2014)
        diastolic_beat{j} = lo .* nan_mask;

        % Activation: frame of max temporal derivative (max upstroke)
        dv           = diff(wn, 1, 3);                           % R x C x (T_win-1)
        [dvMax, af]  = max(dv, [], 3);                           % R x C
        act_beat{j}  = double(af) * dt .* nan_mask;             % ms, NaN outside mask
        dVdt_beat{j} = dvMax / dt .* nan_mask;                  % norm.units/ms

        % Peak location (for "after peak" gate)
        [~, pk] = max(wn, [], 3);                                % R x C
        t_idx   = reshape(1:T_win, 1, 1, T_win);
        after_pk = t_idx > reshape(pk, R, C, 1);                % R x C x T_win

        % APD at each repolarisation level (vectorised cumsum crossing)
        for lv = 1:nLev
            thresh        = 1 - APD_levels(lv) / 100;           % e.g. 0.20 for APD80
            below         = wn <= thresh;
            [hit, rf]     = max(cumsum(below & after_pk, 3) == 1, [], 3);
            apd           = (double(rf) - double(af)) * dt;     % ms
            invalid       = ~logical(hit) | apd <= 0 | isnan(act_beat{j});
            apd(invalid)  = nan;
            APD_beat{j, lv} = apd;
        end
    end

    % ── Per-level reachability diagnostic ──────────────────────────────────
    % Fraction of (in-mask pixel × beat) entries that yield a finite APD at
    % each repolarisation level.  A low value at the deepest level means the AP
    % does not repolarise that far within the beat window (common at fast
    % pacing) -> the corresponding alternans map will be largely NaN.
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
            m   = APD_beat{j, lv};
            cnt = cnt + nnz(isfinite(m) & tissue);
        end
        level_coverage(lv) = cnt / (nT * nBeats);
    end
    fprintf('AP alternans reachability (tissue-beats with measurable APD):\n');
    for lv = 1:nLev
        fprintf('  APD%-2d : %5.1f%%\n', APD_levels(lv), 100*level_coverage(lv));
    end
    if level_coverage(end) < 0.5
        warning('analyzeAPAlternans:lowCoverage', ...
            ['APD%d reached in only %.0f%% of tissue-beats; its alternans map ' ...
             'will be largely NaN. The AP likely does not repolarise to %d%% ' ...
             'within the beat window at this pacing rate -- consider a lower ' ...
             'APD level.'], APD_levels(end), 100*level_coverage(end), APD_levels(end));
    end

    % ── Odd / even beat split ──────────────────────────────────────────────
    odd_idx  = 1:2:nBeats;
    even_idx = 2:2:nBeats;
    n_pairs  = min(numel(odd_idx), numel(even_idx));

    if n_pairs < 1
        warning('analyzeAPAlternans:toofewbeats', ...
            'Need at least 2 beats for alternans. Returning empty struct.');
        alternans = struct('mode','3D'); return;
    end

    % ── APD alternans maps at each level ───────────────────────────────────
    APD_alt_maps = cell(nLev, 1);
    for lv = 1:nLev
        odd_s  = cat(3, APD_beat{odd_idx(1:n_pairs),  lv});   % R x C x n_pairs
        even_s = cat(3, APD_beat{even_idx(1:n_pairs), lv});
        APD_alt_maps{lv} = mean(odd_s - even_s, 3, 'omitnan');
    end

    % ── Amplitude alternans ────────────────────────────────────────────────
    odd_amp  = cat(3, Vamp_beat{odd_idx(1:n_pairs)});
    even_amp = cat(3, Vamp_beat{even_idx(1:n_pairs)});
    amp_alt  = mean(odd_amp - even_amp, 3, 'omitnan');

    % ── dVdt alternans ─────────────────────────────────────────────────────
    odd_dv   = cat(3, dVdt_beat{odd_idx(1:n_pairs)});
    even_dv  = cat(3, dVdt_beat{even_idx(1:n_pairs)});
    dvdt_alt = mean(odd_dv - even_dv, 3, 'omitnan');

    % ── Wang et al. 2014: AP release and load alternans maps (per-pixel L/S) ─
    % Per-pixel classification: at each pixel independently, the group with
    % the larger mean Vamp is L and the smaller is S.  A global classification
    % would invert L/S across minority-phase regions in discordant alternans,
    % producing negative AP_release_alt values.  Per-pixel classification
    % guarantees AP_release_alt_map in [0,1] and AP_load_alt_map >= 0.
    %
    % odd_larger_ap: [R x C] logical — true where odd beats have larger Vamp.
    % This is the Vamp-based phase, which should align with phase_map (APD-based)
    % under normal Ca-driven alternans but can diverge in voltage-driven cases.
    odd_vamp_mean  = mean(odd_amp,  3, 'omitnan');   % R x C
    even_vamp_mean = mean(even_amp, 3, 'omitnan');   % R x C
    odd_larger_ap  = odd_vamp_mean >= even_vamp_mean;  % R x C logical

    odd_diast_stack  = cat(3, diastolic_beat{odd_idx(1:n_pairs)});
    even_diast_stack = cat(3, diastolic_beat{even_idx(1:n_pairs)});

    % Initialise L stacks as odd, then swap pixels where even is locally larger
    L_vamp_stack  = odd_amp;
    S_vamp_stack  = even_amp;
    L_diast_stack = odd_diast_stack;
    S_diast_stack = even_diast_stack;

    swap_ap = repmat(~odd_larger_ap, 1, 1, n_pairs);   % R x C x n_pairs swap mask
    L_vamp_stack(swap_ap)  = even_amp(swap_ap);
    S_vamp_stack(swap_ap)  = odd_amp(swap_ap);
    L_diast_stack(swap_ap) = even_diast_stack(swap_ap);
    S_diast_stack(swap_ap) = odd_diast_stack(swap_ap);

    L_vamp_map = mean(L_vamp_stack, 3, 'omitnan');   % R x C, always >= S_vamp_map
    S_vamp_map = mean(S_vamp_stack, 3, 'omitnan');   % R x C

    % AP release alternans map: 1 - S/L  (Wang et al. 2014) — result in [0,1]
    AP_release_alt_map = 1 - S_vamp_map ./ (L_vamp_map + eps);
    % AP load alternans map: D/L — diastolic Vm alternans / large-beat amplitude
    D_vmap          = mean(abs(L_diast_stack - S_diast_stack), 3, 'omitnan');
    AP_load_alt_map = D_vmap ./ (L_vamp_map + eps);

    % ── Triangulation alternans (APD80 - APD30 if both levels present) ─────
    apd30_col = find(APD_levels == 30, 1);
    apd80_col = find(APD_levels == 80, 1);
    tri_alt   = [];
    if ~isempty(apd30_col) && ~isempty(apd80_col)
        tri_beat = cellfun(@(a,b) a - b, ...
            APD_beat(:, apd80_col), APD_beat(:, apd30_col), ...
            'UniformOutput', false);
        odd_tri  = cat(3, tri_beat{odd_idx(1:n_pairs)});
        even_tri = cat(3, tri_beat{even_idx(1:n_pairs)});
        tri_alt  = mean(odd_tri - even_tri, 3, 'omitnan');
    end

    % ── Paired t-test map on APD80 (requires >= 3 pairs) ──────────────────
    apd80_col = find(APD_levels == 80, 1);
    if isempty(apd80_col), apd80_col = min(3, nLev); end

    pval_map  = nan(R, C);
    tstat_map = nan(R, C);
    if n_pairs >= 3
        odd_80  = cat(3, APD_beat{odd_idx(1:n_pairs),  apd80_col});  % R x C x n_pairs
        even_80 = cat(3, APD_beat{even_idx(1:n_pairs), apd80_col});
        d       = odd_80 - even_80;                                   % R x C x n_pairs
        d_mean  = mean(d, 3, 'omitnan');
        d_std   = std(d,  0, 3, 'omitnan');
        tstat   = d_mean ./ (d_std / sqrt(n_pairs));
        tstat_map = tstat;
        try
            % Two-tailed p-value from t-distribution (Statistics Toolbox)
            pval_map = 2 * (1 - tcdf(abs(tstat), n_pairs - 1));
        catch
            % Approximate via normal distribution if tcdf unavailable
            pval_map = 2 * (1 - normcdf(abs(tstat)));
        end
    end

    % ── Spectral alternans index map (FFT along beat dimension) ───────────
    % Power at 0.5 cycles/beat flags a true beat-alternating pattern.
    spectral_map = nan(R, C);
    if nBeats >= 4
        all_apd = cat(3, APD_beat{:, apd80_col});           % R x C x nBeats
        all_apd(isnan(all_apd)) = 0;                        % replace NaN → 0
        all_apd = all_apd - mean(all_apd, 3);               % detrend (zero-mean)
        Y        = fft(all_apd, [], 3);
        P        = abs(Y / nBeats).^2;
        % Index of 0.5 cycles/beat = beat index N/2+1 (for even N)
        alt_bin  = floor(nBeats/2) + 1;
        if alt_bin <= size(P, 3)
            spectral_map = P(:,:, alt_bin);
        end
    end

    % ── Alternans ratio maps ───────────────────────────────────────────────
    % ratio = |ALT| / mean_metric per pixel  (dimensionless).
    % Normalising removes absolute-value differences across the field and
    % highlights pixels where the alternating fraction is large regardless
    % of local signal amplitude.
    APD_ratio_maps = cell(nLev, 1);
    for lv = 1:nLev
        apd_mean_lv        = mean(cat(3, APD_beat{:, lv}), 3, 'omitnan');
        APD_ratio_maps{lv} = abs(APD_alt_maps{lv}) ./ (apd_mean_lv + eps);
    end

    amp_mean_map  = mean(cat(3, Vamp_beat{:}),  3, 'omitnan');
    amp_ratio_map = abs(amp_alt) ./ (amp_mean_map + eps);

    dvdt_mean_map  = mean(cat(3, dVdt_beat{:}), 3, 'omitnan');
    dvdt_ratio_map = abs(dvdt_alt) ./ (dvdt_mean_map + eps);

    % ── Phase maps ─────────────────────────────────────────────────────────
    % Binary phase map  (+1 / -1 / 0)
    %   +1 : odd beats have longer APD80 than even beats at this pixel
    %   -1 : odd beats have shorter APD80
    %    0 : no APD data (NaN alternans)
    % Regions of opposite sign are DISCORDANT — the spatial boundary between
    % +1 and -1 zones is a nodal line (highest reentry risk).
    apd80_alt_map = APD_alt_maps{apd80_col};
    phase_map     = sign(apd80_alt_map);

    % Continuous phase angle  (radians, −π to +π)
    % Computed as the complex angle of the Fourier coefficient at exactly
    % 0.5 cycles/beat.  Unlike the binary map, this reveals gradual phase
    % gradients across the field — important for identifying slow phase
    % transitions that precede discordant alternans.
    %
    % Interpretation:
    %   pixels with angle ≈ 0 and ≈ ±π are in anti-phase (discordant)
    %   pixels with similar angles are concordant
    %   a smooth spatial angle gradient → travelling alternans wave
    phase_angle_map = nan(R, C);
    if nBeats >= 4
        all_apd80 = cat(3, APD_beat{:, apd80_col});  % R x C x nBeats
        all_apd80(isnan(all_apd80)) = 0;
        all_apd80 = all_apd80 - mean(all_apd80, 3);  % zero-mean detrend
        Y80      = fft(all_apd80, [], 3);
        alt_bin  = floor(nBeats/2) + 1;
        if alt_bin <= size(Y80, 3)
            phase_angle_map = angle(Y80(:,:, alt_bin));  % [-pi, pi]
        end
    end

    % ── Package output ─────────────────────────────────────────────────────
    alternans.mode         = '3D';
    alternans.num_beats    = nBeats;
    alternans.n_pairs      = n_pairs;
    alternans.APD_levels   = APD_levels;
    alternans.dt           = dt;
    alternans.level_coverage = level_coverage;   % fraction of tissue-beats reachable per APD level

    % Named alternans maps for each APD level (e.g. alternans.APD80_alt)
    for lv = 1:nLev
        alternans.(sprintf('APD%d_alt', APD_levels(lv))) = APD_alt_maps{lv};
    end
    alternans.APD_alt_maps  = APD_alt_maps;  % cell array for programmatic access
    alternans.amp_alt_map   = amp_alt;
    alternans.dvdt_alt_map  = dvdt_alt;
    alternans.tri_alt_map   = tri_alt;

    % Statistical and spectral maps
    alternans.pval_map      = pval_map;
    alternans.tstat_map     = tstat_map;
    alternans.spectral_map  = spectral_map;

    % Alternans ratio maps  (|ALT| / mean_metric per pixel)
    alternans.amp_ratio_map  = amp_ratio_map;
    alternans.dvdt_ratio_map = dvdt_ratio_map;
    for lv = 1:nLev
        alternans.(sprintf('APD%d_ratio_map', APD_levels(lv))) = APD_ratio_maps{lv};
    end
    alternans.APD_ratio_maps = APD_ratio_maps;   % cell for programmatic access

    % Phase maps
    alternans.phase_map       = phase_map;        % binary  +1 / -1 / 0
    alternans.phase_angle_map = phase_angle_map;  % continuous  [-pi, pi] rad

    % Wang et al. 2014 load and release alternans maps for AP (per-pixel L/S)
    alternans.AP_release_alt_map = AP_release_alt_map;  % 1 - S/L per pixel, in [0,1]
    alternans.AP_load_alt_map    = AP_load_alt_map;     % D/L   per pixel, >= 0
    alternans.AP_L_vamp_map      = L_vamp_map;          % mean large-beat Vamp map
    alternans.AP_S_vamp_map      = S_vamp_map;          % mean small-beat Vamp map
    alternans.AP_amp_phase_map   = 2*double(odd_larger_ap) - 1;  % +1 odd=L, -1 even=L
    %   AP_amp_phase_map: Vamp-based phase, comparable to phase_map (APD-based)
    %   from assess_arrhythmia_substrate.  Agreement confirms Ca-driven mechanism;
    %   divergence (Vamp phase ≠ APD phase) points to voltage-driven alternans.

    % Per-beat raw maps (for downstream use, e.g. extracting single-pixel traces)
    alternans.APD_beat       = APD_beat;        % {nBeats x nLev} cell of [R x C] maps
    alternans.Vamp_beat      = Vamp_beat;       % {nBeats x 1}   cell of [R x C] maps
    alternans.dVdt_beat      = dVdt_beat;       % {nBeats x 1}   cell of [R x C] maps
    alternans.act_beat       = act_beat;        % {nBeats x 1}   cell of [R x C] maps (ms)
    alternans.diastolic_beat = diastolic_beat;  % {nBeats x 1}   cell of [R x C] maps
end


% =========================================================================
function [peakV, peakLocs, upIdx] = refineAPUpstrokes(voltage, crossIdx, minDist, dt)
    keep = true(size(crossIdx));
    for i = 2:length(crossIdx)
        if crossIdx(i) - crossIdx(i-1) < minDist
            keep(i) = false;
        end
    end
    crossIdx = crossIdx(keep);

    peakV    = zeros(size(crossIdx));
    peakLocs = zeros(size(crossIdx));
    upIdx    = zeros(size(crossIdx));
    searchWin = round(50 / dt);

    for i = 1:length(crossIdx)
        winEnd    = min(crossIdx(i) + searchWin, length(voltage));
        [peakV(i), relIdx] = max(voltage(crossIdx(i):winEnd));
        peakLocs(i) = crossIdx(i) + relIdx - 1;

        dvdtWin = max(1,crossIdx(i)-round(2/dt)) : min(crossIdx(i)+round(5/dt), length(voltage)-1);
        dv      = diff(voltage(dvdtWin)) / dt;
        [~, ri] = max(dv);
        upIdx(i) = dvdtWin(ri);
    end
end


% =========================================================================
function [ALT, DI_ratio, pVal, evenOdd] = calcAlternans(metric, name)
    metric = metric(:);
    valid  = ~isnan(metric);
    m      = metric(valid);
    n      = length(m);

    odd    = m(1:2:end);
    even   = m(2:2:end);
    nPairs = min(length(odd), length(even));
    ALT    = mean(abs(even(1:nPairs) - odd(1:nPairs)));
    meanM  = mean(m);

    DI_ratio = zeros(n-1, 1);
    for i = 1:n-1
        DI_ratio(i) = (m(i+1)-m(i)) / (abs(m(i+1))+abs(m(i)) + eps);
    end

    if nPairs >= 3
        [~, pVal] = ttest(odd(1:nPairs), even(1:nPairs));
    else
        pVal = NaN;
    end

    evenOdd.odd    = odd;
    evenOdd.even   = even;
    evenOdd.nPairs = nPairs;
    evenOdd.ratio  = ALT / (meanM + eps);

    fprintf('  %-18s ALT=%-8.3f  ratio=%.2f%%  p=%.4f\n', ...
            name, ALT, evenOdd.ratio*100, pVal);
end


% =========================================================================
function [SAI, f_peak] = spectralAlternansIndex(metric)
    m = metric(~isnan(metric));
    m = m - mean(m);
    N = length(m);
    if N < 4, SAI=0; f_peak=0.5; return; end

    Y    = fft(m);
    P    = abs(Y/N).^2;
    f    = (0:N-1)/N;
    band = f >= 0.40 & f <= 0.50;
    [SAI, i] = max(P(band));
    fb   = f(band);
    f_peak = fb(i);
end


% =========================================================================
function [shapeALT, PCs] = morphologicalAlternans(voltage, upIdx, pkLocs, nBeats, dt)
    winSamp = round(400 / dt);
    allAPs  = NaN(nBeats, winSamp);
    for b = 1:nBeats
        e   = min(upIdx(b) + winSamp - 1, length(voltage));
        seg = voltage(upIdx(b):e);
        allAPs(b, 1:length(seg)) = seg;
    end
    validRows = ~any(isnan(allAPs), 2);
    APmatrix  = allAPs(validRows, :);
    if size(APmatrix,1) < 4, shapeALT=0; PCs=[]; return; end

    [coeff, score, ~, ~, explained] = pca(APmatrix);
    PCs.coeff    = coeff(:, 1:min(3,size(coeff,2)));
    PCs.score    = score(:, 1:min(3,size(score,2)));
    PCs.explained = explained(1:min(3,end));

    pc1      = score(:,1);
    odd      = pc1(1:2:end);
    even     = pc1(2:2:end);
    nP       = min(length(odd), length(even));
    shapeALT = mean(abs(even(1:nP) - odd(1:nP)));
end


% =========================================================================
function [slope, fitResult] = calcRestitution(APD90, DI)
    valid = ~isnan(APD90) & ~isnan(DI) & DI > 0;
    x     = DI(valid);
    y     = APD90(valid);
    if length(x) < 3, slope = NaN; fitResult = []; return; end

    try
        fitResult = fit(x, y, 'a*(1-exp(-x/b))+c', ...
                        'StartPoint', [50 100 200], ...
                        'Lower', [0 1 0], 'Upper', [500 2000 500]);
        minDI = min(x);
        slope = fitResult.a / fitResult.b * exp(-minDI / fitResult.b);
    catch
        p         = polyfit(x, y, 1);
        slope     = p(1);
        fitResult = p;
    end
end


% =========================================================================
function printSummary(a)
    fprintf('\n--- Per-Beat Metrics (mean ± std) ---\n');
    fprintf('Vamp (mV)          : %.2f ± %.2f\n',  mean(a.Vamp),    std(a.Vamp));
    fprintf('Vmax (mV)          : %.2f ± %.2f\n',  mean(a.Vmax),    std(a.Vmax));
    fprintf('dVdt_max (V/s)     : %.1f ± %.1f\n',  mean(a.dVdtMax), std(a.dVdtMax));
    fprintf('APD80 (ms)         : %.2f ± %.2f\n',  mean(a.APD80,'omitnan'), std(a.APD80,'omitnan'));
    fprintf('APD90 (ms)         : %.2f ± %.2f\n',  mean(a.APD90,'omitnan'), std(a.APD90,'omitnan'));
    fprintf('Triangulation (ms) : %.2f ± %.2f\n',  mean(a.tri,'omitnan'),   std(a.tri,'omitnan'));
    fprintf('DI (ms)            : %.2f ± %.2f\n',  mean(a.DI),      std(a.DI));

    fprintf('\n--- Alternans Summary ---\n');
    fprintf('  %-18s %-10s %-10s %-8s\n','Metric','ALT','ALT%%','p-value');
    fprintf('  %-18s %-10.3f %-10.2f %-8.4f\n','Amplitude (mV)', ...
            a.ALT_amp,   a.ALT_amp/mean(a.Vamp)*100,                        a.pVal_amp);
    fprintf('  %-18s %-10.3f %-10.2f %-8.4f\n','APD80 (ms)', ...
            a.ALT_apd80, a.ALT_apd80/mean(a.APD80,'omitnan')*100,           a.pVal_apd80);
    fprintf('  %-18s %-10.3f %-10.2f %-8.4f\n','APD90 (ms)', ...
            a.ALT_apd90, a.ALT_apd90/mean(a.APD90,'omitnan')*100,           a.pVal_apd90);
    fprintf('  %-18s %-10.3f %-10.2f %-8.4f\n','dVdt (V/s)', ...
            a.ALT_dvdt,  a.ALT_dvdt/mean(a.dVdtMax)*100,                    a.pVal_dvdt);
    fprintf('  %-18s %-10.3f %-10.2f %-8.4f\n','Triangulation', ...
            a.ALT_tri,   a.ALT_tri/mean(a.tri,'omitnan')*100,               a.pVal_tri);

    fprintf('\n--- Spectral Alternans Index ---\n');
    fprintf('  SAI amplitude    : %.4f  (f=%.3f cyc/beat)\n', a.SAI_amp,  a.f_alt_amp);
    fprintf('  SAI APD80        : %.4f  (f=%.3f cyc/beat)\n', a.SAI_apd,  a.f_alt_apd);
    fprintf('  SAI dVdt         : %.4f\n', a.SAI_dvdt);
    fprintf('  Shape ALT (PC1)  : %.4f\n', a.shapeALT);

    fprintf('\n--- Wang et al. 2014 AP Load and Release Alternans ---\n');
    fprintf('  Release ALT (1-S/L): %.4f  (L=%.3f mV, S=%.3f mV)\n', ...
            a.AP_releaseALT_Wang, a.AP_wang_L_mean, a.AP_wang_S_mean);
    fprintf('  Load ALT    (D/L)  : %.4f  (diastolic Vm alternans norm. by L)\n', ...
            a.AP_loadALT_Wang);

    fprintf('\n--- Restitution ---\n');
    fprintf('  Max slope        : %.4f\n', a.restitSlope);
    if a.restitSlope > 1
        fprintf('  >> Slope > 1: UNSTABLE (alternans expected)\n');
    else
        fprintf('  >> Slope < 1: stable\n');
    end

    if a.isAlternans
        fprintf('\n>> AP ALTERNANS DETECTED\n');
    else
        fprintf('\n>> No significant alternans detected.\n');
    end
    fprintf('======================================================\n\n');
end


% =========================================================================
function plotAPAlternans(a)
    nB      = a.nBeats;
    oddIdx  = 1:2:nB;
    evenIdx = 2:2:nB;
    beats   = 1:nB;
    colOdd  = [0.15 0.35 0.75];
    colEven = [0.80 0.15 0.15];

    figure('Name','AP Alternans Analysis','Color','w','Position',[30 30 1500 1000]);

    ax1 = subplot(4,3,[1 2 3]);
    plot(a.time, a.voltage,'Color',[0.1 0.1 0.1],'LineWidth',1.0); hold on;
    plot(a.time(a.peakLocs(oddIdx)),  a.voltage(a.peakLocs(oddIdx)), ...
         'o','Color',colOdd, 'MarkerFaceColor',colOdd, 'MarkerSize',7);
    plot(a.time(a.peakLocs(evenIdx)), a.voltage(a.peakLocs(evenIdx)), ...
         's','Color',colEven,'MarkerFaceColor',colEven,'MarkerSize',7);
    xlabel('Time (ms)'); ylabel('Vm (mV)');
    title('Action Potentials  (blue=odd, red=even beats)','FontWeight','bold');
    grid on; ylim([-100 60]);

    subplot(4,3,4);
    hold on;
    winSamp = round(400 / a.dt);
    for b = 1:nB
        e   = min(a.upstrokeIdx(b)+winSamp-1, length(a.voltage));
        seg = a.voltage(a.upstrokeIdx(b):e);
        t_  = (0:length(seg)-1)*a.dt;
        clr = colOdd*(mod(b,2)==1) + colEven*(mod(b,2)==0);
        plot(t_, seg,'Color',[clr 0.5],'LineWidth',0.9);
    end
    xlabel('Time from upstroke (ms)'); ylabel('Vm (mV)');
    title('Overlaid AP Morphology','FontWeight','bold'); grid on; ylim([-100 60]);

    subplot(4,3,5);
    bar(oddIdx,  a.APD80(oddIdx),  0.4,'FaceColor',colOdd, 'EdgeColor','none'); hold on;
    bar(evenIdx, a.APD80(evenIdx), 0.4,'FaceColor',colEven,'EdgeColor','none');
    yline(mean(a.APD80,'omitnan'),'b--','LineWidth',1.5);
    xlabel('Beat #'); ylabel('APD (ms)');
    title(sprintf('APD80 ALT=%.2f ms  (%.1f%%)  p=%.3f', ...
          a.ALT_apd80, a.ALT_apd80/mean(a.APD80,'omitnan')*100, a.pVal_apd80),'FontWeight','bold');
    grid on;

    subplot(4,3,6);
    bar(oddIdx,  a.Vamp(oddIdx),  0.4,'FaceColor',colOdd, 'EdgeColor','none'); hold on;
    bar(evenIdx, a.Vamp(evenIdx), 0.4,'FaceColor',colEven,'EdgeColor','none');
    xlabel('Beat #'); ylabel('Amplitude (mV)');
    title(sprintf('Amplitude ALT=%.2f mV  p=%.3f', a.ALT_amp, a.pVal_amp),'FontWeight','bold');
    grid on;

    subplot(4,3,7);
    bar(oddIdx,  a.dVdtMax(oddIdx),  0.4,'FaceColor',colOdd, 'EdgeColor','none'); hold on;
    bar(evenIdx, a.dVdtMax(evenIdx), 0.4,'FaceColor',colEven,'EdgeColor','none');
    xlabel('Beat #'); ylabel('dV/dt_{max} (V/s)');
    title(sprintf('dVdt ALT=%.2f  p=%.3f', a.ALT_dvdt, a.pVal_dvdt),'FontWeight','bold'); grid on;

    subplot(4,3,8);
    bar(beats, a.tri, 0.6,'FaceColor',[0.4 0.7 0.3],'EdgeColor','none');
    yline(mean(a.tri,'omitnan'),'k--','LineWidth',1.5);
    xlabel('Beat #'); ylabel('APD90-APD30 (ms)');
    title(sprintf('Triangulation ALT=%.2f ms', a.ALT_tri),'FontWeight','bold'); grid on;

    subplot(4,3,9);
    bar(2:nB, a.DI_apd80,'FaceColor',[0.5 0.3 0.7],'EdgeColor','none','FaceAlpha',0.85);
    yline(0,'k-','LineWidth',1.5);
    yline( 0.05,'r--','LineWidth',1); yline(-0.05,'r--','LineWidth',1);
    xlabel('Beat #'); ylabel('Alternans Index');
    title('Beat-to-Beat APD80 Alternans Index','FontWeight','bold');
    ylim([-1 1]); grid on;

    subplot(4,3,10);
    scatter(a.APD80(1:end-1), a.APD80(2:end), 60, 1:nB-1,'filled','MarkerEdgeColor','k');
    colormap(gca, parula); colorbar;
    rl = linspace(min(a.APD80)*0.9, max(a.APD80)*1.1, 100);
    hold on; plot(rl,rl,'k--','LineWidth',1.5);
    xlabel('APD80_n (ms)'); ylabel('APD80_{n+1} (ms)');
    title('Poincaré Plot (APD80)','FontWeight','bold'); axis equal; grid on;

    subplot(4,3,11);
    N  = sum(~isnan(a.APD80));
    m  = a.APD80(~isnan(a.APD80)); m = m - mean(m);
    Y  = fft(m); P = abs(Y/N).^2; f = (0:N-1)/N;
    stem(f(1:floor(N/2)), P(1:floor(N/2)),'filled','Color',[0.5 0.2 0.7],...
         'LineWidth',1.5,'MarkerSize',6);
    xline(0.5,'r--','0.5 cyc/beat','LineWidth',2,'LabelOrientation','horizontal');
    xlabel('Freq (cyc/beat)'); ylabel('Power');
    title(sprintf('APD80 Spectrum  SAI=%.4f', a.SAI_apd),'FontWeight','bold'); grid on;

    subplot(4,3,12);
    scatter(a.DI, a.APD90(1:end-1), 50, beats(1:end-1),'filled','MarkerEdgeColor','k');
    colormap(gca, cool); colorbar;
    xlabel('Diastolic Interval (ms)'); ylabel('APD90 (ms)');
    title(sprintf('APD Restitution  slope=%.3f', a.restitSlope),'FontWeight','bold'); grid on;

    sgtitle('Action Potential Alternans Analysis','FontSize',14,'FontWeight','bold');
end


% =========================================================================
function [releaseALT, loadALT, Lmean, Smean, Lidx, Sidx] = ...
        wangLoadReleaseAP(Vamp, Vrest)
% WANGLOADRELEASEAP  Wang et al. 2014 load/release framework applied to APs.
%
%   AP amplitude (Vamp) maps to Ca2+ release: large vs small beat amplitude.
%   Diastolic Vm (Vrest) maps to SR Ca2+ load: pre-beat membrane potential.
%
%   Release alternans = 1 - S/L  (AP amplitude ratio, Wang Eq.)
%   Load alternans    = D/L      (|Vrest_L - Vrest_S| / L, Wang Eq.)

    Vamp  = Vamp(:);
    Vrest = Vrest(:);
    nB    = length(Vamp);

    odd_idx  = (1:2:nB)';
    even_idx = (2:2:nB)';
    nP       = min(numel(odd_idx), numel(even_idx));
    odd_amp  = Vamp(odd_idx(1:nP));
    even_amp = Vamp(even_idx(1:nP));

    if mean(odd_amp) >= mean(even_amp)
        Lidx = odd_idx(1:nP);
        Sidx = even_idx(1:nP);
    else
        Lidx = even_idx(1:nP);
        Sidx = odd_idx(1:nP);
    end

    Lmean = mean(Vamp(Lidx));
    Smean = mean(Vamp(Sidx));

    releaseALT = 1 - Smean / (Lmean + eps);

    nD       = min(numel(Lidx), numel(Sidx));
    D        = mean(abs(Vrest(Lidx(1:nD)) - Vrest(Sidx(1:nD))));
    loadALT  = D / (Lmean + eps);
end


% =========================================================================
function demoAPAlternans()
% DEMOAPALTERNANS  Simulate AP alternans and run full analysis.
    dt   = 0.05;
    tSim = 0:dt:10000;
    BCL  = 300;
    Vm   = zeros(size(tSim)); Vm(1) = -85;
    nStim = floor(tSim(end)/BCL);
    APamp_base = [115, 95];
    APD_base   = [220, 160];
    for s = 1:nStim
        t0  = (s-1)*BCL;
        idx = mod(s-1,2)+1;
        Amp = APamp_base(idx) + 2*randn;
        apd = APD_base(idx)   + 3*randn;
        tau_r = 1.5; tau_d = apd/3;
        tB  = tSim - t0 - 2;
        ap  = Amp*(1-exp(-max(tB,0)/tau_r)).*exp(-max(tB,0)/tau_d) - 85;
        ap(tB<0) = 0; ap = ap + 85;
        mask = tB >= 0 & tB < apd*2;
        Vm(mask) = Vm(mask) + ap(mask);
    end
    Vm = Vm + 0.5*randn(size(Vm));
    analyzeAPAlternans(tSim, Vm, 'PacingRate', 1000/BCL, ...
        'Threshold', -20, 'APD_levels', [30 50 80 90]);
end
