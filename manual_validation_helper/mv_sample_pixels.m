function manifest = mv_sample_pixels(conditioned_file, cam, n_pixels, varargin)
%MV_SAMPLE_PIXELS  Draw the candidate pixels a lead reviewer turns into the pool.
%
%   manifest = mv_sample_pixels(conditioned_file, cam, n_pixels)
%   manifest = mv_sample_pixels(..., 'Sampling', 'grid')
%   manifest = mv_sample_pixels(..., 'Name', value)
%
%   Builds the pixel list for one recording x channel.  Run this ONCE per
%   recording: all reviewers must mark the SAME pixels, or the inter-observer
%   comparison in mv_compare has nothing to pair up.
%
%   CANDIDATES, THEN A POOL
%
%   A pixel can clear an SNR floor and still be unmarkable — a distorted
%   morphology, a motion artefact, a double hump — and the harness must not
%   decide that for itself, because "is this a usable action potential" is a
%   judgement that belongs with the human, not with code that would then be
%   pre-screening the software's own input.  So the draw is oversampled by
%   'Candidates' (default 2x n_pixels) and written with the pool OPEN.  The
%   LEAD reviewer then runs
%
%       mv_mark(manifest_file, 'AB', 'Lead', true)
%
%   which presents the candidates in manifest order and stops when n_pixels
%   have been accepted.  Accepted pixels become the pool; the lead's skips are
%   recorded (and reportable as a skip rate); candidates never reached are
%   discarded.  Every other reviewer, and mv_compare, sees only the pool.
%   Blinding is unchanged: the lead sees no CADENCE value and no pixel position.
%
%   Pass 'Candidates', 1 for the older behaviour, where the draw IS the pool
%   and reviewers skip what they cannot mark.
%
%   SNR FLOOR
%
%   Candidates are drawn from pixels at or above an SNR floor.  The default is
%   the same ADAPTIVE rule Feature Extraction masks its maps with,
%
%       floor = max(1.5, 0.5 * median(snr(snr > 1.5)))
%
%   re-implemented here (one line) rather than imported, so the harness stays
%   independent of the pipeline.  Sampling from the same floor makes the
%   comparison like-for-like: pixels below it have no CADENCE value to compare
%   against and were being dropped in mv_compare anyway.  One consequence to
%   state in Methods: the SNR envelope then spans above-floor pixels only, so
%   the lowest stratum means "usable but dim", not "noise".
%
%   TWO SAMPLING MODES
%
%     'stratified' (default)  Random draw, equal numbers from each SNR
%       tercile, fixed seed.  A uniform random draw over tissue pixels is
%       dominated by whatever SNR happens to be typical for that heart, so it
%       says nothing about behaviour at the margins.  Sampling equally from
%       terciles lets mv_compare report agreement AS A FUNCTION OF SNR, which
%       is the operating envelope a user actually needs ("APD80 agrees within
%       5 ms above SNR 8") and pre-empts the obvious criticism that the
%       validation only used clean signals.  Candidates are interleaved
%       round-robin across strata, so wherever the lead stops the pool stays
%       balanced.
%
%     'grid'  Regular lattice over the tissue.  Costs the same number of marks
%       but places them evenly, so they can be INTERPOLATED INTO A MANUAL MAP.
%       Scattered pixels cannot be, and the manual-vs-automated map panels are
%       half the point of the validation figure.  The same marks then do triple
%       duty: map, agreement statistics, and the SNR envelope.  SNR strata are
%       assigned AFTER the pool is final (post-stratified) rather than driving
%       the draw.  Candidates come in TIERS: tier 1 is the lattice at the
%       chosen stride, tier 2 the lattice offset by half a stride (the cell
%       centres), then the row- and column-offset lattices.  The lead marks
%       tier 1 first and fills from later tiers, so the pool stays evenly
%       spread wherever it stops.
%
%       Honesty rule that goes with this mode: interpolate for the MAP PANELS
%       only.  Every statistic is computed at the marked points, and the
%       difference map is plotted as dots at those points, never interpolated.
%
%   Either way the list is SEEDED and fixed before anyone looks at an answer.
%   Hand-picking pixels, or re-drawing until the agreement looks good,
%   invalidates the whole exercise.
%
%   Inputs
%     conditioned_file  path to a *-conditioned.mat (signal-conditioning output)
%     cam               "CAM1" | "CAM2" | "CAM3" | "CAM4"
%     n_pixels          pixels wanted in the POOL.  Exact for 'stratified' when
%                       the lead reaches it; a TARGET for 'grid', where tier-1
%                       counts are quantized by the lattice (the chosen stride
%                       and actual counts are printed and recorded).
%
%   Name-value
%     'Candidates' oversampling factor (default 2).  Candidates drawn =
%                ceil(Candidates * n_pixels), capped by what the floor allows.
%                1 = the draw is the pool (no lead session).
%     'SNRFloor' 'adaptive' (default) | numeric.  0 reproduces the old
%                "SNR > 0 is tissue" convention.
%     'Sampling' 'stratified' (default) | 'grid'.
%     'Stride'   lattice spacing in pixels for 'grid'.  Default [] = choose the
%                stride whose tier-1 yield comes closest to n_pixels.
%     'MetricsFile'  path to the *-metrics.mat this recording will be compared
%                against.  STRONGLY RECOMMENDED for the raw arm.  Takes the
%                marking window from the analysis window CADENCE actually used,
%                so reviewers mark the same beat the software measured.  Without
%                it the window is chosen independently from the pacing trace,
%                which is fine for DURATIONS (APD, CaTD, rise time) but makes
%                any window-relative TIME (activation, repolarization)
%                impossible to compare at all.  Reads only the window: no metric
%                value is opened, so blinding is unaffected.
%     'Window'   [t0 t1] in ms — the beat(s) to mark.  All reviewers mark the
%                same window.  Default: from 'MetricsFile' if given, else auto
%                from the pacing trace (analog1), taking the beat starting at
%                the 'BeatIndex'-th stimulus.
%     'BeatIndex' which paced beat to use when Window is auto (default 5, to
%                skip the rate transient at the start of the train).
%     'NumBeats' beats to include in the window (default 1; use 2 for alternans).
%     'NumStrata' SNR strata (default 3).
%     'Seed'     RNG seed (default 42).
%     'Output'   manifest path (default alongside the conditioned file).
%
%   Output
%     manifest   struct, also saved to disk, consumed by mv_mark.  The pixel
%                arrays (pixels, linear_index, snr, stratum, tier) list ALL
%                candidates; manifest.active says which are in the pool and
%                manifest.pool_status is 'open' until the lead finalises it.
%
%   Deliberately standalone: this file loads the .mat with plain load() and does
%   its own field checks rather than calling load_cmos, and computes nothing that
%   feature_extraction_helper computes.  See README.md.

    p = inputParser;
    p.addParameter('Source',     'ensemble', @(v) any(strcmpi(v, {'ensemble','raw'})));
    p.addParameter('Sampling',   'stratified', @(v) any(strcmpi(v, {'stratified','grid'})));
    p.addParameter('Candidates', 2,   @(v) isnumeric(v) && isscalar(v) && v >= 1);
    p.addParameter('SNRFloor',   'adaptive', @(v) (isnumeric(v) && isscalar(v) && v >= 0) || ...
                                              ((ischar(v) || isstring(v)) && strcmpi(v, 'adaptive')));
    p.addParameter('Stride',     [],  @(v) isempty(v) || (isnumeric(v) && isscalar(v) && v >= 1));
    p.addParameter('MetricsFile', '', @(v) ischar(v) || isstring(v));
    p.addParameter('Window',     [],  @(v) isempty(v) || (isnumeric(v) && numel(v) == 2));
    p.addParameter('BeatIndex',  5,   @(v) isnumeric(v) && isscalar(v) && v >= 1);
    p.addParameter('NumBeats',   1,   @(v) isnumeric(v) && isscalar(v) && v >= 1);
    p.addParameter('NumStrata',  3,   @(v) isnumeric(v) && isscalar(v) && v >= 1);
    p.addParameter('Seed',       42,  @(v) isnumeric(v) && isscalar(v));
    p.addParameter('Output',     '',  @(v) ischar(v) || isstring(v));
    p.parse(varargin{:});
    o = p.Results;
    o.Sampling = lower(o.Sampling);
    o.Source   = lower(o.Source);

    cam = char(cam);

    % ---- load and check only what this harness needs ----
    S = load(conditioned_file);
    if isfield(S, 'cmos_all_data')
        d = S.cmos_all_data;
    else
        fn = fieldnames(S);
        is_struct = structfun(@isstruct, S);
        if ~any(is_struct)
            error('mv_sample_pixels:noStruct', '%s contains no cmos_all_data struct.', conditioned_file);
        end
        d = S.(fn{find(is_struct, 1)});
    end

    % Which array do the reviewers see?  Feature extraction runs on the
    % ENSEMBLE AVERAGE whenever that checkbox is on, which is the recommended
    % workflow, and every metric except alternans is derived from it.  Marking
    % the same array is what makes the comparison like-for-like: it removes
    % beat-to-beat variation as a source of apparent disagreement, and it gives
    % both sides the SAME time origin, which is the only way a window-relative
    % time (activation, repolarization) can be validated at all.
    if strcmp(o.Source, 'ensemble')
        data_field = [cam '_average'];
        if ~isfield(d, data_field) || isempty(d.(data_field))
            error('mv_sample_pixels:noEnsemble', ...
                  ['%s has no %s. Ensemble averaging was not run in Signal Conditioning. ' ...
                   'Re-condition with averaging on, or pass ''Source'', ''raw'' to mark a ' ...
                   'single beat (durations only — window-relative times will not be ' ...
                   'comparable).'], conditioned_file, data_field);
        end
    else
        data_field = cam;
    end

    require_field(d, data_field,    conditioned_file);
    require_field(d, [cam '_SNR'],  conditioned_file);
    require_field(d, 'acqFreq',     conditioned_file);

    X   = d.(data_field);
    % NB the SNR map is computed on the RAW stack, before averaging (Signal
    % Conditioning runs it first, "since later stages (mask, inversion)"), so
    % it describes single-beat signal quality even when the reviewers are
    % marking a much cleaner averaged trace.  That is the useful reading — it
    % reports agreement as a function of the underlying data quality — but it
    % must be stated that way, not as the SNR of the trace on screen.
    snr = double(d.([cam '_SNR']));
    fs  = double(d.acqFreq);
    if ~isfinite(fs) || fs <= 0
        error('mv_sample_pixels:acqFreq', 'acqFreq is missing or non-finite in %s.', conditioned_file);
    end
    [nr, nc, nt] = size(X);
    if ~isequal(size(snr), [nr nc])
        error('mv_sample_pixels:snrSize', ...
              '%s_SNR is %s but %s is %dx%d — SNR map and stack disagree.', ...
              cam, mat2str(size(snr)), cam, nr, nc);
    end

    % Alternans is the one metric that CANNOT come from the average: the whole
    % point is the beat-to-beat difference, and averaging removes it.
    if strcmp(o.Source, 'ensemble') && o.NumBeats > 1
        error('mv_sample_pixels:ensembleAlternans', ...
              ['NumBeats = %d with Source = ''ensemble''. The ensemble average is a single ' ...
               'representative beat, so alternans cannot be marked on it — averaging is ' ...
               'exactly what removes the alternation. Use ''Source'', ''raw'' for the ' ...
               'alternans arm.'], o.NumBeats);
    end

    % ---- marking window ----
    if ~isempty(o.Window)
        win_ms = double(o.Window(:))';
        win_src = 'user-specified';
    elseif strcmp(o.Source, 'ensemble')
        % The averaged array IS one beat, so there is no beat to choose and no
        % origin to reconcile: frame 1 here is frame 1 for feature extraction.
        win_ms  = [0, (nt - 1) / fs * 1000];
        win_src = 'ensemble-averaged beat (full array)';
    elseif ~isempty(o.MetricsFile)
        [win_ms, win_src] = window_from_metrics(o.MetricsFile, fs);
    else
        [win_ms, win_src] = window_from_pacing(d, fs, nt, o.BeatIndex, o.NumBeats);
        warning('mv_sample_pixels:independentWindow', ...
                ['Marking window chosen from the pacing trace, independently of the window ' ...
                 'CADENCE analysed. Durations (APD, CaTD, rise time) are still comparable, ' ...
                 'but window-relative TIMES (act_times, rep_data) are not. Pass ' ...
                 '''MetricsFile'' to mark the same beat the software measured.']);
    end
    win_fr = round(win_ms / 1000 * fs) + 1;                 % ms -> 1-based frame
    win_fr = [max(1, win_fr(1)), min(nt, win_fr(2))];
    if diff(win_fr) < 3
        error('mv_sample_pixels:window', ...
              'Marking window [%.1f %.1f] ms spans %d frames — too short to mark.', ...
              win_ms(1), win_ms(2), diff(win_fr) + 1);
    end

    % ---- tissue pixels: SNR floor ----
    % Mirrors create_snr_mask's adaptive rule without calling it (see header).
    % No morphological cleanup: that needs the Image Processing Toolbox and
    % would only matter at the mask border, which the lattice/draw rarely hits.
    [floor_snr, floor_mode] = resolve_floor(snr, o.SNRFloor);
    if floor_snr > 0
        tissue = isfinite(snr) & snr >= floor_snr;
    else
        tissue = isfinite(snr) & snr > 0;      % legacy convention
    end
    lin = find(tissue);
    if isempty(lin)
        error('mv_sample_pixels:noTissue', ...
              'No pixels at or above the SNR floor %.2f (%s) in %s_SNR.', floor_snr, floor_mode, cam);
    end

    % Also require the trace to be finite and non-flat inside the marking
    % window: a pixel that is flat there is unmarkable and would only ever
    % produce a skip.
    W  = double(reshape(X, nr * nc, nt));
    W  = W(lin, win_fr(1):win_fr(2));
    ok = all(isfinite(W), 2) & (max(W, [], 2) - min(W, [], 2)) > 0;
    lin = lin(ok);
    if isempty(lin)
        error('mv_sample_pixels:noUsable', ...
              'No tissue pixel has finite, non-flat data inside the marking window.');
    end
    s_ok = snr(lin);

    % ---- draw the candidates ----
    % pctile/randperm rather than quantile/randsample: base MATLAB only, so the
    % harness runs without the Statistics Toolbox.
    n_cand = ceil(o.Candidates * n_pixels);
    rng(o.Seed, 'twister');

    switch o.Sampling
        case 'stratified'
            edges    = [-inf, pctile(s_ok, (1:o.NumStrata-1) / o.NumStrata * 100), inf];
            per      = floor(n_pixels / o.NumStrata);
            per_cand = ceil(o.Candidates * per);

            lists = cell(1, o.NumStrata);
            for k = 1:o.NumStrata
                in_k = find(s_ok > edges(k) & s_ok <= edges(k+1));
                take = min(per_cand, numel(in_k));
                if numel(in_k) < per
                    warning('mv_sample_pixels:thinStratum', ...
                            'SNR stratum %d has only %d pixels; the pool wants %d from it.', ...
                            k, numel(in_k), per);
                end
                lists{k} = lin(in_k(randperm(numel(in_k), take)));
            end
            % Round-robin across strata: the lead marks in this order and may
            % stop anywhere, so any prefix must be as balanced as the whole.
            [sel, stratum] = interleave(lists);
            tier        = nan(numel(sel), 1);
            stride_used = NaN;

        case 'grid'
            usable = false(nr, nc);
            usable(lin) = true;
            % With Candidates = 1 the draw is the pool, so only the base
            % lattice is used (a later tier would break the single-lattice
            % property the map interpolation relies on).
            if o.Candidates > 1; want = n_cand; else; want = 0; end
            [tiers, stride_used] = grid_tiers(usable, n_pixels, o.Stride, want);
            if isempty(tiers{1})
                error('mv_sample_pixels:emptyGrid', ...
                      'No lattice point at stride %d fell on usable tissue.', stride_used);
            end
            % Random order WITHIN each tier: a partly-marked tier is then still
            % spatially spread rather than filled from one corner.
            sel = []; tier = [];
            for t = 1:numel(tiers)
                pts = tiers{t};
                pts = pts(randperm(numel(pts)));
                sel  = [sel;  pts(:)];                       %#ok<AGROW>
                tier = [tier; repmat(t, numel(pts), 1)];     %#ok<AGROW>
            end
            % Strata are assigned by rank over the FINAL pool (mv_finalize_pool),
            % since the lattice fixes where the pixels are, not which SNR they
            % span.  Until then they are undefined.
            stratum = nan(numel(sel), 1);
            edges   = [];
    end

    [row, col] = ind2sub([nr nc], sel);
    n_sel = numel(sel);

    manifest = struct();
    manifest.conditioned_file = char(conditioned_file);
    manifest.cam              = cam;
    manifest.source           = o.Source;      % 'ensemble' | 'raw'
    manifest.data_field       = data_field;    % the array the reviewers see
    manifest.acqFreq          = fs;
    manifest.size             = [nr nc nt];
    manifest.window_ms        = win_ms;
    manifest.window_frames    = win_fr;
    manifest.window_source    = win_src;
    manifest.num_beats        = o.NumBeats;
    manifest.seed             = o.Seed;
    manifest.sampling         = o.Sampling;
    manifest.stride           = stride_used;      % NaN unless sampling == 'grid'
    manifest.num_strata       = o.NumStrata;
    manifest.snr_floor        = floor_snr;
    manifest.snr_floor_mode   = floor_mode;
    manifest.snr_edges        = edges;
    manifest.n_target         = n_pixels;
    manifest.candidates_factor = o.Candidates;
    manifest.pixels           = [row col];        % ALL candidates, n x 2, [row col]
    manifest.linear_index     = sel;
    manifest.snr              = snr(sel);         % positional — order matches pixels
    manifest.stratum          = stratum;
    manifest.tier             = tier;             % grid: 1 = base lattice; NaN otherwise
    manifest.pool_status      = 'open';
    manifest.active           = false(n_sel, 1);
    manifest.lead_reviewer    = '';
    manifest.lead_marks_file  = '';
    manifest.lead_skipped     = zeros(0, 1);
    manifest.pool_finalized   = '';
    manifest.created          = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

    if isempty(o.Output)
        [pdir, base] = fileparts(conditioned_file);
        o.Output = fullfile(pdir, sprintf('%s_%s_manifest.mat', base, cam));
    end
    o.Output = char(o.Output);
    save(o.Output, 'manifest');
    manifest.manifest_file = o.Output;

    fprintf('Drew %d candidate pixels from %s (%s) for a pool of %d\n', ...
            n_sel, cam, manifest.window_source, n_pixels);
    fprintf('  window    : %.1f-%.1f ms (frames %d-%d), %d beat(s)\n', ...
            win_ms(1), win_ms(2), win_fr(1), win_fr(2), o.NumBeats);
    fprintf('  SNR floor : %.2f (%s) -> %d usable pixels\n', floor_snr, floor_mode, numel(lin));
    if strcmp(o.Sampling, 'grid')
        cnt = arrayfun(@(t) sum(tier == t), 1:max(tier));
        fprintf('  sampling  : grid, stride %d px, tiers %s (target %d)\n', ...
                stride_used, mat2str(cnt), n_pixels);
    else
        fprintf('  sampling  : SNR-stratified random (seed %d), round-robin across %d strata\n', ...
                o.Seed, o.NumStrata);
    end
    fprintf('  SNR range : %.1f-%.1f\n', min(manifest.snr), max(manifest.snr));
    fprintf('  manifest  : %s\n', o.Output);

    if o.Candidates > 1
        fprintf('\nPool is OPEN. The lead reviewer builds it:\n');
        fprintf('  mv_mark(''%s'', ''<initials>'', ''Lead'', true)\n', o.Output);
        fprintf('Other reviewers can start once the pool is final. Do not re-run with a\n');
        fprintf('new seed after marking has started.\n');
    else
        % The draw is the pool: finalise it now so downstream tools agree.
        manifest = mv_finalize_pool(o.Output, 'Active', true(n_sel, 1), 'Quiet', true);
        manifest.manifest_file = o.Output;
        fprintf('\nPool is final (Candidates = 1). Give this manifest to every reviewer.\n');
        fprintf('Do not re-run with a new seed after marking has started.\n');
    end
end


% ======================================================================
function require_field(d, f, src)
    if ~isfield(d, f) || isempty(d.(f))
        error('mv_sample_pixels:missingField', ...
              'Required field ''%s'' missing or empty in %s.', f, src);
    end
end


function [floor_snr, mode] = resolve_floor(snr, spec)
%RESOLVE_FLOOR  The SNR floor candidates must clear.
%   'adaptive' reproduces create_snr_mask's per-recording rule: at least half
%   as bright as typical tissue, never below 1.5.  A fixed floor does not
%   transfer across preparations (a dim heart keeps only its brightest corner;
%   a bright one admits border noise), which is why the pipeline uses this
%   rule for the analysis mask and why the harness samples from the same one.
    MIN_FLOOR = 1.5;
    if ischar(spec) || isstring(spec)
        base = isfinite(snr) & snr > MIN_FLOOR;
        if any(base(:))
            floor_snr = max(MIN_FLOOR, 0.5 * median(snr(base)));
        else
            floor_snr = MIN_FLOOR;
        end
        mode = 'adaptive';
    else
        floor_snr = double(spec);
        mode = 'absolute';
    end
end


function [sel, stratum] = interleave(lists)
%INTERLEAVE  Round-robin merge of per-stratum candidate lists.
    n_k = cellfun(@numel, lists);
    sel = zeros(sum(n_k), 1); stratum = zeros(sum(n_k), 1);
    pos = 0;
    for j = 1:max(n_k)
        for k = 1:numel(lists)
            if j <= n_k(k)
                pos = pos + 1;
                sel(pos) = lists{k}(j);
                stratum(pos) = k;
            end
        end
    end
end


function [win_ms, src] = window_from_metrics(metrics_file, fs)
%WINDOW_FROM_METRICS  Take the marking window from CADENCE's analysis window.
%
%   Reads ONLY cmos_all_data.window, the frame range Feature Extraction was run
%   over.  It opens no ep_metrics field and returns no metric value, so this
%   does not weaken the blinding: the window says WHICH BEAT to look at, never
%   where any feature lies inside it.
%
%   Marking the same beat the software measured removes beat-to-beat variation
%   as a source of apparent disagreement, and is the only way a window-relative
%   time (activation, repolarization) can be compared at all.

    S = load(char(metrics_file), 'cmos_all_data');
    if ~isfield(S, 'cmos_all_data')
        error('mv_sample_pixels:badMetrics', ...
              '%s contains no cmos_all_data struct.', metrics_file);
    end
    w = S.cmos_all_data;
    if ~isfield(w, 'window') || isempty(w.window) || ~iscell(w.window) || isempty(w.window{1,1})
        error('mv_sample_pixels:noAnalysisWindow', ...
              ['%s records no analysis window. Feature Extraction was run on the ' ...
               'ENSEMBLE AVERAGE, whose frames do not correspond to any stretch of the ' ...
               'recording — so there is no beat to point the reviewers at, and ' ...
               'window-relative times from that run cannot be validated at all. Either ' ...
               're-extract with a drawn window, or validate durations only and omit ' ...
               '''MetricsFile''.'], metrics_file);
    end

    fr = double(w.window{1,1});
    win_ms = ([fr(1) fr(2)] - 1) / fs * 1000;
    src = sprintf('CADENCE analysis window, frames %d-%d', fr(1), fr(2));
end


function [win_ms, src] = window_from_pacing(d, fs, nt, beat_index, num_beats)
%WINDOW_FROM_PACING  Beat window from the pacing stimulus channel.
%
%   Reads rising edges on analog1 only.  This touches the STIMULUS channel, not
%   the optical signal — it decides which beat everyone marks, never where any
%   feature lies within that beat.  Detecting the beat is not the thing under
%   test; locating fiducials inside it is, and that stays with the human.

    if ~isfield(d, 'analog1') || isempty(d.analog1)
        error('mv_sample_pixels:noPacing', ...
              ['No analog1 pacing trace, so the marking window cannot be chosen ' ...
               'automatically. Pass ''Window'', [t0 t1] in ms explicitly.']);
    end

    a  = double(d.analog1(:));
    mx = max(a); mn = min(a);
    if numel(a) < 10 || mx <= mn
        error('mv_sample_pixels:flatPacing', ...
              'analog1 is flat or too short to find stimuli. Pass ''Window'' explicitly.');
    end

    stim = find(diff(a > (mx + mn) / 2) > 0);   % rising edges
    if numel(stim) < beat_index + num_beats
        error('mv_sample_pixels:tooFewStim', ...
              'analog1 has %d stimuli; need at least %d for BeatIndex %d and %d beat(s).', ...
              numel(stim), beat_index + num_beats, beat_index, num_beats);
    end

    i0 = stim(beat_index);
    i1 = stim(beat_index + num_beats);
    % Start a little before the stimulus so diastole is visible for the baseline
    % click, and stop just short of the next stimulus.
    pad   = round(0.10 * (i1 - i0));
    i0    = max(1,  i0 - pad);
    i1    = min(nt, i1 - 1);
    win_ms = ([i0 i1] - 1) / fs * 1000;
    src    = sprintf('auto from analog1, beat %d', beat_index);
end


function [tiers, stride] = grid_tiers(usable, n_target, stride_req, n_cand)
%GRID_TIERS  Lattice candidates in tiers of decreasing priority.
%
%   Tier 1 is the lattice whose yield comes closest to n_target (or the
%   requested stride).  Later tiers are the same lattice shifted by half a
%   stride — both axes (cell centres), then rows only, then columns only —
%   and are added until the candidate count reaches n_cand.  Because each tier
%   sits between the points of the previous ones, a pool made of tier 1 plus
%   any part of tier 2 is still evenly spread, which is what lets the lead
%   stop at the target and still hand mv_compare something that interpolates
%   into a map.
%
%   Points falling outside the tissue mask are dropped rather than nudged to a
%   neighbour: moving them would bend the lattice and cost the even spacing
%   that is the whole reason for it.

    [rr, cc] = find(usable);
    r0 = min(rr); r1 = max(rr);
    c0 = min(cc); c1 = max(cc);

    if ~isempty(stride_req)
        stride = max(1, round(stride_req));
    else
        % Auto: start from the area estimate, then walk outward and keep
        % whichever stride yields closest to the request.  Yield depends on the
        % mask's shape, not just its area, so the estimate alone is not reliable.
        est      = max(1, round(sqrt(nnz(usable) / max(n_target, 1))));
        stride   = est;
        best_err = inf;
        for s = max(1, est - 4) : (est + 4)
            err = abs(numel(lattice(usable, r0, r1, c0, c1, s, [0 0])) - n_target);
            if err < best_err
                best_err = err;
                stride   = s;
            end
        end
    end

    h = floor(stride / 2);
    offsets = {[0 0], [h h], [h 0], [0 h]};
    if h == 0
        offsets = offsets(1);            % stride 1: nothing lies between points
    end

    tiers = {};
    taken = false(size(usable));
    total = 0;
    for t = 1:numel(offsets)
        pts = lattice(usable, r0, r1, c0, c1, stride, offsets{t});
        pts = pts(~taken(pts));
        taken(pts) = true;
        tiers{end+1} = pts;              %#ok<AGROW>
        total = total + numel(pts);
        if total >= n_cand; break; end
    end
    if total < n_cand && numel(tiers) > 1
        warning('mv_sample_pixels:fewCandidates', ...
                'Lattice tiers yield %d candidates; %d were requested.', total, n_cand);
    end
end


function sel = lattice(usable, r0, r1, c0, c1, stride, offset)
%LATTICE  Grid points inside the tissue bounding box, centred on it.
%   Centring spreads the leftover margin evenly instead of piling it against
%   the bottom-right edge, which would bias coverage away from one corner.
%   'offset' shifts the whole lattice by [drow dcol] pixels.

    rows = r0:stride:r1;
    cols = c0:stride:c1;
    rows = rows + floor(((r1 - r0) - (rows(end) - rows(1))) / 2) + offset(1);
    cols = cols + floor(((c1 - c0) - (cols(end) - cols(1))) / 2) + offset(2);
    rows = rows(rows >= 1 & rows <= size(usable, 1));
    cols = cols(cols >= 1 & cols <= size(usable, 2));

    [RR, CC] = ndgrid(rows, cols);
    lin = sub2ind(size(usable), RR(:), CC(:));
    sel = sort(lin(usable(lin)));
end


function v = pctile(x, p)
%PCTILE  Linear-interpolated percentiles. Base MATLAB (no Statistics Toolbox).
    x = sort(x(isfinite(x)));
    n = numel(x);
    if n == 0
        v = nan(size(p)); return;
    elseif n == 1
        v = repmat(x, size(p)); return;
    end
    pos = (p(:)' / 100) * n + 0.5;               % midpoint convention
    v   = interp1(1:n, x(:)', min(max(pos, 1), n), 'linear');
end
