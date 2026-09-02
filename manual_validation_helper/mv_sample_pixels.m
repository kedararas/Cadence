function manifest = mv_sample_pixels(conditioned_file, cam, n_pixels, varargin)
%MV_SAMPLE_PIXELS  Draw the pixel sample every reviewer will mark.
%
%   manifest = mv_sample_pixels(conditioned_file, cam, n_pixels)
%   manifest = mv_sample_pixels(..., 'Sampling', 'grid')
%   manifest = mv_sample_pixels(..., 'Name', value)
%
%   Builds the pixel list that every reviewer will mark.  Run this ONCE per
%   recording: all reviewers must mark the SAME pixels, or the inter-observer
%   comparison in mv_compare has nothing to pair up.
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
%       validation only used clean signals.
%
%     'grid'  Regular lattice over the tissue.  Costs the same number of marks
%       but places them evenly, so they can be INTERPOLATED INTO A MANUAL MAP.
%       Scattered pixels cannot be, and the manual-vs-automated map panels are
%       half the point of the validation figure.  The same marks then do triple
%       duty: map, agreement statistics, and the SNR envelope.  SNR strata are
%       assigned AFTER selection (post-stratified) rather than driving it.
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
%     n_pixels          pixels to sample.  Exact for 'stratified'; a TARGET for
%                       'grid', where achievable counts are quantized by the
%                       lattice (the chosen stride and actual count are printed
%                       and recorded in the manifest).
%
%   Name-value
%     'Sampling' 'stratified' (default) | 'grid'.
%     'Stride'   lattice spacing in pixels for 'grid'.  Default [] = choose the
%                stride whose yield comes closest to n_pixels.
%     'MetricsFile'  path to the *-metrics.mat this recording will be compared
%                against.  STRONGLY RECOMMENDED.  Takes the marking window from
%                the analysis window CADENCE actually used, so reviewers mark
%                the same beat the software measured.  Without it the window is
%                chosen independently from the pacing trace, which is fine for
%                DURATIONS (APD, CaTD, rise time — beat-to-beat variation adds
%                scatter, not bias) but makes any window-relative TIME
%                (activation, repolarization) impossible to compare at all.
%                Reads only the window: no metric value is opened, so blinding
%                is unaffected.
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
%     manifest   struct, also saved to disk, consumed by mv_mark.
%
%   Deliberately standalone: this file loads the .mat with plain load() and does
%   its own field checks rather than calling load_cmos, and computes nothing that
%   feature_extraction_helper computes.  See README.md.

    p = inputParser;
    p.addParameter('Source',    'ensemble', @(v) any(strcmpi(v, {'ensemble','raw'})));
    p.addParameter('Sampling',  'stratified', @(v) any(strcmpi(v, {'stratified','grid'})));
    p.addParameter('Stride',    [],  @(v) isempty(v) || (isnumeric(v) && isscalar(v) && v >= 1));
    p.addParameter('MetricsFile', '', @(v) ischar(v) || isstring(v));
    p.addParameter('Window',    [],  @(v) isempty(v) || (isnumeric(v) && numel(v) == 2));
    p.addParameter('BeatIndex', 5,   @(v) isnumeric(v) && isscalar(v) && v >= 1);
    p.addParameter('NumBeats',  1,   @(v) isnumeric(v) && isscalar(v) && v >= 1);
    p.addParameter('NumStrata', 3,   @(v) isnumeric(v) && isscalar(v) && v >= 1);
    p.addParameter('Seed',      42,  @(v) isnumeric(v) && isscalar(v));
    p.addParameter('Output',    '',  @(v) ischar(v) || isstring(v));
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

    % ---- tissue pixels ----
    % Tissue is SNR > 0, matching the convention validate_conditioned uses.
    % Also require the trace to be finite and non-flat inside the marking
    % window: a pixel that is flat there is unmarkable and would only ever
    % produce a skip.
    tissue = isfinite(snr) & snr > 0;
    lin    = find(tissue);
    if isempty(lin)
        error('mv_sample_pixels:noTissue', 'No pixels with SNR > 0 in %s_SNR.', cam);
    end

    W  = double(reshape(X, nr * nc, nt));
    W  = W(lin, win_fr(1):win_fr(2));
    ok = all(isfinite(W), 2) & (max(W, [], 2) - min(W, [], 2)) > 0;
    lin = lin(ok);
    if isempty(lin)
        error('mv_sample_pixels:noUsable', ...
              'No tissue pixel has finite, non-flat data inside the marking window.');
    end
    s_ok = snr(lin);

    % ---- select the pixels ----
    % pctile/randperm rather than quantile/randsample: base MATLAB only, so the
    % harness runs without the Statistics Toolbox.
    switch o.Sampling
        case 'stratified'
            edges = [-inf, pctile(s_ok, (1:o.NumStrata-1) / o.NumStrata * 100), inf];
            per   = floor(n_pixels / o.NumStrata);

            rng(o.Seed, 'twister');
            sel = []; stratum = [];
            for k = 1:o.NumStrata
                in_k = find(s_ok > edges(k) & s_ok <= edges(k+1));
                take = min(per, numel(in_k));
                if take < per
                    warning('mv_sample_pixels:thinStratum', ...
                            'SNR stratum %d has only %d pixels; requested %d.', k, numel(in_k), per);
                end
                pick    = in_k(randperm(numel(in_k), take));
                sel     = [sel;     lin(pick)];      %#ok<AGROW>
                stratum = [stratum; repmat(k, take, 1)]; %#ok<AGROW>
            end
            stride_used = NaN;

        case 'grid'
            usable = false(nr, nc);
            usable(lin) = true;
            [sel, stride_used] = grid_sample(usable, n_pixels, o.Stride);
            if isempty(sel)
                error('mv_sample_pixels:emptyGrid', ...
                      'No lattice point at stride %d fell on usable tissue.', stride_used);
            end

            % Post-stratify.  The lattice fixes WHERE the pixels are, so strata
            % are assigned afterwards from the SNR the sample happens to span.
            % mv_compare's by-stratum reporting is unchanged; only the sampling
            % feeding it differs, and post-stratification is the honest way to
            % keep that reporting when coverage, not SNR balance, sets the draw.
            %
            % Assigned by RANK, not by percentile value.  An SNR map with few
            % distinct values (or heavy ties at one level) collapses percentile
            % edges onto each other and silently empties a stratum.  Equal-count
            % groups are what "terciles" means operationally here, and they are
            % what the stratified mode produces by construction, so a stratum
            % index means the same thing in both modes.
            if numel(sel) < o.NumStrata
                error('mv_sample_pixels:tooFewForStrata', ...
                      'Grid yielded %d pixels, fewer than the %d requested strata.', ...
                      numel(sel), o.NumStrata);
            end
            s_sel      = snr(sel);
            [~, ord]   = sort(s_sel, 'ascend');
            cuts       = round(linspace(0, numel(sel), o.NumStrata + 1));
            stratum    = zeros(numel(sel), 1);
            for k = 1:o.NumStrata
                stratum(ord(cuts(k)+1 : cuts(k+1))) = k;
            end
            % Record the boundaries the ranking actually landed on.
            edges = [-inf, arrayfun(@(k) s_sel(ord(cuts(k+1))), 1:o.NumStrata-1), inf];
    end

    [row, col] = ind2sub([nr nc], sel);

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
    manifest.snr_edges        = edges;
    manifest.pixels           = [row col];        % n x 2, [row col]
    manifest.linear_index     = sel;
    manifest.snr              = snr(sel);         % positional — order matches pixels
    manifest.stratum          = stratum;
    manifest.created          = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

    if isempty(o.Output)
        [pdir, base] = fileparts(conditioned_file);
        o.Output = fullfile(pdir, sprintf('%s_%s_manifest.mat', base, cam));
    end
    save(char(o.Output), 'manifest');
    manifest.manifest_file = char(o.Output);

    fprintf('Sampled %d pixels from %s (%s)\n', size(manifest.pixels, 1), cam, manifest.window_source);
    fprintf('  window   : %.1f-%.1f ms (frames %d-%d), %d beat(s)\n', ...
            win_ms(1), win_ms(2), win_fr(1), win_fr(2), o.NumBeats);
    if strcmp(o.Sampling, 'grid')
        fprintf('  sampling : grid, stride %d px (target was %d)\n', stride_used, n_pixels);
        fprintf('  SNR      : %.1f-%.1f, post-stratified into %d strata\n', ...
                min(manifest.snr), max(manifest.snr), o.NumStrata);
    else
        fprintf('  sampling : SNR-stratified random (seed %d)\n', o.Seed);
        fprintf('  SNR      : %.1f-%.1f across %d strata\n', ...
                min(manifest.snr), max(manifest.snr), o.NumStrata);
    end
    fprintf('  manifest : %s\n', o.Output);
    fprintf('\nGive this manifest to every reviewer. Do not re-run with a new seed\n');
    fprintf('after marking has started.\n');
end


% ======================================================================
function require_field(d, f, src)
    if ~isfield(d, f) || isempty(d.(f))
        error('mv_sample_pixels:missingField', ...
              'Required field ''%s'' missing or empty in %s.', f, src);
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


function [sel, stride] = grid_sample(usable, n_target, stride_req)
%GRID_SAMPLE  Regular lattice over the tissue, keeping only usable points.
%
%   Returns linear indices, in raster order.  Points falling outside the tissue
%   mask are dropped rather than nudged to a neighbour: moving them would bend
%   the lattice and cost the even spacing that is the whole reason for it.
%
%   With a square lattice the achievable counts are quantized — 8x8, 9x9, 11x11
%   — so n_target is a target, not a quota.  The caller reports what was
%   actually drawn.

    [rr, cc] = find(usable);
    r0 = min(rr); r1 = max(rr);
    c0 = min(cc); c1 = max(cc);

    if ~isempty(stride_req)
        stride = max(1, round(stride_req));
        sel    = lattice(usable, r0, r1, c0, c1, stride);
        return;
    end

    % Auto: start from the area estimate, then walk outward and keep whichever
    % stride yields closest to the request.  Yield depends on the mask's shape,
    % not just its area, so the closed-form estimate alone is not reliable.
    est         = max(1, round(sqrt(nnz(usable) / max(n_target, 1))));
    sel         = [];
    stride      = est;
    best_err    = inf;
    for s = max(1, est - 4) : (est + 4)
        cand = lattice(usable, r0, r1, c0, c1, s);
        err  = abs(numel(cand) - n_target);
        if err < best_err
            best_err = err;
            sel      = cand;
            stride   = s;
        end
    end
end


function sel = lattice(usable, r0, r1, c0, c1, stride)
%LATTICE  Grid points inside the tissue bounding box, centred on it.
%   Centring spreads the leftover margin evenly instead of piling it against
%   the bottom-right edge, which would bias coverage away from one corner.

    rows = r0:stride:r1;
    cols = c0:stride:c1;
    rows = rows + floor(((r1 - r0) - (rows(end) - rows(1))) / 2);
    cols = cols + floor(((c1 - c0) - (cols(end) - cols(1))) / 2);

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
