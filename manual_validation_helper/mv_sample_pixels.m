function manifest = mv_sample_pixels(conditioned_file, cam, n_pixels, varargin)
%MV_SAMPLE_PIXELS  Draw an SNR-stratified random pixel sample for manual marking.
%
%   manifest = mv_sample_pixels(conditioned_file, cam, n_pixels)
%   manifest = mv_sample_pixels(..., 'Name', value)
%
%   Builds the pixel list that every reviewer will mark.  Run this ONCE per
%   recording: all reviewers must mark the SAME pixels, or the inter-observer
%   comparison in mv_compare has nothing to pair up.
%
%   The sample is stratified by SNR and drawn with a fixed seed.  Both matter:
%
%     Stratified — a uniform random draw over tissue pixels is dominated by
%       whatever SNR happens to be typical for that heart, so the result says
%       nothing about behaviour at the margins.  Sampling equally from SNR
%       terciles lets mv_compare report agreement AS A FUNCTION OF SNR, which
%       is the operating envelope a user actually needs ("APD80 agrees within
%       5 ms above SNR 8") and pre-empts the obvious criticism that the
%       validation only used clean signals.
%
%     Seeded — the pixel list is a pre-registered choice made before anyone
%       looks at an answer.  Hand-picking pixels, or re-drawing until the
%       agreement looks good, invalidates the whole exercise.
%
%   Inputs
%     conditioned_file  path to a *-conditioned.mat (signal-conditioning output)
%     cam               "CAM1" | "CAM2" | "CAM3" | "CAM4"
%     n_pixels          total pixels to sample (divided evenly across strata)
%
%   Name-value
%     'Window'   [t0 t1] in ms — the beat(s) to mark.  All reviewers mark the
%                same window.  Default: auto from the pacing trace (analog1),
%                taking the beat starting at the 'BeatIndex'-th stimulus.
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
    p.addParameter('Window',    [],  @(v) isempty(v) || (isnumeric(v) && numel(v) == 2));
    p.addParameter('BeatIndex', 5,   @(v) isnumeric(v) && isscalar(v) && v >= 1);
    p.addParameter('NumBeats',  1,   @(v) isnumeric(v) && isscalar(v) && v >= 1);
    p.addParameter('NumStrata', 3,   @(v) isnumeric(v) && isscalar(v) && v >= 1);
    p.addParameter('Seed',      42,  @(v) isnumeric(v) && isscalar(v));
    p.addParameter('Output',    '',  @(v) ischar(v) || isstring(v));
    p.parse(varargin{:});
    o = p.Results;

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

    require_field(d, cam,           conditioned_file);
    require_field(d, [cam '_SNR'],  conditioned_file);
    require_field(d, 'acqFreq',     conditioned_file);

    X   = d.(cam);
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

    % ---- marking window ----
    if ~isempty(o.Window)
        win_ms = double(o.Window(:))';
        win_src = 'user-specified';
    else
        [win_ms, win_src] = window_from_pacing(d, fs, nt, o.BeatIndex, o.NumBeats);
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

    % ---- stratify by SNR, sample evenly ----
    % pctile/randperm rather than quantile/randsample: base MATLAB only, so the
    % harness runs without the Statistics Toolbox.
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

    [row, col] = ind2sub([nr nc], sel);

    manifest = struct();
    manifest.conditioned_file = char(conditioned_file);
    manifest.cam              = cam;
    manifest.acqFreq          = fs;
    manifest.size             = [nr nc nt];
    manifest.window_ms        = win_ms;
    manifest.window_frames    = win_fr;
    manifest.window_source    = win_src;
    manifest.num_beats        = o.NumBeats;
    manifest.seed             = o.Seed;
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
    fprintf('  SNR      : %.1f-%.1f across %d strata (seed %d)\n', ...
            min(manifest.snr), max(manifest.snr), o.NumStrata, o.Seed);
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
