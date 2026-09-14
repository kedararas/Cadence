function results = mv_paced_df(root, varargin)
%MV_PACED_DF  Validate dominant frequency against the pacing rate, at scale.
%
%   results = mv_paced_df('/path/to/data')
%   results = mv_paced_df(root, 'Cam', 'CAM1', 'Stride', 4, 'Output', 'df.csv')
%
%   Under S1S1 pacing the tissue's dominant frequency MUST equal the pacing
%   frequency.  The pacing frequency is recorded in the stimulus channel of
%   every file, so every paced recording you own carries its own ground truth —
%   exact, external to the optical signal, and free.  This walks a directory
%   tree, computes DF with the same engine CADENCE ships (cardiacSpectralMetrics)
%   and compares it against the stimulus rate.
%
%   Unlike the marking harness in this folder, this script DOES call CADENCE's
%   own code, and should: the reference here is the stimulus channel, not
%   another algorithm.  What is being validated is the spectral engine; what
%   provides the truth is the pacing hardware.
%
%   PERFORMANCE
%   Conditioned/metrics files run 2-8 GB each, so a full load per file is not
%   viable across thousands of recordings.  These are v7.3 (HDF5) files, so the
%   camera stack is read as a STRIDED SLAB — every Stride-th pixel in each
%   spatial axis, all frames — which cuts a 2 GB read to ~130 MB at Stride 4.
%   The spatial-mean trace is then cached in a small sidecar (~50 KB), so the
%   expensive pass happens once and every later re-run is instant.  Delete the
%   cache directory or pass 'Refresh' to force re-extraction.
%
%   Subsampling is safe here: DF is a property of the rhythm, not of fine
%   spatial detail, and the median of the full per-pixel DF map agrees with the
%   spatial-mean DF to the sample resolution.
%
%   WHAT IS VALIDATED: THE PER-PIXEL ENGINE, AS THE APP RUNS IT
%   The strided stack is masked the way Feature Extraction masks it (adaptive
%   SNR mask from CAM<n>_SNR, background set to NaN) and passed to the 3-D
%   path of cardiacSpectralMetrics -- the same call stage_complexity makes.
%   The reported DF is the median of the per-pixel DF map over the mask,
%   which is what the DF figure and mv_figure_stats report.  Two extra
%   columns say how the map behaved: frac_px_at_pace, the fraction of tissue
%   pixels whose own DF is within tolerance of the pacing rate, and
%   frac_px_floor, the fraction whose peak sat at the lower edge of the auto
%   band.  A recording is "floor-pinned" when more than half its tissue
%   pixels are.
%
%   Until 2026-09-13 the harness used a spatial-MEAN trace instead.  On 21
%   of 40 human slice recordings that mean was dominated by slow drift (the
%   mask kept most of the field because the SNR map was elevated everywhere)
%   and the auto band locked onto 0.6-0.7 Hz with the paced peak visible but
%   outranked.  The tissue-mean trace is still stored (trace) and its DF is
%   reported as df_trace_hz, so the two can be compared.  Cache keys changed
%   with this ('_px'), so the first run after it re-extracts.
%
%   TRUTH WITHOUT A STIMULUS CHANNEL
%   With 'TruthFromName', true a file that has no usable analog1 takes its
%   pacing rate from the cycle length in the filename instead, and
%   truth_source says so per row.  Weaker than the hardware, but justified
%   by the same run: across 103 recordings carrying both, label and stimulus
%   channel agreed every time.  The label never touches the DF calculation,
%   so a mislabelled file reads as an engine failure, the conservative
%   direction.  Off by default.
%
%   CAPTURE, AND WHY NOTHING IS EXCLUDED
%   When the tissue does not follow the drive 1:1 -- 2:1 block, spontaneous
%   rhythm, an arrhythmia -- DF legitimately differs from the pacing rate, and
%   that is a property of the heart, not an error in the software.  Excluding
%   those recordings using a capture test would be circular, because the
%   capture test is itself built on DF.  So nothing is excluded: every paced
%   recording is reported, the DF/pacing ratio is classified, and the agreement
%   statistics are computed on the 1:1 population with the rest reported
%   alongside.  A histogram of the ratio should show a sharp spike at 1.0 with
%   small satellites at 0.5 and 2.0 -- which validates the DF engine and
%   demonstrates capture detection in the same figure.
%
%   Name-value
%     'Cam'       camera field (default 'CAM1')
%     'Stride'    spatial subsample step (default 4)
%     'Pattern'   file pattern (default '*.mat')
%     'CacheDir'  sidecar location (default <root>/.mv_df_cache)
%     'MaxFiles'  cap for a pilot run (default Inf)
%     'Refresh'   re-extract even if cached (default false)
%     'Output'    write a CSV of per-recording rows
%     'Tolerance' ratio tolerance for the capture classes (default 0.05).
%                 Applied as the LARGER of Tolerance x pacing rate and 0.55
%                 FFT bins: below ~5 Hz a 5% window is narrower than one bin
%                 (0.244 Hz at 1 kHz / NFFT 4096), so a correct-bin answer
%                 at 1.33 Hz would otherwise be classed 'other'.
%     'HeartMap'  CSV with columns  file, heart  and optionally  slice
%                 (file = path or basename, with or without '-metrics.mat';
%                 column names are case-insensitive; a title line above the
%                 header is skipped) giving the biological replicate for each
%                 recording.  Without it a filename heuristic (strip the CL /
%                 voltage tokens) is used and marked as such; do not report
%                 heart counts from the heuristic.
%     'Band'      explicit DF search band [lo hi] in Hz.  Default [] uses the
%                 shipped adaptive band, which is what a user actually gets and
%                 is therefore the primary result.  A fixed wide band such as
%                 [0.5 50] is a legitimate SECONDARY run: it is the same band for
%                 every recording, so it tests the spectral estimator without
%                 telling it the answer.  Never set a band bracketing a specific
%                 recording's known pacing rate -- that is circular and would
%                 manufacture the agreement this script exists to measure.
%
%   Output
%     results  table, one row per file: species, labelled CL, pacing frequency,
%              DF, ratio, class, RI, OI, status.

    p = inputParser;
    p.addParameter('Cam',       'CAM1', @(v) ischar(v) || isstring(v));
    p.addParameter('Stride',    4,      @(v) isnumeric(v) && isscalar(v) && v >= 1);
    p.addParameter('Pattern',   '*.mat',@(v) ischar(v) || isstring(v));
    p.addParameter('CacheDir',  '',     @(v) ischar(v) || isstring(v));
    p.addParameter('MaxFiles',  Inf,    @isnumeric);
    p.addParameter('Refresh',   false,  @islogical);
    p.addParameter('Output',    '',     @(v) ischar(v) || isstring(v));
    p.addParameter('Tolerance', 0.05,   @isnumeric);
    p.addParameter('Band',      [],     @(v) isempty(v) || numel(v) == 2);
    p.addParameter('HeartMap',  '',     @(v) ischar(v) || isstring(v));
    p.addParameter('TruthFromName', false, @(v) islogical(v) || isnumeric(v));
    p.parse(varargin{:});
    o = p.Results;
    o.Cam = char(o.Cam);

    % This harness needs the engine's sub-bin refinement fields (dfBin,
    % peakAtEdge).  Check once, loudly, rather than failing on every file with
    % a bare "nonExistentField" -- the usual cause is a stale copy of
    % cardiacSpectralMetrics.m earlier on the path.
    [~, ~, ~, probe] = cardiacSpectralMetrics(sin(2*pi*5*(0:1999)'/1000), 1000);
    if ~isfield(probe, 'dfBin') || ~isfield(probe, 'peakAtEdge')
        error('mv_paced_df:staleEngine', ...
              ['The cardiacSpectralMetrics on the path (%s) predates the sub-bin DF ' ...
               'refinement (no dfBin/peakAtEdge). Update feature_extraction_helper/' ...
               'cardiacSpectralMetrics.m to the current version before running.'], ...
              which('cardiacSpectralMetrics'));
    end

    if isempty(o.CacheDir); o.CacheDir = fullfile(char(root), '.mv_df_cache'); end
    o.CacheDir = char(o.CacheDir);
    if ~isfolder(o.CacheDir); mkdir(o.CacheDir); end

    files = dir(fullfile(char(root), '**', char(o.Pattern)));
    files = files(~[files.isdir]);
    % Skip this harness's own artefacts and obvious non-recordings.  The folder
    % test matters: dir('**') recurses into the sidecar cache, whose files are
    % named after their source and so do not match on name alone.
    drop  = contains({files.name},   {'_manifest', '_marks'}) | ...
            contains({files.folder}, '.mv_df_cache');
    files = files(~drop);
    if numel(files) > o.MaxFiles; files = files(1:o.MaxFiles); end

    fprintf('Scanning %d candidate files under %s\n', numel(files), root);
    fprintf('Camera %s, stride %d, cache %s\n\n', o.Cam, o.Stride, o.CacheDir);

    n = numel(files);
    [fpath, species, status, heart, slice, heart_src, mask_src, truth_src] = deal(strings(n,1));
    [npix, npix_all, df_trace, frac_pace, frac_floor] = deal(nan(n,1));
    [cl_label, f_pace, df, df_bin, bin_hz, err_hz, err_bins, ratio, ri, oi, nstim] = deal(nan(n,1));
    class_ = strings(n,1);
    hmap = load_heart_map(o.HeartMap);

    t_all = tic;
    for k = 1:n
        f = fullfile(files(k).folder, files(k).name);
        fpath(k)   = string(f);
        species(k) = species_from_path(f);
        cl_label(k)= cl_from_name(files(k).name);
        class_(k)  = "";
        [heart(k), slice(k), heart_src(k)] = heart_of(files(k).name, species(k), hmap);

        try
            sc = get_trace(f, o);
            if isempty(sc)
                status(k) = "no-" + string(o.Cam); continue;
            end
            if ~isfinite(sc.fs) || sc.fs <= 0
                status(k) = "bad-acqFreq"; continue;
            end
            npix(k) = sc.npix;
            if isfield(sc, 'npix_all'), npix_all(k) = sc.npix_all; end
            if isfield(sc, 'mask_source'), mask_src(k) = string(sc.mask_source); else, mask_src(k) = "all"; end

            [fp, ns] = pacing_frequency(sc.analog1, sc.fs);
            nstim(k) = ns;
            truth_src(k) = "stimulus";
            if ~isfinite(fp) && o.TruthFromName && isfinite(cl_label(k)) && cl_label(k) > 0
                fp = 1000 / cl_label(k); truth_src(k) = "filename";
            end
            if ~isfinite(fp)
                status(k) = "unpaced"; continue;
            end
            f_pace(k) = fp;

            % Trace DF (the pre-2026-09-13 quantity), kept for comparison.
            if isempty(o.Band)
                [dt_, ~, ~, ~] = cardiacSpectralMetrics(sc.trace(:), sc.fs);
            else
                [dt_, ~, ~, ~] = cardiacSpectralMetrics(sc.trace(:), sc.fs, 'Band', o.Band);
            end
            df_trace(k) = dt_;

            % Per-pixel engine run over the app mask: the validated quantity.
            if ~isfield(sc, 'df_map')
                status(k) = "stale-cache"; continue;
            end
            v  = isfinite(sc.df_map);
            d_ = median(sc.df_map(v)); r_ = median(sc.ri_map(v), 'omitnan'); o_ = median(sc.oi_map(v), 'omitnan');
            df(k) = d_; ri(k) = r_; oi(k) = o_;
            df_bin(k)   = median(sc.dfbin_map(v));
            bin_hz(k)   = sc.bin_hz;
            err_hz(k)   = d_ - fp;
            err_bins(k) = err_hz(k) / bin_hz(k);
            ratio(k)    = d_ / fp;
            w = max(o.Tolerance * fp, 0.55 * bin_hz(k));
            frac_pace(k)  = mean(abs(sc.df_map(v) - fp) <= w);
            frac_floor(k) = mean(sc.edge_map(v) < 0);
            % A peak pinned to the LOWER edge of an auto band is not a
            % measurement: sub-band power (drift, a pump line) outranked the
            % rhythm and the band locked onto it.  Reported separately, never
            % as a capture class -- a 3.3 Hz recording "reading" 0.73 Hz is
            % not 2:1 block.
            if isempty(o.Band) && frac_floor(k) > 0.5
                class_(k) = "floor-pinned";
            else
                class_(k) = classify_ratio(d_, fp, bin_hz(k), o.Tolerance);
            end
            status(k) = "ok";

        catch ME
            status(k) = "error:" + string(ME.identifier);
        end

        if mod(k, 10) == 0 || k == n
            fprintf('  %4d/%d  (%.0f s elapsed)\n', k, n, toc(t_all));
        end
    end

    results = table(fpath, species, heart, slice, heart_src, cl_label, f_pace, truth_src, df, df_bin, bin_hz, ...
                    err_hz, err_bins, ratio, class_, frac_pace, frac_floor, df_trace, ri, oi, nstim, ...
                    npix, npix_all, mask_src, status, ...
        'VariableNames', {'file','species','heart','slice','heart_source','cl_label_ms','f_pace_hz', ...
                          'truth_source','df_hz','df_bin_hz','bin_hz','err_hz','err_bins','ratio', ...
                          'class','frac_px_at_pace','frac_px_floor','df_trace_hz','ri','oi','n_stim', ...
                          'n_pix_tissue','n_pix_all','mask_source','status'});

    report(results);

    if ~isempty(o.Output)
        writetable(results, char(o.Output));
        fprintf('\nWrote %s\n', o.Output);
    end
end


% ======================================================================
function sc = get_trace(f, o)
%GET_TRACE  Spatial-mean trace + stimulus channel, from cache or by strided read.

    key   = sprintf('%s_%s_s%d_px.mat', matlab.lang.makeValidName(f), o.Cam, o.Stride);
    if numel(key) > 200; key = [key(1:150) '_' sprintf('%08x', string_hash(f)) '.mat']; end
    cache = fullfile(o.CacheDir, key);

    if ~o.Refresh && isfile(cache)
        C = load(cache); sc = C.sc; return;
    end

    gn = h5_root_group(f);
    if isempty(gn)
        sc = load_fallback(f, o);
        if ~isempty(sc); save(cache, 'sc'); end
        return;
    end

    dsets = h5_dataset_names(f, gn);
    if ~ismember(o.Cam, dsets); sc = []; return; end

    szc = h5_size(f, [gn '/' o.Cam]);
    if numel(szc) < 3; sc = []; return; end

    st  = o.Stride;
    cnt = [floor(szc(1)/st), floor(szc(2)/st), szc(3)];
    A   = double(h5read(f, [gn '/' o.Cam], [1 1 1], cnt, [st st 1]));

    snr = [];
    if ismember([o.Cam '_SNR'], dsets)
        try
            snr = double(h5read(f, [gn '/' o.Cam '_SNR'], [1 1], cnt(1:2), [st st]));
        catch
            snr = [];
        end
    end

    fs = double(h5read(f, [gn '/acqFreq']));
    sc = tissue_mean(A, snr, fs);
    if isempty(sc); return; end
    sc.fs      = fs;
    sc.analog1 = get_optional(f, gn, dsets, 'analog1');
    sc.size    = szc;
    sc.source  = f;
    save(cache, 'sc');
end


function sc = tissue_mean(A, snr, fs)
%TISSUE_MEAN  Per-pixel DF/RI/OI maps and the tissue-mean trace of a (strided) stack.
%   Masking follows stage_complexity in cadence_extract_features: the
%   adaptive SNR mask (create_snr_mask(snr, [], true) when available, else
%   the same rule inline without cleanup, else every finite pixel), with the
%   background set to NaN, then the 3-D path of cardiacSpectralMetrics.
%   Pixels must be finite across ALL frames to count as tissue.
    nt   = size(A, 3);
    P    = reshape(A, [], nt);
    good = all(isfinite(P), 2);
    if ~any(good); sc = []; return; end

    sc = struct();
    sc.trace_all = mean(P(good, :), 1)';
    sc.npix_all  = sum(good);

    tissue = false(size(good));
    if ~isempty(snr) && numel(snr) == numel(good)
        try
            [m, info] = create_snr_mask(snr, [], true);
            tissue = good & m(:); sc.mask_source = "snr-adaptive"; sc.mask_floor = info.threshold;
        catch
            s = snr(:); base = isfinite(s) & s > 1.5;
            if any(base)
                floor_snr = max(1.5, 0.5 * median(s(base)));
                tissue = good & isfinite(s) & s >= floor_snr;
                sc.mask_source = "snr-adaptive"; sc.mask_floor = floor_snr;
            end
        end
    end
    if ~any(tissue)
        tissue = good; sc.mask_source = "all-finite"; sc.mask_floor = NaN;
    end
    sc.trace = mean(P(tissue, :), 1)';
    sc.npix  = sum(tissue);

    % The engine, exactly as the app calls it: masked stack, 3-D path.
    Am = A; Pm = reshape(Am, [], nt); Pm(~tissue, :) = NaN; Am = reshape(Pm, size(A));
    ws = warning('off', 'all'); c = onCleanup(@() warning(ws));
    [sc.df_map, sc.ri_map, sc.oi_map, sp] = cardiacSpectralMetrics(Am, fs);
    sc.dfbin_map = sp.dfBin;
    sc.edge_map  = sp.peakAtEdge;
    sc.band      = sp.band;
    sc.bin_hz    = sp.f(2) - sp.f(1);
end


function sc = load_fallback(f, o)
%LOAD_FALLBACK  Full load for non-v7.3 files (slow; expected to be rare).
    S = load(f);
    vn = fieldnames(S);
    d  = [];
    for k = 1:numel(vn)
        if isstruct(S.(vn{k})); d = S.(vn{k}); break; end
    end
    if isempty(d) || ~isfield(d, o.Cam); sc = []; return; end

    A  = double(d.(o.Cam));
    A  = A(1:o.Stride:end, 1:o.Stride:end, :);
    snr = [];
    if isfield(d, [o.Cam '_SNR'])
        snr = double(d.([o.Cam '_SNR'])); snr = snr(1:o.Stride:end, 1:o.Stride:end);
    end
    sc = tissue_mean(A, snr, double(d.acqFreq));
    if isempty(sc); return; end
    sc.fs      = double(d.acqFreq);
    if isfield(d, 'analog1'); sc.analog1 = double(d.analog1); else; sc.analog1 = []; end
    sc.size    = size(d.(o.Cam));
    sc.source  = f;
end


function gn = h5_root_group(f)
%H5_ROOT_GROUP  Top-level variable group of a v7.3 MAT file ('/metrics',
%   '/cmos_all_data', ...).  Empty if the file is not HDF5-based.
    gn = '';
    try
        i = h5info(f);
    catch
        return;
    end
    for k = 1:numel(i.Groups)
        if ~startsWith(i.Groups(k).Name, '/#')
            gn = i.Groups(k).Name; return;
        end
    end
end


function names = h5_dataset_names(f, gn)
    names = {};
    try
        g = h5info(f, gn);
        names = {g.Datasets.Name};
    catch
    end
end


function sz = h5_size(f, ds)
    sz = [];
    try
        i = h5info(f, ds);
        sz = i.Dataspace.Size;
    catch
    end
end


function v = get_optional(f, gn, dsets, name)
    if ismember(name, dsets)
        v = double(h5read(f, [gn '/' name]));
    else
        v = [];
    end
end


function [fp, nstim] = pacing_frequency(analog1, fs)
%PACING_FREQUENCY  Stimulus rate from rising edges on the pacing channel.
%   This is the ground truth: it comes from the stimulator, not the optics.

    fp = NaN; nstim = 0;
    if isempty(analog1); return; end
    a  = double(analog1(:));
    mx = max(a); mn = min(a);
    if numel(a) < 10 || mx <= mn; return; end

    stim  = find(diff(a > (mx + mn) / 2) > 0);
    nstim = numel(stim);
    if nstim < 3; return; end

    isi = diff(stim);
    % A pacing train should be near-perfectly regular; a wandering interval
    % means this is not a clean S1S1 protocol and the "known" rate is not known.
    if median(isi) <= 0 || (std(isi) / median(isi)) > 0.05; return; end
    fp = fs / median(isi);
end


function c = classify_ratio(df, fp, bin_hz, tol)
%CLASSIFY_RATIO  Capture class from DF vs pacing rate, tolerance in Hz.
%   Window = max(tol * fp, 0.55 bin).  The bin term is what quantisation can
%   cost even when the right bin wins (half a bin, plus a little for ties);
%   the percentage term takes over once the rate is high enough for it to
%   exceed a bin.
    if ~isfinite(df) || ~isfinite(fp), c = "n/a"; return; end
    w = max(tol * fp, 0.55 * bin_hz);
    if     abs(df - fp)     <= w,      c = "1:1";
    elseif abs(df - fp / 2) <= w,      c = "2:1 block";
    elseif abs(df - 2 * fp) <= 2 * w,  c = "2x harmonic";
    else,                              c = "other";
    end
end


function hmap = load_heart_map(f)
%LOAD_HEART_MAP  file -> "heart|slice" lookup, tolerant of how the CSV was made.
    hmap = [];
    if isempty(f), return; end
    f = char(f);
    assert(isfile(f), 'mv_paced_df:heartMap', 'HeartMap file not found: %s', f);
    L = readlines(f); L = strtrim(L); L = L(L ~= "");
    % A spreadsheet export often carries a title line above the header: skip
    % leading lines that do not look like a header/row (fewer than 2 fields).
    % The header is the first line that names a 'file' column; anything above
    % it (a title such as "hearts,,") is skipped.
    first = find(arrayfun(@(l) any(lower(strtrim(split(l, ","))) == "file"), L), 1);
    assert(~isempty(first), 'mv_paced_df:heartMap', 'HeartMap %s has no header line with a ''file'' column.', f);
    hdr = lower(strtrim(split(L(first), ",")));
    ic = find(hdr == "file", 1); ih = find(hdr == "heart", 1); is = find(hdr == "slice", 1);
    assert(~isempty(ic) && ~isempty(ih), 'mv_paced_df:heartMap', ...
           'HeartMap needs columns ''file'' and ''heart'' (found: %s).', strjoin(hdr, ', '));
    keys = {}; vals = {};
    for r = first + 1:numel(L)
        c = strtrim(split(L(r), ","));
        if numel(c) < max([ic ih]), continue; end
        [~, nm, ext] = fileparts(c(ic));
        k = regexprep(lower(nm + ext), '(-metrics|-conditioned)?\.mat$', '');
        sl = ""; if ~isempty(is) && numel(c) >= is, sl = c(is); end
        keys{end+1} = char(k); vals{end+1} = char(c(ih) + "|" + sl); %#ok<AGROW>
    end
    hmap = containers.Map(keys, vals);
end


function [h, sl, src] = heart_of(name, species, hmap)
%HEART_OF  Biological replicate id (and slice, if mapped) for a recording.
    base = regexprep(lower(name), '(-metrics|-conditioned)?\.mat$', '');
    sl = "";
    if ~isempty(hmap) && isKey(hmap, base)
        v = split(string(hmap(base)), "|");
        h = v(1); if numel(v) > 1, sl = v(2); end
        src = "map"; return;
    end
    % Heuristic: drop the CL / voltage tokens and hope the rest names the heart.
    % Wrong whenever several slices come from one heart, so it is labelled.
    h = regexprep(base, '[_-]?\d+(\.\d+)?ms[_-]?\d*(\.\d+)?[vV]?', '');
    h = regexprep(h, '[_-]?\d+bpm', '');
    h = regexprep(h, '[_-]?\d+(\.\d+)?[vV]\b', '');
    h = species + "|" + string(h);
    src = "heuristic";
end


function s = species_from_path(f)
    l = lower(f);
    if     contains(l, 'rat'),    s = "rat";
    elseif contains(l, 'rabbit'), s = "rabbit";
    elseif contains(l, 'human'),  s = "human";
    elseif contains(l, 'pig') || contains(l, 'swine'), s = "pig";
    else,  s = "unknown";
    end
end


function cl = cl_from_name(name)
%CL_FROM_NAME  Cycle length advertised by the filename, for cross-checking the
%   stimulus-derived rate.  Label/hardware mismatches are worth knowing about.
    cl = NaN;
    t = lower(name);
    m = regexp(t, '(\d+)\s*ms', 'tokens', 'once');
    if ~isempty(m); cl = str2double(m{1}); return; end
    m = regexp(t, '(\d+)\s*bpm', 'tokens', 'once');
    if ~isempty(m); cl = 60000 / str2double(m{1}); return; end
    if ~isempty(regexp(t, '[^a-z0-9]1s[^a-z0-9]', 'once')); cl = 1000; end
end


function h = string_hash(s)
    h = uint32(5381);
    for k = 1:numel(s)
        h = uint32(mod(double(h) * 33 + double(s(k)), 2^32));
    end
end


function report(R)
    ok = R.status == "ok";
    fprintf('\n================ paced DF validation ================\n');
    fprintf('files scanned      : %d\n', height(R));
    fprintf('paced & analysable : %d\n', sum(ok));
    st = unique(R.status(~ok));
    for k = 1:numel(st)
        fprintf('   %-18s %d\n', st(k), sum(R.status == st(k)));
    end
    if ~any(ok); fprintf('\nNothing to summarise.\n'); return; end

    Rk = R(ok, :);
    npin = sum(Rk.class == "floor-pinned");
    if npin > 0
        fprintf('\nnot resolved: %d recording(s) with the peak pinned to the auto-band floor\n', npin);
        fprintf('   (sub-band power outranked the rhythm; excluded from the capture classes below)\n');
        sp = unique(Rk.species(Rk.class == "floor-pinned"));
        for k = 1:numel(sp)
            s = Rk.class == "floor-pinned" & Rk.species == sp(k);
            fprintf('     %-8s %3d   pacing %.1f-%.1f Hz\n', sp(k), sum(s), min(Rk.f_pace_hz(s)), max(Rk.f_pace_hz(s)));
        end
    end
    Rk = Rk(Rk.class ~= "floor-pinned", :);
    cls = ["1:1", "2:1 block", "2x harmonic", "other"];
    fprintf('\ncapture classification (n = %d resolved)\n', height(Rk));
    for k = 1:numel(cls)
        nk = sum(Rk.class == cls(k));
        fprintf('   %-12s %4d  (%.1f%%)\n', cls(k), nk, 100 * nk / max(height(Rk), 1));
    end

    one = Rk(Rk.class == "1:1", :);
    if ~isempty(one)
        err_hz  = one.df_hz - one.f_pace_hz;
        err_pct = 100 * err_hz ./ one.f_pace_hz;
        fprintf('\n1:1 population (n = %d)\n', height(one));
        fprintf('   pacing range     : %.2f - %.2f Hz (CL %.0f - %.0f ms)\n', ...
                min(one.f_pace_hz), max(one.f_pace_hz), ...
                1000/max(one.f_pace_hz), 1000/min(one.f_pace_hz));
        fprintf('   median |error|   : %.4f Hz  (%.3f%%)\n', ...
                median(abs(err_hz)), median(abs(err_pct)));
        fprintf('   95th pct |error| : %.4f Hz\n', prctile_(abs(err_hz), 95));
        fprintf('   bias             : %+.4f Hz\n', mean(err_hz));
        fprintf('   |error| in bins  : median %.2f, 95th pct %.2f  (bin %.4f Hz)\n', ...
                median(abs(one.err_bins)), prctile_(abs(one.err_bins), 95), median(one.bin_hz));
        fprintf('   correct bin      : %.1f%%   (bin peak within half a bin of the pacing rate)\n', ...
                100 * mean(abs(one.df_bin_hz - one.f_pace_hz) <= 0.5 * one.bin_hz + 1e-9));
        fprintf('   within 1%% / 2%%   : %.1f%% / %.1f%%\n', ...
                100 * mean(abs(err_pct) <= 1), 100 * mean(abs(err_pct) <= 2));
        fprintf('   tissue pixels at the pacing rate: median %.0f%% of pixels per recording (min %.0f%%)\n', ...
                100 * median(one.frac_px_at_pace), 100 * min(one.frac_px_at_pace));
        if any(one.truth_source == "filename")
            fprintf('   truth from FILENAME (no stimulus channel) in %d of these recordings\n', sum(one.truth_source == "filename"));
        end

        sp = unique(one.species);
        fprintf('\n   by species (hearts from %s)\n', strjoin(cellstr(unique(one.heart_source)), '+'));
        for k = 1:numel(sp)
            s = one.species == sp(k);
            fprintf('     %-8s n=%4d  hearts %3d  median |err| %.4f Hz  pacing %.1f-%.1f Hz\n', ...
                    sp(k), sum(s), numel(unique(one.heart(s))), median(abs(err_hz(s))), ...
                    min(one.f_pace_hz(s)), max(one.f_pace_hz(s)));
        end
        if any(one.heart_source == "heuristic")
            fprintf('   heart counts above use a FILENAME HEURISTIC: pass ''HeartMap'' before reporting them.\n');
        end
    end

    % Filename label vs stimulus channel — catches mislabelled recordings.
    hasl = ok & isfinite(R.cl_label_ms) & isfinite(R.f_pace_hz);
    if any(hasl)
        lab = R.cl_label_ms(hasl);
        act = 1000 ./ R.f_pace_hz(hasl);
        bad = abs(lab - act) > 0.05 * act;
        fprintf('\nfilename CL vs stimulus channel: %d checked, %d mismatch >5%%\n', ...
                sum(hasl), sum(bad));
        if any(bad)
            fprintf('   (mislabelled recordings — worth fixing before they reach a figure)\n');
            idx = find(hasl); idx = idx(bad);
            for k = 1:min(5, numel(idx))
                [~, nm] = fileparts(R.file(idx(k)));
                fprintf('     %-55s label %.0f ms, actual %.0f ms\n', ...
                        extractBefore(nm + "                    ", 56), ...
                        R.cl_label_ms(idx(k)), 1000 / R.f_pace_hz(idx(k)));
            end
        end
    end
    fprintf('=====================================================\n');
end


function v = prctile_(x, p)
    x = sort(x(isfinite(x)));
    if isempty(x); v = NaN; return; end
    if isscalar(x); v = x; return; end
    pos = (p / 100) * numel(x) + 0.5;
    v   = interp1(1:numel(x), x, min(max(pos, 1), numel(x)), 'linear');
end
