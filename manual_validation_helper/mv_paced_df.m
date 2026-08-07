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
%     'Tolerance' ratio tolerance for the 1:1 class (default 0.05)
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
    p.parse(varargin{:});
    o = p.Results;
    o.Cam = char(o.Cam);

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
    [fpath, species, status] = deal(strings(n,1));
    [cl_label, f_pace, df, ratio, ri, oi, nstim] = deal(nan(n,1));
    class_ = strings(n,1);

    t_all = tic;
    for k = 1:n
        f = fullfile(files(k).folder, files(k).name);
        fpath(k)   = string(f);
        species(k) = species_from_path(f);
        cl_label(k)= cl_from_name(files(k).name);
        class_(k)  = "";

        try
            sc = get_trace(f, o);
            if isempty(sc)
                status(k) = "no-" + string(o.Cam); continue;
            end
            if ~isfinite(sc.fs) || sc.fs <= 0
                status(k) = "bad-acqFreq"; continue;
            end

            [fp, ns] = pacing_frequency(sc.analog1, sc.fs);
            nstim(k) = ns;
            if ~isfinite(fp)
                status(k) = "unpaced"; continue;
            end
            f_pace(k) = fp;

            if isempty(o.Band)
                [d_, r_, o_] = cardiacSpectralMetrics(sc.trace(:), sc.fs);
            else
                [d_, r_, o_] = cardiacSpectralMetrics(sc.trace(:), sc.fs, 'Band', o.Band);
            end
            df(k) = d_; ri(k) = r_; oi(k) = o_;
            ratio(k)  = d_ / fp;
            class_(k) = classify_ratio(ratio(k), o.Tolerance);
            status(k) = "ok";

        catch ME
            status(k) = "error:" + string(ME.identifier);
        end

        if mod(k, 10) == 0 || k == n
            fprintf('  %4d/%d  (%.0f s elapsed)\n', k, n, toc(t_all));
        end
    end

    results = table(fpath, species, cl_label, f_pace, df, ratio, class_, ri, oi, nstim, status, ...
        'VariableNames', {'file','species','cl_label_ms','f_pace_hz','df_hz','ratio', ...
                          'class','ri','oi','n_stim','status'});

    report(results);

    if ~isempty(o.Output)
        writetable(results, char(o.Output));
        fprintf('\nWrote %s\n', o.Output);
    end
end


% ======================================================================
function sc = get_trace(f, o)
%GET_TRACE  Spatial-mean trace + stimulus channel, from cache or by strided read.

    key   = sprintf('%s_%s_s%d.mat', matlab.lang.makeValidName(f), o.Cam, o.Stride);
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

    % Tissue-mean trace: average only over pixels finite across ALL frames, so
    % the trace's spatial composition does not change frame to frame.
    P    = reshape(A, [], szc(3));
    good = all(isfinite(P), 2);
    if ~any(good); sc = []; return; end

    sc = struct();
    sc.trace   = mean(P(good, :), 1)';
    sc.npix    = sum(good);
    sc.fs      = double(h5read(f, [gn '/acqFreq']));
    sc.analog1 = get_optional(f, gn, dsets, 'analog1');
    sc.size    = szc;
    sc.source  = f;
    save(cache, 'sc');
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
    nt = size(A, 3);
    P  = reshape(A, [], nt);
    good = all(isfinite(P), 2);
    if ~any(good); sc = []; return; end

    sc = struct();
    sc.trace   = mean(P(good, :), 1)';
    sc.npix    = sum(good);
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


function c = classify_ratio(r, tol)
    if ~isfinite(r);                 c = "n/a";
    elseif abs(r - 1)   <= tol,      c = "1:1";
    elseif abs(r - 0.5) <= tol,      c = "2:1 block";
    elseif abs(r - 2)   <= 2 * tol,  c = "2x harmonic";
    else,                            c = "other";
    end
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
    cls = ["1:1", "2:1 block", "2x harmonic", "other"];
    fprintf('\ncapture classification\n');
    for k = 1:numel(cls)
        nk = sum(Rk.class == cls(k));
        fprintf('   %-12s %4d  (%.1f%%)\n', cls(k), nk, 100 * nk / height(Rk));
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
        fprintf('   within 0.25 Hz   : %.1f%%\n', 100 * mean(abs(err_hz) <= 0.25));
        fprintf('   within 1%%        : %.1f%%\n', 100 * mean(abs(err_pct) <= 1));

        sp = unique(one.species);
        fprintf('\n   by species\n');
        for k = 1:numel(sp)
            s = one.species == sp(k);
            fprintf('     %-8s n=%4d  median |err| %.4f Hz  pacing %.1f-%.1f Hz\n', ...
                    sp(k), sum(s), median(abs(err_hz(s))), ...
                    min(one.f_pace_hz(s)), max(one.f_pace_hz(s)));
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
