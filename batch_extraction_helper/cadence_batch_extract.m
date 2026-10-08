function T = cadence_batch_extract(input_root, output_root, varargin)
%CADENCE_BATCH_EXTRACT  Re-extract CADENCE features for a whole study tree.
%
%   T = cadence_batch_extract(input_root, output_root)
%   T = cadence_batch_extract(input_root, output_root, Name, Value, ...)
%
%   Walks every sub-folder of <input_root> (any depth) for conditioned
%   recordings (*.mat, skipping *-metrics.mat and dot-files), runs the
%   headless Feature Extraction pipeline (cadence_extract_features) on each,
%   writes  <output_root>/<same relative folder>/<name>-metrics.mat  (same
%   layout the app's SAVE DATA writes, so Signal Analysis / Conduction
%   Velocity load them), computes the per-recording metric medians
%   (cadence_recording_medians) and appends them to a CSV after EVERY file,
%   so a crash or a stopped session loses nothing.  Finally builds the Excel
%   summary (cadence_build_summary).
%
%   Name/Value options
%     'Folders'      cellstr of relative folders to process, e.g.
%                    {'ZT2_8am/R2'} or {'ZT2_8am'} (a folder selects
%                    everything below it).  Default {} = the whole tree.
%     'ArrhythmiaPattern' regexp on the file name that marks an arrhythmia
%                    recording.  A HINT only: every recording is classified
%                    from its data (cadence_extract_features, RHYTHM) and gets
%                    the paced or the arrhythmia metric set; the tag decides
%                    only a regular unpaced recording, and a tag that disagrees
%                    with the data is noted in the row.  Default
%                    '(?i)a[r]{1,2}[rh]?h?y+t?h?m', deliberately loose because
%                    the corpus contains misspellings ("Arhhythmia"); it does
%                    not match "alternans" or "EAD".
%     'RhythmOverrides' CSV with columns File and Rhythm (paced | arrhythmia):
%                    recordings whose file name (with or without .mat) is
%                    listed skip the classification.  Default
%                    <output_root>/cadence_rhythm_overrides.csv, if present.
%     'LightsOnHour' clock hour of lights-on, default 6.  When the label has no
%                    explicit ZT (no ZT<n> folder or ZT<n> in the file name), ZT
%                    is derived from the clock in the file label: ZT = clock
%                    hour - LightsOnHour (mod 24), e.g. 8AM -> ZT2, 8PM -> ZT14.
%                    ZT_source says which.  Acquired is the labelled date and
%                    clock ("2025-11-10 08:00").  Both are (re)filled for rows
%                    already in the CSVs, from their file names.
%                    Acquired_camera is the camera's own time stamp
%                    ("2025-12-08 20:45:39"), read from the raw header by
%                    cadence_batch_pipeline (job field acquired_at; header
%                    only, no conversion) and filled for existing rows too.
%     'ExcludePattern' regexp on the file name; matching recordings are NOT
%                    extracted but are listed (ZT, experiment, condition, CL,
%                    tag) in <output_root>/cadence_excluded_recordings.csv.
%                    Default '' (extract everything; arrhythmia recordings
%                    get the arrhythmia metric set).
%     'NormalizeZT'  true (default): an input folder named ZT<n>_<anything>
%                    is written as ZT<n> in the output tree (matches the
%                    lab's existing Metrics layout).  false: mirror names.
%     'Files'        cellstr of file-name patterns (regexp) to keep, default {}.
%     'Resume'       true (default): skip recordings already listed in the
%                    medians CSV whose metrics file exists.  A recording whose
%                    metrics file exists but has NO row (the CSV was lost or
%                    truncated) is ADOPTED: its medians are recomputed from the
%                    saved metrics file instead of extracting it again.
%                    false: redo all.
%
%   DISK SPACE.  The free space is checked before each metrics file is saved;
%   when it is short the run STOPS with one message ('cadence:diskFull')
%   instead of failing every remaining recording, and the CSVs are written
%   through cadence_write_table so a full disk cannot truncate them.
%     'DryRun'       list what would be processed and return (default false).
%     'SaveMetrics'  write the *-metrics.mat files (default true).
%     'ExtractOpts'  struct passed to cadence_extract_features (FOV_mm, ...).
%     'MediansFile'  CSV of per-recording medians
%                    (default <output_root>/cadence_recording_medians.csv)
%     'SummaryFile'  Excel workbook
%                    (default <output_root>/CADENCE_median_summary.xlsx)
%     'CameraFile'   CSV of per-camera medians, one row per camera of every
%                    recording (default <output_root>/cadence_camera_medians.csv;
%                    see cadence_camera_medians).  The workbook gets a
%                    'Per camera' sheet when a recording has more than one
%                    camera of the same signal type.
%     'BuildSummary' rebuild the workbook at the end (default true).
%     'LogFile'      text log (default <output_root>/cadence_batch_log.txt)
%     'MaxFiles'     stop after this many new recordings (default Inf).
%     'Jobs'         struct array of recordings to process instead of walking
%                    <input_root> (fields as built by the local discover:
%                    source, file, zt_folder, rat_folder, zt_out, out_dir,
%                    metrics_file, rel_folder; extra fields are allowed).
%                    Used by cadence_batch_pipeline.  Folders/Files are then
%                    ignored (the caller has already selected).
%     'Loader'       function handle  d = Loader(job, logf)  returning the
%                    conditioned cmos_all_data for a job, in place of
%                    load(job.source).  Lets a caller convert and condition
%                    on demand (cadence_batch_pipeline).  Default: load.
%     'LogFcn'       function handle LogFcn(line) called with every console
%                    line in addition to stdout / the log file (a UI console).
%     'ShouldStop'   function handle returning true to stop before the next
%                    recording (a UI STOP button).  Everything done so far is
%                    kept; the summary is still built.
%     'SaveMapsPDF'  true: also write <output_root>/<rel>/<name>-maps.pdf, the
%                    cardiac maps of every feature and camera (cadence_maps_pdf,
%                    the app's map windows as PDF pages).  Recordings that are
%                    skipped by Resume but have no PDF yet get one rendered
%                    from their saved metrics file, without re-extraction.
%                    Default false.
%     'ExpectFilter' medians 'Filter' label the run is meant to produce, e.g.
%                    "hybrid 100/50 Hz" or "50 Hz" (cadence_batch_pipeline
%                    sets it).  Resume then redoes recordings extracted with a
%                    different filter instead of skipping them; rows from
%                    before the Filter column count as single-band.  Default
%                    '' (no check).
%
%   Returns the medians table (one row per recording).
%
%   Example -- test on one rat, then the whole study:
%     cadence_batch_paths();
%     cadence_batch_extract(in, out, 'Folders', {'ZT2_8am/R2'});
%     cadence_batch_extract(in, out);            % resumes, skips R2

    p = inputParser;
    p.addParameter('Folders', {}, @(x) iscellstr(x) || isstring(x) || ischar(x));
    p.addParameter('Files', {}, @(x) iscellstr(x) || isstring(x) || ischar(x));
    p.addParameter('Resume', true, @islogical);
    p.addParameter('DryRun', false, @islogical);
    p.addParameter('SaveMetrics', true, @islogical);
    p.addParameter('ExtractOpts', struct(), @isstruct);
    p.addParameter('MediansFile', '', @(x) ischar(x) || isstring(x));
    p.addParameter('SummaryFile', '', @(x) ischar(x) || isstring(x));
    p.addParameter('CameraFile', '', @(x) ischar(x) || isstring(x));
    p.addParameter('BuildSummary', true, @islogical);
    p.addParameter('LogFile', '', @(x) ischar(x) || isstring(x));
    p.addParameter('MaxFiles', Inf, @isnumeric);
    p.addParameter('NormalizeZT', true, @islogical);
    p.addParameter('ExcludePattern', '', @(x) ischar(x) || isstring(x));
    p.addParameter('ArrhythmiaPattern', '(?i)a[r]{1,2}[rh]?h?y+t?h?m', @(x) ischar(x) || isstring(x));
    p.addParameter('RhythmOverrides', '', @(x) ischar(x) || isstring(x));
    p.addParameter('LightsOnHour', 6, @isnumeric);
    p.addParameter('Jobs', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Loader', [], @(x) isempty(x) || isa(x, 'function_handle'));
    p.addParameter('LogFcn', [], @(x) isempty(x) || isa(x, 'function_handle'));
    p.addParameter('ShouldStop', [], @(x) isempty(x) || isa(x, 'function_handle'));
    p.addParameter('SaveMapsPDF', false, @islogical);
    p.addParameter('ExpectFilter', '', @(x) ischar(x) || isstring(x));
    p.parse(varargin{:});
    o = p.Results;
    o.Folders = cellstr(o.Folders);  o.Files = cellstr(o.Files);
    o.Folders = o.Folders(~cellfun(@isempty, o.Folders));
    o.Files   = o.Files(~cellfun(@isempty, o.Files));

    input_root  = char(input_root);
    output_root = char(output_root);
    note = @(s) log_line(-1, s, o.LogFcn);      % console line before the log file is open
    if isempty(o.MediansFile), o.MediansFile = fullfile(output_root, 'cadence_recording_medians.csv'); end
    if isempty(o.SummaryFile), o.SummaryFile = fullfile(output_root, 'CADENCE_median_summary.xlsx'); end
    if isempty(o.CameraFile),  o.CameraFile  = fullfile(output_root, 'cadence_camera_medians.csv'); end
    if isempty(o.LogFile),     o.LogFile     = fullfile(output_root, 'cadence_batch_log.txt'); end
    if isempty(o.RhythmOverrides), o.RhythmOverrides = fullfile(output_root, 'cadence_rhythm_overrides.csv'); end
    overrides = read_overrides(o.RhythmOverrides);
    if ~isempty(overrides)
        note(sprintf('%d rhythm override(s) from %s', numel(overrides.file), o.RhythmOverrides));
    end

    if ~exist('compute_lat_50', 'file')
        cadence_batch_paths();
    end
    if isempty(o.Jobs) && ~isfolder(input_root)
        error('cadence_batch_extract:input',  'Input root not found: %s', input_root);
    end
    if ~isfolder(output_root), mkdir(output_root); end

    % ---- discover recordings --------------------------------------------------
    if isempty(o.Jobs)
        jobs = discover(input_root, output_root, o.Folders, o.Files, o.NormalizeZT);
        note(sprintf('%d recording(s) found under %s', numel(jobs), input_root));
    else
        jobs = o.Jobs(:)';
        need = {'source', 'file', 'zt_folder', 'rat_folder', 'zt_out', 'out_dir', 'metrics_file', 'rel_folder'};
        miss = need(~isfield(jobs, need));
        if ~isempty(miss)
            error('cadence_batch_extract:jobs', 'Jobs are missing field(s): %s', strjoin(miss, ', '));
        end
        note(sprintf('%d recording(s) supplied by the caller', numel(jobs)));
    end
    if isempty(jobs), T = table(); return; end

    % ---- excluded recordings (arrhythmia etc.): listed, not extracted -----------------
    o.ExcludePattern = char(o.ExcludePattern);
    if ~isempty(o.ExcludePattern)
        ex = arrayfun(@(j) ~isempty(regexp(j.file, o.ExcludePattern, 'once')), jobs);
        if any(ex)
            write_excluded(jobs(ex), fullfile(output_root, 'cadence_excluded_recordings.csv'));
            note(sprintf('%d recording(s) match ''%s'' -> listed in cadence_excluded_recordings.csv, not extracted', ...
                nnz(ex), o.ExcludePattern));
            jobs = jobs(~ex);
        end
    end

    % ---- resume state ----------------------------------------------------------
    % The medians CSV is always read so that a partial re-run ('Folders',
    % 'Resume' false) replaces its own rows and keeps everyone else's; Resume
    % only decides whether recordings already listed are skipped.
    rows = {};
    if isfile(o.MediansFile)
        try
            [sn, st] = row_schema();
            rows = table_to_rows(cadence_read_table(o.MediansFile, sn(strcmp(st, 'string'))));
            if o.Resume
                note(sprintf('Resuming: %d recording(s) already in %s', numel(rows), o.MediansFile));
            else
                note(sprintf('%d recording(s) already in %s (kept; selected ones are redone)', numel(rows), o.MediansFile));
            end
        catch ME
            warning('cadence_batch_extract:resume', 'Could not read %s (%s); starting fresh.', o.MediansFile, ME.message);
            rows = {};
        end
    end
    Tcam = read_camera_table(o.CameraFile);
    [rows, Tcam] = refresh_meta(rows, Tcam, o.LightsOnHour);   % ZT / clock / acquired from the file label
    [rows, Tcam] = refresh_camera_time(rows, Tcam, jobs);      % camera time stamp, when the caller read it
    done_sources = cellfun(@(r) char(r.Source), rows, 'UniformOutput', false);
    done_metrics = cellfun(@(r) char(r.Metrics_file), rows, 'UniformOutput', false);
    find_row = @(job) find(strcmp(done_sources, job.source) | strcmp(done_metrics, job.metrics_file), 1);

    % With ExpectFilter set (the pipeline sets it), a recording extracted with a
    % different temporal filter is redone rather than skipped, so a resumed run
    % never mixes filters.  Rows from before the Filter column are single-band
    % (filter unknown): kept for a single-band run, redone for a hybrid one.
    expect = string(o.ExpectFilter);
    filter_ok = @(f) expect == "" || f == expect || (f == "" && ~startsWith(expect, "hybrid"));
    todo = false(numel(jobs), 1);  refilter = 0;
    adopt = false(numel(jobs), 1);        % metrics file on disk but no row: summarize it, do not re-extract
    for j = 1:numel(jobs)
        r = find_row(jobs(j));
        already = ~isempty(r) && isfile(jobs(j).metrics_file);
        if already && ~filter_ok(rows{r}.Filter)
            already = false;  refilter = refilter + 1;
        end
        todo(j) = ~(o.Resume && already);
        adopt(j) = o.Resume && isempty(r) && isfile(jobs(j).metrics_file);
    end
    note(sprintf('%d to process, %d skipped (already done)', nnz(todo), nnz(~todo)));
    if any(adopt)
        note(sprintf(['%d of them already have a metrics file but no row in the medians CSV: ' ...
                      'their medians are read from the saved file (no re-extraction).'], nnz(adopt)));
    end
    if o.Resume && refilter > 0
        note(sprintf('%d of them were extracted with a different filter than %s and are redone', refilter, expect));
    end

    if o.DryRun
        for j = find(todo)'
            note(sprintf('  %-24s %s', jobs(j).rel_folder, jobs(j).source));
            note(sprintf('     -> %s', jobs(j).metrics_file));
        end
        T = rows_to_table(rows);
        return;
    end

    % ---- main loop -------------------------------------------------------------
    fid = fopen(o.LogFile, 'a');
    logf = @(s) log_line(fid, s, o.LogFcn);
    logf(sprintf('==== cadence_batch_extract started %s ====', datestr(now, 'yyyy-mm-dd HH:MM:SS')));
    logf(sprintf('input  : %s', input_root));
    logf(sprintf('output : %s', output_root));

    % ---- maps PDF for recordings already extracted (no re-extraction) -----------
    if o.SaveMapsPDF
        for j = find(~todo)'
            job = jobs(j);
            pdf = maps_pdf_name(job.metrics_file);
            if isfile(pdf) || ~isfile(job.metrics_file), continue; end
            if ~isempty(o.ShouldStop) && o.ShouldStop(), break; end
            t0 = tic;
            try
                S = load(job.metrics_file);
                if isfield(S, 'cmos_all_data'), d = S.cmos_all_data;
                else, fn = fieldnames(S); d = S.(fn{1}); end
                clear S
                n = cadence_maps_pdf(d, pdf, 'Name', job.file, 'Log', @(s) logf(['   ' s]));
                clear d
                logf(sprintf('[maps only] %s: %d page(s) in %.0f s', pdf, n, toc(t0)));
                k = find_row(job);
                if ~isempty(k) && n > 0
                    rows{k}.Maps_pdf = string(pdf);
                    cadence_write_table(rows_to_table(rows), o.MediansFile, logf);
                end
            catch ME
                logf(sprintf('[maps only] %s FAILED: %s', job.metrics_file, ME.message));
            end
        end
    end

    n_new = 0;  idx = find(todo)';
    for jj = 1:numel(idx)
        j = idx(jj);
        job = jobs(j);
        if n_new >= o.MaxFiles, break; end
        if ~isempty(o.ShouldStop) && o.ShouldStop()
            logf('STOP requested: stopping before the next recording (done so far is kept).');
            break;
        end
        t_all = tic;
        logf(sprintf('[%d/%d] %s', jj, numel(idx), job.source));

        row = new_row(job, o.LightsOnHour);
        Tc = read_camera_table('');
        row.Extracted_on = string(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        row.Name_tag = double(~isempty(o.ArrhythmiaPattern) && ~isempty(regexp(job.file, o.ArrhythmiaPattern, 'once')));
        disk_full = false;
        if adopt(j)
            [done, row, Tc] = adopt_metrics(job, row, filter_ok, logf);
            if done
                row.Elapsed_s = round(toc(t_all));
                [rows, done_sources, done_metrics, Tcam] = keep_row(rows, done_sources, done_metrics, Tcam, row, Tc, job, find_row);
                cadence_write_table(rows_to_table(rows), o.MediansFile, logf);
                cadence_write_table(Tcam, o.CameraFile, logf);
                n_new = n_new + 1;
                continue;
            end
            row = new_row(job, o.LightsOnHour);                 % could not be adopted: extract it
            Tc = read_camera_table('');
            row.Extracted_on = string(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        end
        try
            % load (or build, when the caller supplied a Loader) + schema gate
            t0 = tic;
            if isempty(o.Loader)
                S = load(job.source);
                if isfield(S, 'cmos_all_data'), d = S.cmos_all_data;
                else, fn = fieldnames(S); d = S.(fn{1}); end
                clear S
            else
                d = o.Loader(job, logf);
            end
            [ok, msgs] = schema_gate(d, "conditioned", job.file);
            for k = 1:numel(msgs), logf(['   ' msgs{k}]); end
            if ~ok, error('cadence_batch_extract:schema', 'schema FAIL'); end
            logf(sprintf('   loaded in %.0f s (%d cams, %d frames, %.0f Hz)', toc(t0), d.num_files, size(d.CAM1,3), d.acqFreq));

            % extract
            eo = o.ExtractOpts;  eo.Log = @(s) logf(['   ' s]);
            eo.NameTag = ~isempty(o.ArrhythmiaPattern) && ~isempty(regexp(job.file, o.ArrhythmiaPattern, 'once'));
            ov = override_for(overrides, job.file);
            if ~isempty(ov), eo.Rhythm = ov; logf(sprintf('   rhythm override: %s', ov)); end
            [d, st] = cadence_extract_features(d, eo);
            row.Features = string(strjoin(st.features, ' '));
            if isfield(d, 'rhythm')
                row.Rhythm = string(d.rhythm.class);  row.Capture = string(d.rhythm.capture);
                row.Rhythm_basis = string(d.rhythm.basis);
                if ~isempty(d.rhythm.note), row.Rhythm_basis = row.Rhythm_basis + " | NOTE: " + string(d.rhythm.note); end
            end
            row.Name_tag = double(eo.NameTag);
            row.Filter   = filter_label(d);
            row.Errors   = string(strjoin(st.errors, ' | '));

            % save (write to a temp name, then rename, so a killed run never
            % leaves a half-written file with the final name)
            if o.SaveMetrics
                t0 = tic;
                if ~isfolder(job.out_dir), mkdir(job.out_dir); end
                tmp = [job.metrics_file '.part'];
                w = whos('d');                       % v7.3 compresses to roughly half of this
                cadence_require_space(job.out_dir, 0.7 * w.bytes + 0.5e9, 'metrics file');
                cmos_all_data = d; %#ok<NASGU>
                try
                    save(tmp, 'cmos_all_data', '-v7.3');
                catch ME
                    if isfile(tmp), delete(tmp); end      % never leave a truncated .part behind
                    free = cadence_require_space(job.out_dir);
                    if free < 1e9
                        error('cadence:diskFull', 'The disk is full (%.2f GB free on the volume holding %s): %s', ...
                              free / 1e9, job.out_dir, ME.message);
                    end
                    rethrow(ME);
                end
                clear cmos_all_data
                if isfile(job.metrics_file), delete(job.metrics_file); end
                movefile(tmp, job.metrics_file);
                logf(sprintf('   saved %s in %.0f s', job.metrics_file, toc(t0)));
                row.Metrics_file = string(job.metrics_file);
            end

            % cardiac maps PDF (after the metrics are safe on disk)
            if o.SaveMapsPDF
                t0 = tic;
                pdf = maps_pdf_name(job.metrics_file);
                try
                    n = cadence_maps_pdf(d, pdf, 'Name', job.file, 'Log', @(s) logf(['   ' s]));
                    if n > 0, row.Maps_pdf = string(pdf); end
                    logf(sprintf('   maps PDF: %d page(s) in %.0f s', n, toc(t0)));
                catch ME
                    logf(sprintf('   maps PDF FAILED: %s', ME.message));
                    if strlength(row.Errors) == 0, row.Errors = string(['maps PDF: ' ME.message]);
                    else, row.Errors = string([char(row.Errors) ' | maps PDF: ' ME.message]); end
                end
            end

            % medians: per recording (first voltage / first calcium camera) ...
            vals = cadence_recording_medians(d);
            fn = fieldnames(vals);
            for k = 1:numel(fn), row.(fn{k}) = vals.(fn{k}); end
            % ... and per camera, so every camera of a multi-camera rig is kept
            try
                Tc = cadence_camera_medians(d, 'Meta', camera_meta(row));
            catch ME
                logf(sprintf('   per-camera medians FAILED: %s', ME.message));
            end
            clear d
        catch ME
            row.Errors = string(sprintf('%s%s', char(row.Errors), [' FATAL: ' ME.message]));
            logf(sprintf('   FAILED: %s', ME.message));
            disk_full = strcmp(ME.identifier, 'cadence:diskFull');
        end
        if disk_full
            % Nothing was saved for this recording and its row is left as it
            % was, so the next run picks up exactly here.
            logf('DISK FULL: stopping the run. Free space on the output volume and run again; Resume continues from this recording.');
            break;
        end
        row.Elapsed_s = round(toc(t_all));
        logf(sprintf('   done in %.0f s', row.Elapsed_s));

        % replace or append, then persist
        [rows, done_sources, done_metrics, Tcam] = keep_row(rows, done_sources, done_metrics, Tcam, row, Tc, job, find_row);
        T = rows_to_table(rows);
        cadence_write_table(T, o.MediansFile, logf);
        cadence_write_table(Tcam, o.CameraFile, logf);
        n_new = n_new + 1;
    end

    T = rows_to_table(rows);
    if ~isempty(rows), cadence_write_table(T, o.MediansFile, logf); end
    if height(Tcam) > 0, cadence_write_table(Tcam, o.CameraFile, logf); end   % refreshed label / time columns
    have_cam = ismember(T.Source, Tcam.Source) | ismember(T.Metrics_file, Tcam.Metrics_file);
    if any(~have_cam)
        logf(sprintf(['%d recording(s) in the medians CSV have no per-camera rows (extracted before ' ...
                      'per-camera medians existed); cadence_recompute_medians fills them.'], nnz(~have_cam)));
    end
    if o.BuildSummary && ~isempty(rows)
        try
            cadence_build_summary(T, o.SummaryFile, 'PerCamera', Tcam);
            logf(sprintf('summary written: %s', o.SummaryFile));
        catch ME
            logf(sprintf('summary FAILED: %s', ME.message));
        end
    end
    logf(sprintf('==== finished %s: %d processed ====', datestr(now, 'yyyy-mm-dd HH:MM:SS'), n_new));
    fclose(fid);
end


% =========================================================================
function jobs = discover(input_root, output_root, folders, file_pats, normalize_zt)
% Recursive walk of input_root.  Each *.mat that is not a metrics / manifest
% / marks file and not a dot-file becomes a job.  The output folder mirrors
% the relative folder of the source (ZT<n>_xxx -> ZT<n> when normalize_zt).
    jobs = struct('source', {}, 'file', {}, 'zt_folder', {}, 'rat_folder', {}, ...
                  'zt_out', {}, 'out_dir', {}, 'metrics_file', {}, 'rel_folder', {});
    folders = strrep(folders, '\', '/');
    folders = regexprep(folders, '/+$', '');

    all_files = dir(fullfile(input_root, '**', '*.mat'));
    all_files = all_files(~[all_files.isdir]);
    for c = 1:numel(all_files)
        nm = all_files(c).name;
        if startsWith(nm, '.'), continue; end
        if ~isempty(regexp(nm, '(-metrics|_manifest|_marks)\.mat$', 'once')), continue; end
        rel = relative_folder(input_root, all_files(c).folder);
        if any(cellfun(@(q) startsWith(q, '.'), strsplit(rel, '/'))), continue; end   % hidden dirs
        if ~isempty(folders)
            hit = false;
            for f = 1:numel(folders)
                if strcmp(rel, folders{f}) || startsWith(rel, [folders{f} '/']), hit = true; break; end
            end
            if ~hit, continue; end
        end
        if ~isempty(file_pats) && ~any(cellfun(@(pt) ~isempty(regexp(nm, pt, 'once')), file_pats))
            continue;
        end
        m = cadence_parse_recording_name(fullfile(all_files(c).folder, nm));
        out_rel = rel;
        if normalize_zt
            out_rel = regexprep(out_rel, '(^|/)(ZT\d+)[^/]*', '$1$2');
        end
        parts = strsplit(rel, '/'); parts = parts(~cellfun(@isempty, parts));
        j = struct();
        j.source       = fullfile(all_files(c).folder, nm);
        j.file         = nm;
        j.zt_folder    = char(m.zt_folder);
        if isempty(parts), j.rat_folder = ''; else, j.rat_folder = parts{end}; end
        if isnan(m.ZT), j.zt_out = ''; else, j.zt_out = sprintf('ZT%d', m.ZT); end
        j.out_dir      = fullfile(output_root, out_rel);
        j.metrics_file = fullfile(j.out_dir, [m.base '-metrics.mat']);
        j.rel_folder   = rel;
        jobs(end+1) = j; %#ok<AGROW>
    end
    % order: ZT number, experiment number, relative folder, run number
    if ~isempty(jobs)
        key = zeros(numel(jobs), 3);  relk = cell(numel(jobs), 1);
        for j = 1:numel(jobs)
            m = cadence_parse_recording_name(jobs(j).source);
            en = regexp(char(m.experiment), '\d+', 'match', 'once');
            key(j,:) = [nz(m.ZT), nz(str2double(en)), nz(m.run)];
            relk{j} = jobs(j).rel_folder;
        end
        [~, ord] = sortrows([num2cell(key(:,1)), num2cell(key(:,2)), relk, num2cell(key(:,3))]);
        jobs = jobs(ord);
    end
end

function write_excluded(jobs, csvfile)
% Merge the excluded recordings into the CSV (keyed on source path).
    n = numel(jobs);
    ZT = nan(n,1); Experiment = strings(n,1); Condition = strings(n,1); CL_ms = nan(n,1);
    Tag = strings(n,1); Run = nan(n,1); File = strings(n,1); Source = strings(n,1); Rel_folder = strings(n,1);
    for j = 1:n
        m = cadence_parse_recording_name(jobs(j).source);
        ZT(j) = m.ZT; Experiment(j) = m.experiment; Condition(j) = m.condition; CL_ms(j) = m.CL_ms;
        Tag(j) = m.tag; Run(j) = m.run; File(j) = string(m.file); Source(j) = string(jobs(j).source);
        Rel_folder(j) = string(jobs(j).rel_folder);
    end
    E = table(ZT, Experiment, Condition, CL_ms, Tag, Run, File, Source, Rel_folder);
    if isfile(csvfile)
        try
            old = readtable(csvfile, 'TextType', 'string', 'Delimiter', ',');
            old = old(~ismember(string(old.Source), E.Source), :);
            E = [old(:, E.Properties.VariableNames); E];
        catch
        end
    end
    E = sortrows(E, {'ZT', 'Experiment', 'Condition', 'CL_ms'}, {'ascend', 'ascend', 'ascend', 'descend'});
    writetable(E, csvfile);
end

function ov = read_overrides(csvfile)
% {file -> rhythm} from cadence_rhythm_overrides.csv (columns File, Rhythm).
    ov = [];
    if isempty(csvfile) || ~isfile(csvfile), return; end
    R = readtable(csvfile, 'TextType', 'string', 'Delimiter', ',');
    if ~all(ismember({'File', 'Rhythm'}, R.Properties.VariableNames))
        warning('cadence_batch_extract:overrides', '%s needs columns File and Rhythm; ignored.', csvfile);
        return;
    end
    r = lower(strtrim(R.Rhythm));
    ok = ismember(r, ["paced", "arrhythmia"]);
    if any(~ok)
        warning('cadence_batch_extract:overrides', '%d row(s) of %s have a Rhythm other than paced/arrhythmia; ignored.', nnz(~ok), csvfile);
    end
    ov = struct('file', {regexprep(strtrim(R.File(ok)), '\.mat$', '')}, 'rhythm', {r(ok)});
end

function r = override_for(ov, file)
    r = '';
    if isempty(ov), return; end
    k = find(ov.file == string(regexprep(file, '\.mat$', '')), 1);
    if ~isempty(k), r = char(ov.rhythm(k)); end
end

function pdf = maps_pdf_name(metrics_file)
    pdf = regexprep(char(metrics_file), '-metrics\.mat$', '-maps.pdf');
    if strcmp(pdf, char(metrics_file)), pdf = [regexprep(pdf, '\.mat$', '') '-maps.pdf']; end
end

function rel = relative_folder(root, folder)
    root = char(root); folder = char(folder);
    if ~endsWith(root, filesep), root = [root filesep]; end
    if startsWith(folder, root), rel = folder(numel(root)+1:end); else, rel = folder; end
    rel = strrep(rel, filesep, '/');
    rel = regexprep(rel, '^/+|/+$', '');
end

function v = nz(v)
    if isempty(v) || isnan(v), v = 1e9; end
end


% =========================================================================
%  Row schema: fixed meta + every metric in cadence_metric_labels + QC.
% =========================================================================
function [names, types] = row_schema()
    meta  = {'ZT','double'; 'ZT_folder','string'; 'Experiment','string'; 'Condition','string'; ...
             'CL_ms','double'; 'Tag','string'; 'Run','double'; 'Rat','double'; 'Date','string'; ...
             'Clock','string'; 'Acquired','string'; 'Acquired_camera','string'; 'ZT_source','string'; ...
             'Stim','string'; 'File','string'; 'Source','string'; 'Metrics_file','string'};
    L = cadence_metric_labels();
    mets = [L(:,1), repmat({'double'}, size(L,1), 1)];
    tail  = {'Features','string'; 'Errors','string'; 'Elapsed_s','double'; 'Extracted_on','string'; ...
             'Maps_pdf','string'; 'Filter','string'; ...
             'Rhythm','string'; 'Capture','string'; 'Rhythm_basis','string'; 'Name_tag','double'};
    all = [meta; mets; tail];
    names = all(:,1);  types = all(:,2);
end

function s = filter_label(d)
% "50 Hz" for a single-band run; "hybrid 100/50 Hz" when rise times and
% DF/RI/OI came from the conditioned band and everything else from the low band.
    s = "";
    if ~isfield(d, 'filter_bands') || ~isstruct(d.filter_bands), return; end
    fb = d.filter_bands;
    hz = @(x) strtrim(sprintf('%g', x));
    if strcmp(fb.mode, 'hybrid')
        s = string(sprintf('hybrid %s/%s Hz', hz(fb.conditioned_hz), hz(fb.lowband_hz)));
    elseif ~isempty(fb.conditioned_hz) && fb.conditioned_hz > 0
        s = string(sprintf('%s Hz', hz(fb.conditioned_hz)));
    elseif ~isempty(fb.conditioned_hz)
        s = "none";
    end
end

function m = camera_meta(row)
% Identity columns carried into the per-camera table (Source / Metrics_file
% are the keys a re-run uses to replace its own rows).
    f = {'ZT','Experiment','Condition','CL_ms','Tag','Run','Date','Clock','Acquired','Acquired_camera','File','Source','Metrics_file'};
    m = struct();
    for k = 1:numel(f), m.(f{k}) = row.(f{k}); end
end

function Tcam = read_camera_table(csvfile)
% Per-camera table from disk, coerced to the current schema; empty when the
% file is absent or unreadable (it is rebuilt as recordings are processed).
    [names, types] = row_schema();
    meta = struct();
    for k = 1:numel(names)
        if strcmp(names{k}, 'Metrics_file'), break; end
        if strcmp(types{k}, 'double'), meta.(names{k}) = NaN; else, meta.(names{k}) = ""; end
    end
    meta.Metrics_file = "";
    meta = camera_meta(meta);
    [Tcam, mnames] = cadence_camera_medians([], 'Meta', meta);
    if isempty(csvfile) || ~isfile(csvfile), return; end
    try
        mf = fieldnames(meta);
        R = cadence_read_table(csvfile, [mf(~structfun(@isnumeric, meta)); {'Camera'; 'Signal'}]);
    catch ME
        warning('cadence_batch_extract:cameraCSV', 'Could not read %s (%s); starting it fresh.', csvfile, ME.message);
        return;
    end
    n = height(R);
    cols = Tcam.Properties.VariableNames;
    out = table();
    for k = 1:numel(cols)
        c = cols{k};
        numeric = ismember(c, mnames) || (isfield(meta, c) && isnumeric(meta.(c)));
        if ~ismember(c, R.Properties.VariableNames)
            if numeric, out.(c) = nan(n, 1); else, out.(c) = strings(n, 1); end
        elseif numeric
            v = R.(c); if ~isnumeric(v), v = str2double(string(v)); end
            out.(c) = double(v);
        else
            v = string(R.(c)); v(ismissing(v)) = ""; out.(c) = v;
        end
    end
    Tcam = out;
end

function row = new_row(job, lights_on)
    if nargin < 2, lights_on = 6; end
    [names, types] = row_schema();
    row = struct();
    for k = 1:numel(names)
        if strcmp(types{k}, 'double'), row.(names{k}) = NaN; else, row.(names{k}) = ""; end
    end
    m = cadence_parse_recording_name(job.source, lights_on);
    row.Acquired = m.acquired;  row.ZT_source = m.zt_source;
    if isfield(job, 'acquired_at') && ~isempty(job.acquired_at), row.Acquired_camera = string(job.acquired_at); end
    row.ZT = m.ZT;  row.ZT_folder = string(job.zt_folder);  row.Experiment = string(m.experiment);
    row.Condition = m.condition;  row.CL_ms = m.CL_ms;  row.Tag = m.tag;  row.Run = m.run;  row.Rat = m.rat;
    row.Date = m.date;  row.Clock = m.clock;  row.Stim = m.stim;  row.File = string(m.file);
    row.Source = string(job.source);  row.Metrics_file = string(job.metrics_file);
end

function [rows, Tcam] = refresh_meta(rows, Tcam, lights_on)
% (Re)derive ZT, ZT_source, Clock and Acquired from the file label for rows
% read back from the CSVs, so recordings extracted before these columns
% existed (or before ZT was taken from the clock) carry them too.
    for i = 1:numel(rows)
        src = char(rows{i}.Source);  if isempty(src), src = char(rows{i}.File); end
        if isempty(src), continue; end
        m = cadence_parse_recording_name(src, lights_on);
        if isnan(rows{i}.ZT) || startsWith(string(rows{i}.ZT_source), "clock") || rows{i}.ZT_source == ""
            if ~isnan(m.ZT), rows{i}.ZT = m.ZT;  rows{i}.ZT_source = m.zt_source; end
        end
        if strlength(rows{i}.Clock) == 0, rows{i}.Clock = m.clock; end
        rows{i}.Acquired = m.acquired;
    end
    if isempty(Tcam) || height(Tcam) == 0, return; end
    [usrc, ~, ic] = unique(Tcam.Source);
    for u = 1:numel(usrc)
        src = char(usrc(u));  sel = ic == u;
        if isempty(src), continue; end
        m = cadence_parse_recording_name(src, lights_on);
        if ~isnan(m.ZT), Tcam.ZT(sel & isnan(Tcam.ZT)) = m.ZT; end
        Tcam.Clock(sel & strlength(Tcam.Clock) == 0) = m.clock;
        Tcam.Acquired(sel) = m.acquired;
    end
end

function [rows, Tcam] = refresh_camera_time(rows, Tcam, jobs)
% Fill Acquired_camera from the jobs' acquired_at (the raw header's time stamp).
    if isempty(jobs) || ~isfield(jobs, 'acquired_at'), return; end
    has = ~cellfun(@isempty, {jobs.acquired_at});
    if ~any(has), return; end
    src = string({jobs(has).source});  met = string({jobs(has).metrics_file});  at = string({jobs(has).acquired_at});
    for i = 1:numel(rows)
        k = find(src == string(rows{i}.Source) | met == string(rows{i}.Metrics_file), 1);
        if ~isempty(k), rows{i}.Acquired_camera = at(k); end
    end
    if isempty(Tcam) || height(Tcam) == 0, return; end
    for k = 1:numel(at)
        sel = Tcam.Source == src(k) | Tcam.Metrics_file == met(k);
        if any(sel), Tcam.Acquired_camera(sel) = at(k); end
    end
end

function [rows, done_sources, done_metrics, Tcam] = keep_row(rows, done_sources, done_metrics, Tcam, row, Tc, job, find_row)
% Replace or append a recording's row and its per-camera rows.
    k = find_row(job);
    if isempty(k)
        rows{end+1} = row;
        done_sources{end+1} = job.source;
        done_metrics{end+1} = job.metrics_file;
    else
        rows{k} = row;
    end
    old = Tcam.Source == string(job.source) | Tcam.Metrics_file == string(job.metrics_file);
    Tcam = [Tcam(~old, :); Tc];
end

function [done, row, Tc] = adopt_metrics(job, row, filter_ok, logf)
% A metrics file exists but the medians CSV has no row for it (the CSV was
% lost, or truncated by a full disk).  Read the saved file and fill the row
% from it: seconds instead of a re-extraction.  done = false (extract it
% after all) when the file cannot be read, holds no metrics, or was made with
% another temporal filter than this run asks for.
    done = false;  Tc = table();
    t0 = tic;
    try
        S = load(job.metrics_file);
        if isfield(S, 'cmos_all_data'), d = S.cmos_all_data; else, fn = fieldnames(S); d = S.(fn{1}); end
        clear S
        if ~isfield(d, 'ep_metrics') || ~isstruct(d.ep_metrics) || isempty(fieldnames(d.ep_metrics))
            logf('   existing metrics file holds no metrics: extracting again');  return;
        end
        row.Filter = filter_label(d);
        if ~filter_ok(row.Filter)
            logf(sprintf('   existing metrics file was made with filter "%s": extracting again', row.Filter));  return;
        end
        row.Features = string(strjoin(fieldnames(d.ep_metrics)', ' '));
        if isfield(d, 'rhythm') && isstruct(d.rhythm)
            row.Rhythm = string(d.rhythm.class);  row.Capture = string(d.rhythm.capture);
            row.Rhythm_basis = string(d.rhythm.basis);
            if ~isempty(d.rhythm.note), row.Rhythm_basis = row.Rhythm_basis + " | NOTE: " + string(d.rhythm.note); end
        end
        vals = cadence_recording_medians(d);
        fn = fieldnames(vals);
        for k = 1:numel(fn), row.(fn{k}) = vals.(fn{k}); end
        row.Metrics_file = string(job.metrics_file);
        f = dir(job.metrics_file);
        row.Extracted_on = string(datestr(f.datenum, 'yyyy-mm-dd HH:MM:SS'));
        pdf = maps_pdf_name(job.metrics_file);
        if isfile(pdf), row.Maps_pdf = string(pdf); end
        Tc = cadence_camera_medians(d, 'Meta', camera_meta(row));
        done = true;
        logf(sprintf('   adopted the existing metrics file in %.0f s (medians read, no re-extraction)', toc(t0)));
    catch ME
        logf(sprintf('   existing metrics file could not be used (%s): extracting again', ME.message));
    end
end

function T = rows_to_table(rows)
    [names, types] = row_schema();
    n = numel(rows);
    T = table();
    for k = 1:numel(names)
        if strcmp(types{k}, 'double')
            col = nan(n, 1);
            for i = 1:n
                if isfield(rows{i}, names{k}) && ~isempty(rows{i}.(names{k}))
                    col(i) = double(rows{i}.(names{k}));
                end
            end
        else
            col = strings(n, 1);
            for i = 1:n
                if isfield(rows{i}, names{k}), col(i) = string(rows{i}.(names{k})); end
            end
            col(ismissing(col)) = "";
        end
        T.(names{k}) = col;
    end
end

function rows = table_to_rows(T)
    [names, types] = row_schema();
    rows = cell(1, height(T));
    for i = 1:height(T)
        r = struct();
        for k = 1:numel(names)
            if ismember(names{k}, T.Properties.VariableNames)
                v = T.(names{k})(i);
                if strcmp(types{k}, 'double')
                    if isnumeric(v), r.(names{k}) = double(v);
                    else, r.(names{k}) = str2double(string(v)); end
                else
                    v = string(v); if ismissing(v), v = ""; end
                    r.(names{k}) = v;
                end
            else
                if strcmp(types{k}, 'double'), r.(names{k}) = NaN; else, r.(names{k}) = ""; end
            end
        end
        rows{i} = r;
    end
end

function log_line(fid, s, log_fcn)
    if nargin < 3, log_fcn = []; end
    if fid > 0, line = sprintf('%s  %s', datestr(now, 'HH:MM:SS'), s); else, line = s; end
    fprintf('%s\n', line);
    if fid > 0, fprintf(fid, '%s\n', line); end
    if ~isempty(log_fcn)
        try, log_fcn(line); catch, end     % a broken console must not stop the batch
    end
end
