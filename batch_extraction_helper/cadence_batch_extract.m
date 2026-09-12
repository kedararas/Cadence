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
%     'ExcludePattern' regexp on the file name; matching recordings are NOT
%                    extracted but are listed (ZT, experiment, condition, CL,
%                    tag) in <output_root>/cadence_excluded_recordings.csv so
%                    the onset analysis still knows where arrhythmia occurred.
%                    Default '(?i)a[r]{1,2}[rh]?h?y+t?h?m' (arrhythmia
%                    recordings need the arrhythmia-dynamics stages, not the
%                    regular metrics).  The pattern is deliberately loose
%                    because the corpus contains misspellings ("Arhhythmia");
%                    it matches those without matching "alternans" or "EAD".
%                    Pass '' to extract everything.
%     'NormalizeZT'  true (default): an input folder named ZT<n>_<anything>
%                    is written as ZT<n> in the output tree (matches the
%                    lab's existing Metrics layout).  false: mirror names.
%     'Files'        cellstr of file-name patterns (regexp) to keep, default {}.
%     'Resume'       true (default): skip recordings already listed in the
%                    medians CSV whose metrics file exists.  false: redo all.
%     'DryRun'       list what would be processed and return (default false).
%     'SaveMetrics'  write the *-metrics.mat files (default true).
%     'ExtractOpts'  struct passed to cadence_extract_features (FOV_mm, ...).
%     'MediansFile'  CSV of per-recording medians
%                    (default <output_root>/cadence_recording_medians.csv)
%     'SummaryFile'  Excel workbook
%                    (default <output_root>/CADENCE_median_summary.xlsx)
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
    p.addParameter('BuildSummary', true, @islogical);
    p.addParameter('LogFile', '', @(x) ischar(x) || isstring(x));
    p.addParameter('MaxFiles', Inf, @isnumeric);
    p.addParameter('NormalizeZT', true, @islogical);
    p.addParameter('ExcludePattern', '(?i)a[r]{1,2}[rh]?h?y+t?h?m', @(x) ischar(x) || isstring(x));
    p.addParameter('Jobs', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Loader', [], @(x) isempty(x) || isa(x, 'function_handle'));
    p.addParameter('LogFcn', [], @(x) isempty(x) || isa(x, 'function_handle'));
    p.addParameter('ShouldStop', [], @(x) isempty(x) || isa(x, 'function_handle'));
    p.addParameter('SaveMapsPDF', false, @islogical);
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
    if isempty(o.LogFile),     o.LogFile     = fullfile(output_root, 'cadence_batch_log.txt'); end

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
            rows = table_to_rows(readtable(o.MediansFile, 'TextType', 'string', 'Delimiter', ','));
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
    done_sources = cellfun(@(r) char(r.Source), rows, 'UniformOutput', false);
    done_metrics = cellfun(@(r) char(r.Metrics_file), rows, 'UniformOutput', false);
    find_row = @(job) find(strcmp(done_sources, job.source) | strcmp(done_metrics, job.metrics_file), 1);

    todo = false(numel(jobs), 1);
    for j = 1:numel(jobs)
        already = ~isempty(find_row(jobs(j))) && isfile(jobs(j).metrics_file);
        todo(j) = ~(o.Resume && already);
    end
    note(sprintf('%d to process, %d skipped (already done)', nnz(todo), nnz(~todo)));

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
                    writetable(rows_to_table(rows), o.MediansFile);
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

        row = new_row(job);
        row.Extracted_on = string(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
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
            [d, st] = cadence_extract_features(d, eo);
            row.Features = string(strjoin(st.features, ' '));
            row.Errors   = string(strjoin(st.errors, ' | '));

            % save (write to a temp name, then rename, so a killed run never
            % leaves a half-written file with the final name)
            if o.SaveMetrics
                t0 = tic;
                if ~isfolder(job.out_dir), mkdir(job.out_dir); end
                tmp = [job.metrics_file '.part'];
                cmos_all_data = d; %#ok<NASGU>
                save(tmp, 'cmos_all_data', '-v7.3');
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

            % medians
            vals = cadence_recording_medians(d);
            fn = fieldnames(vals);
            for k = 1:numel(fn), row.(fn{k}) = vals.(fn{k}); end
            clear d
        catch ME
            row.Errors = string(sprintf('%s%s', char(row.Errors), [' FATAL: ' ME.message]));
            logf(sprintf('   FAILED: %s', ME.message));
        end
        row.Elapsed_s = round(toc(t_all));
        logf(sprintf('   done in %.0f s', row.Elapsed_s));

        % replace or append, then persist
        k = find_row(job);
        if isempty(k)
            rows{end+1} = row; %#ok<AGROW>
            done_sources{end+1} = job.source; %#ok<AGROW>
            done_metrics{end+1} = job.metrics_file; %#ok<AGROW>
        else
            rows{k} = row;
        end
        T = rows_to_table(rows);
        writetable(T, o.MediansFile);
        n_new = n_new + 1;
    end

    T = rows_to_table(rows);
    if ~isempty(rows), writetable(T, o.MediansFile); end
    if o.BuildSummary && ~isempty(rows)
        try
            cadence_build_summary(T, o.SummaryFile);
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
             'Clock','string'; 'Stim','string'; 'File','string'; 'Source','string'; 'Metrics_file','string'};
    L = cadence_metric_labels();
    mets = [L(:,1), repmat({'double'}, size(L,1), 1)];
    tail  = {'Features','string'; 'Errors','string'; 'Elapsed_s','double'; 'Extracted_on','string'; ...
             'Maps_pdf','string'};
    all = [meta; mets; tail];
    names = all(:,1);  types = all(:,2);
end

function row = new_row(job)
    [names, types] = row_schema();
    row = struct();
    for k = 1:numel(names)
        if strcmp(types{k}, 'double'), row.(names{k}) = NaN; else, row.(names{k}) = ""; end
    end
    m = cadence_parse_recording_name(job.source);
    row.ZT = m.ZT;  row.ZT_folder = string(job.zt_folder);  row.Experiment = string(m.experiment);
    row.Condition = m.condition;  row.CL_ms = m.CL_ms;  row.Tag = m.tag;  row.Run = m.run;  row.Rat = m.rat;
    row.Date = m.date;  row.Clock = m.clock;  row.Stim = m.stim;  row.File = string(m.file);
    row.Source = string(job.source);  row.Metrics_file = string(job.metrics_file);
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
