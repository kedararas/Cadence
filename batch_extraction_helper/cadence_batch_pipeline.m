function T = cadence_batch_pipeline(raw_root, processed_root, metrics_root, varargin)
%CADENCE_BATCH_PIPELINE  Raw camera files -> .mat -> conditioned -> metrics + summary.
%
%   T = cadence_batch_pipeline(raw_root, processed_root, metrics_root)
%   T = cadence_batch_pipeline(raw_root, processed_root, metrics_root, Name, Value, ...)
%
%   Runs the three CADENCE modules headlessly over a whole study tree:
%     1  Data Conversion      every folder under <raw_root> (any depth) that
%                             holds .tif volumes or SciMedia .gsh/.gsd pairs
%                             is a recording (cadence_raw_plan /
%                             cadence_convert_raw, the app's rules for cameras
%                             and legacy multi-recording folders)
%     2  Signal Conditioning  cadence_condition_data with the app's default
%                             settings: SVD denoising, temporal filter
%                             [0 50] Hz, drift correction, normalization,
%                             polarity check, ensemble averaging, SNR mask
%     3  Feature Extraction   cadence_batch_extract with its defaults
%                             (adaptive SNR mask per camera, ensemble beat),
%                             per-recording medians CSV and the Excel summary
%
%   Layout (the recording folder's PARENT path is mirrored, so
%   raw_root/ZT2_8am/R2/<rec>/cam files gives):
%     <processed_root>/converted/ZT2_8am/R2/<rec>.mat                raw cmos_all_data
%     <processed_root>/conditioned/ZT2_8am/R2/<rec>-conditioned.mat
%     <processed_root>/cadence_conditioning_log.csv                  one row per conditioned recording
%     <metrics_root>/ZT2/R2/<rec>-metrics.mat                        (ZT<n>_x -> ZT<n>, see NormalizeZT)
%     <metrics_root>/cadence_recording_medians.csv, CADENCE_median_summary.xlsx,
%     cadence_batch_log.txt, cadence_excluded_recordings.csv          as cadence_batch_extract writes them
%   The conditioned tree is exactly what run_cadence_batch(<processed_root>/conditioned,
%   <metrics_root>) would re-extract, so the last stage can be re-run alone.
%
%   Each recording is converted, conditioned and extracted in turn, and every
%   product is saved as it is made, so the run can be stopped and restarted:
%   a recording whose metrics are already in the medians CSV is skipped; one
%   with a conditioned file is only re-extracted; one with a converted file is
%   only re-conditioned and re-extracted.  Files are written to a .part name
%   and renamed, so an interrupted save never leaves a truncated file.
%
%   Name/Value options
%     'Folders'      cellstr of folders relative to raw_root to process; a
%                    folder selects everything below it, and may name a single
%                    recording folder.  Default {} = whole tree.
%     'Files'        cellstr of regexps on the recording (.mat base) name.
%     'SamplingHz'   frame rate used when the raw files carry none (.tif);
%                    default 1000 (the Data Conversion app's field default).
%                    .gsh headers always carry it.
%     'ConditionOpts' struct for cadence_condition_data (Drift, SVD, SVDRank,
%                    Binning, BinSize, FilterHz, Motion, Normalize, Ensemble,
%                    MaskFloor).  Default struct() = the app's defaults.
%     'ExtractOpts'  struct for cadence_extract_features (FOV_mm, ...).
%     'SaveConverted' keep the converted (raw .mat) files, default true.  They
%                    are what a re-condition with other settings starts from
%                    (a conditioned file must not be conditioned again); pass
%                    false to save disk when the raw files stay available.
%     'Resume'       true (default): skip what is already done (see above).
%                    false: re-extract every recording; existing conditioned
%                    and converted files are still reused.
%     'Recondition'  true: rebuild the conditioned files (from the converted
%                    files when present) and re-extract.  Default false.
%     'Reconvert'    true: rebuild everything from the raw files.  Default false.
%     'DryRun'       list every recording with what exists for it, then return.
%     'ExcludePattern' regexp on the conditioned file name; matching
%                    recordings (default: arrhythmia tags, see
%                    cadence_batch_extract) are converted and conditioned but
%                    NOT extracted, and listed in cadence_excluded_recordings.csv.
%                    '' = extract everything.
%     'ConditionExcluded' also convert + condition the excluded recordings
%                    (default true) so they are ready for the arrhythmia
%                    dynamics module.  false: skip them entirely.
%     'StrictConditioning' true (default): a conditioning stage that throws
%                    fails the recording (nothing saved, FATAL in the row).
%                    false: mimic the app, save and extract anyway with the
%                    stage skipped; the error stays in the conditioning log.
%     'SaveMapsPDF'  true: write <metrics_root>/<rel>/<name>-maps.pdf with the
%                    cardiac maps of every feature and camera (default false).
%     'NormalizeZT', 'MaxFiles', 'BuildSummary', 'SaveMetrics', 'MediansFile',
%     'SummaryFile', 'LogFile'      passed to cadence_batch_extract.
%     'ConditioningLog' CSV (default <processed_root>/cadence_conditioning_log.csv)
%     'LogFcn'       function handle called with every console line (a UI
%                    console); 'ShouldStop' function handle returning true to
%                    stop before the next recording (a UI STOP button).
%     'ConvertedSubdir', 'ConditionedSubdir'   default 'converted', 'conditioned'.
%
%   Returns the medians table (one row per extracted recording).
%
%   Example -- one rat first, then the whole study (resumes, skips the rat):
%     cadence_batch_paths();
%     cadence_batch_pipeline(raw, proc, met, 'Folders', {'ZT2_8am/R2'});
%     cadence_batch_pipeline(raw, proc, met);
%
%   See also run_cadence_pipeline, cadence_batch_extract, cadence_condition_data,
%   cadence_convert_raw, cadence_raw_plan.

    p = inputParser;
    p.addParameter('Folders', {}, @(x) iscellstr(x) || isstring(x) || ischar(x));
    p.addParameter('Files', {}, @(x) iscellstr(x) || isstring(x) || ischar(x));
    p.addParameter('SamplingHz', 1000, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('ConditionOpts', struct(), @isstruct);
    p.addParameter('ExtractOpts', struct(), @isstruct);
    p.addParameter('SaveConverted', true, @islogical);
    p.addParameter('Resume', true, @islogical);
    p.addParameter('Recondition', false, @islogical);
    p.addParameter('Reconvert', false, @islogical);
    p.addParameter('DryRun', false, @islogical);
    p.addParameter('ExcludePattern', '(?i)a[r]{1,2}[rh]?h?y+t?h?m', @(x) ischar(x) || isstring(x));
    p.addParameter('ConditionExcluded', true, @islogical);
    p.addParameter('StrictConditioning', true, @islogical);
    p.addParameter('NormalizeZT', true, @islogical);
    p.addParameter('MaxFiles', Inf, @isnumeric);
    p.addParameter('BuildSummary', true, @islogical);
    p.addParameter('SaveMetrics', true, @islogical);
    p.addParameter('SaveMapsPDF', false, @islogical);
    p.addParameter('MediansFile', '', @(x) ischar(x) || isstring(x));
    p.addParameter('SummaryFile', '', @(x) ischar(x) || isstring(x));
    p.addParameter('LogFile', '', @(x) ischar(x) || isstring(x));
    p.addParameter('ConditioningLog', '', @(x) ischar(x) || isstring(x));
    p.addParameter('ConvertedSubdir', 'converted', @(x) ischar(x) || isstring(x));
    p.addParameter('ConditionedSubdir', 'conditioned', @(x) ischar(x) || isstring(x));
    p.addParameter('LogFcn', [], @(x) isempty(x) || isa(x, 'function_handle'));
    p.addParameter('ShouldStop', [], @(x) isempty(x) || isa(x, 'function_handle'));
    p.parse(varargin{:});
    o = p.Results;
    o.Folders = cellstr(o.Folders);  o.Files = cellstr(o.Files);
    o.Folders = o.Folders(~cellfun(@isempty, o.Folders));
    o.Files   = o.Files(~cellfun(@isempty, o.Files));
    o.ExcludePattern = char(o.ExcludePattern);
    if o.Reconvert, o.Recondition = true; end
    if o.Recondition, o.Resume = false; end

    raw_root       = char(raw_root);
    processed_root = char(processed_root);
    metrics_root   = char(metrics_root);
    if isempty(o.LogFile),         o.LogFile         = fullfile(metrics_root, 'cadence_batch_log.txt'); end
    if isempty(o.ConditioningLog), o.ConditioningLog = fullfile(processed_root, 'cadence_conditioning_log.csv'); end
    o.LogFile = char(o.LogFile);  o.ConditioningLog = char(o.ConditioningLog);

    if ~exist('compute_lat_50', 'file'), cadence_batch_paths(); end
    if ~isfolder(raw_root), error('cadence_batch_pipeline:input', 'Raw root not found: %s', raw_root); end

    % ---- discover raw recordings ------------------------------------------------
    jobs = discover_raw(raw_root, processed_root, metrics_root, o);
    note = @(s) log_line(-1, s, o.LogFcn);
    note(sprintf('%d raw recording(s) found under %s', numel(jobs), raw_root));
    if isempty(jobs), T = table(); return; end

    excluded = false(size(jobs));
    if ~isempty(o.ExcludePattern)
        excluded = arrayfun(@(j) ~isempty(regexp(j.file, o.ExcludePattern, 'once')), jobs);
    end

    if o.DryRun
        for j = 1:numel(jobs)
            st = sprintf('[converted %s] [conditioned %s] [metrics %s]', ...
                yn(isfile(jobs(j).converted_file)), yn(isfile(jobs(j).conditioned_file)), yn(isfile(jobs(j).metrics_file)));
            tag = '';  if excluded(j), tag = '  (excluded from extraction)'; end
            note(sprintf('  %-24s %s  %s, %d cam(s)  %s%s', jobs(j).rel_folder, jobs(j).mat_base, ...
                jobs(j).raw_entry.ext, numel(jobs(j).raw_entry.cam_names), st, tag));
        end
        T = table();
        return;
    end

    if ~isfolder(processed_root), mkdir(processed_root); end
    if ~isfolder(metrics_root),   mkdir(metrics_root);   end
    loader = @(job, logf) pipeline_loader(job, logf, o);

    % ---- excluded recordings: convert + condition only ----------------------------
    fid = fopen(o.LogFile, 'a');
    logf = @(s) log_line(fid, s, o.LogFcn);
    logf(sprintf('==== cadence_batch_pipeline started %s ====', datestr(now, 'yyyy-mm-dd HH:MM:SS')));
    logf(sprintf('raw       : %s', raw_root));
    logf(sprintf('processed : %s', processed_root));
    logf(sprintf('metrics   : %s', metrics_root));
    if any(excluded) && o.ConditionExcluded
        idx = find(excluded);
        logf(sprintf('%d excluded recording(s): converting + conditioning only', numel(idx)));
        for jj = 1:numel(idx)
            if ~isempty(o.ShouldStop) && o.ShouldStop()
                logf('STOP requested: stopping before the next recording (done so far is kept).');
                break;
            end
            job = jobs(idx(jj));
            logf(sprintf('[excluded %d/%d] %s', jj, numel(idx), job.source));
            try
                loader(job, logf);
            catch ME
                logf(sprintf('   FAILED: %s', ME.message));
            end
        end
    end
    fclose(fid);

    % ---- everything else: convert + condition on demand, then extract -------------
    T = cadence_batch_extract(raw_root, metrics_root, 'Jobs', jobs, 'Loader', loader, ...
        'Resume', o.Resume, 'ExcludePattern', o.ExcludePattern, 'ExtractOpts', o.ExtractOpts, ...
        'NormalizeZT', o.NormalizeZT, 'MaxFiles', o.MaxFiles, 'BuildSummary', o.BuildSummary, ...
        'SaveMetrics', o.SaveMetrics, 'MediansFile', o.MediansFile, 'SummaryFile', o.SummaryFile, ...
        'LogFile', o.LogFile, 'LogFcn', o.LogFcn, 'ShouldStop', o.ShouldStop, 'SaveMapsPDF', o.SaveMapsPDF);
end


% =========================================================================
function d = pipeline_loader(job, logf, o)
% Return the conditioned cmos_all_data for a job: load it if it exists,
% otherwise condition the converted .mat if it exists, otherwise convert the
% raw folder first.  Every product is saved as it is made.
    name = job.file;
    if ~o.Recondition && isfile(job.conditioned_file)
        t0 = tic;
        d = load_struct(job.conditioned_file);
        logf(sprintf('   conditioned file exists, loaded in %.0f s: %s', toc(t0), job.conditioned_file));
        return;
    end

    t_all = tic;
    row = cond_row(job);
    try
        % ---- raw data: converted file or the raw camera files ----------------
        if ~o.Reconvert && isfile(job.converted_file)
            t0 = tic;
            raw = load_struct(job.converted_file);
            logf(sprintf('   converted file exists, loaded in %.0f s: %s', toc(t0), job.converted_file));
            row.Converted = "reused";
        else
            t0 = tic;
            raw = cadence_convert_raw(job.raw_entry, o.SamplingHz, @(s) logf(['   ' s]));
            logf(sprintf('   converted %d camera file(s) from %s in %.0f s', ...
                numel(job.raw_entry.cam_names), job.raw_folder, toc(t0)));
            [ok, msgs] = schema_gate(raw, "raw", name);
            for k = 1:numel(msgs), logf(['   ' msgs{k}]); end
            if ~ok, error('cadence_batch_pipeline:rawSchema', 'raw schema FAIL after conversion'); end
            if o.SaveConverted
                t0 = tic;
                save_atomic(job.converted_file, raw);
                logf(sprintf('   saved %s in %.0f s', job.converted_file, toc(t0)));
                row.Converted = "saved";
            else
                row.Converted = "not saved";
            end
        end
        row.Cameras    = double(raw.num_files);
        row.Frames     = size(raw.CAM1, 3);
        row.acqFreq_Hz = double(raw.acqFreq);
        row.Pacing     = double(isfield(raw, 'analog1') && ~isempty(raw.analog1));

        % ---- condition ---------------------------------------------------------
        co = o.ConditionOpts;  co.Name = name;  co.Log = @(s) logf(['   ' s]);
        t0 = tic;
        [d, st] = cadence_condition_data(raw, co);
        clear raw
        logf(sprintf('   conditioned in %.0f s', toc(t0)));
        row.Polarity = string(strjoin(arrayfun(@(q) sprintf('CAM%d %s (%.2f, %s)', q.cam, ...
            tern(q.inverted, 'inverted', 'ok'), q.confidence, q.method), st.polarity, 'UniformOutput', false), '; '));
        if ~isempty(st.qc)
            nonpass = st.qc(arrayfun(@(x) x.status ~= "PASS", st.qc));
            row.QC_fail = nnz(arrayfun(@(x) x.status == "FAIL", st.qc));
            row.QC_warn = nnz(arrayfun(@(x) x.status == "WARN", st.qc));
            row.QC = string(strjoin(arrayfun(@(x) sprintf('%s %s: %s', x.status, x.check, x.message), ...
                nonpass, 'UniformOutput', false), ' | '));
        end
        row.Stage_errors = string(strjoin(st.errors, ' | '));
        if ~isempty(st.errors) && o.StrictConditioning
            error('cadence_batch_pipeline:conditioning', 'conditioning stage error(s): %s', strjoin(st.errors, ' | '));
        end

        t0 = tic;
        save_atomic(job.conditioned_file, d);
        logf(sprintf('   saved %s in %.0f s', job.conditioned_file, toc(t0)));
        row.Saved = 1;
    catch ME
        row.Errors = string(ME.message);
        row.Elapsed_s = round(toc(t_all));
        row.Conditioned_on = string(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        append_cond_log(o.ConditioningLog, row);
        rethrow(ME);
    end
    row.Elapsed_s = round(toc(t_all));
    row.Conditioned_on = string(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    append_cond_log(o.ConditioningLog, row);
end


% =========================================================================
function jobs = discover_raw(raw_root, processed_root, metrics_root, o)
% Every folder under raw_root holding .tif or .gsh files is a recording
% folder; cadence_raw_plan splits legacy multi-recording folders.  The
% recording folder's parent path is mirrored under the output roots.
    jobs = struct('source', {}, 'file', {}, 'zt_folder', {}, 'rat_folder', {}, 'zt_out', {}, ...
                  'out_dir', {}, 'metrics_file', {}, 'rel_folder', {}, ...
                  'raw_folder', {}, 'raw_entry', {}, 'mat_base', {}, 'converted_file', {}, 'conditioned_file', {});
    folders = strrep(o.Folders, '\', '/');
    folders = regexprep(folders, '/+$', '');

    listing = [dir(fullfile(raw_root, '**', '*.gsh')); dir(fullfile(raw_root, '**', '*.tif'))];
    listing = listing(~[listing.isdir]);
    listing = listing(~startsWith({listing.name}, '.'));
    rec_folders = unique({listing.folder});
    for c = 1:numel(rec_folders)
        folder = rec_folders{c};
        rel = relative_folder(raw_root, folder);
        if any(cellfun(@(q) startsWith(q, '.'), strsplit(rel, '/'))), continue; end   % hidden dirs
        if ~isempty(folders)
            hit = false;
            for f = 1:numel(folders)
                if strcmp(rel, folders{f}) || startsWith(rel, [folders{f} '/']), hit = true; break; end
            end
            if ~hit, continue; end
        end
        rel_parent = regexprep(rel, '(^|/)[^/]*$', '');
        out_rel = rel_parent;
        if o.NormalizeZT
            out_rel = regexprep(out_rel, '(^|/)(ZT\d+)[^/]*', '$1$2');
        end

        plan = cadence_raw_plan(folder);
        for e = 1:numel(plan)
            base = plan(e).mat_base;
            if ~isempty(o.Files) && ~any(cellfun(@(pt) ~isempty(regexp(base, pt, 'once')), o.Files))
                continue;
            end
            j = struct();
            j.file             = [base '-conditioned.mat'];
            j.conditioned_file = fullfile(processed_root, char(o.ConditionedSubdir), rel_parent, j.file);
            j.source           = j.conditioned_file;
            j.converted_file   = fullfile(processed_root, char(o.ConvertedSubdir), rel_parent, [base '.mat']);
            j.raw_folder       = folder;
            j.raw_entry        = plan(e);
            j.mat_base         = base;
            m = cadence_parse_recording_name(j.source);
            j.zt_folder  = char(m.zt_folder);
            j.rat_folder = char(m.rat_folder);
            if isnan(m.ZT), j.zt_out = ''; else, j.zt_out = sprintf('ZT%d', m.ZT); end
            j.out_dir      = fullfile(metrics_root, out_rel);
            j.metrics_file = fullfile(j.out_dir, [base '-metrics.mat']);
            j.rel_folder   = rel_parent;
            jobs(end+1) = orderfields(j, jobs); %#ok<AGROW>
        end
    end
    % order: ZT number, experiment number, relative folder, run number (as cadence_batch_extract)
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

function rel = relative_folder(root, folder)
    root = char(root); folder = char(folder);
    if ~endsWith(root, filesep), root = [root filesep]; end
    if startsWith(folder, root), rel = folder(numel(root)+1:end); else, rel = ''; end
    rel = strrep(rel, filesep, '/');
    rel = regexprep(rel, '^/+|/+$', '');
end

function v = nz(v)
    if isempty(v) || isnan(v), v = 1e9; end
end

function s = yn(b)
    if b, s = 'yes'; else, s = 'no '; end
end

function v = tern(c, a, b)
    if c, v = a; else, v = b; end
end


% =========================================================================
function d = load_struct(file)
    S = load(file);
    if isfield(S, 'cmos_all_data'), d = S.cmos_all_data;
    else, fn = fieldnames(S); d = S.(fn{1}); end
end

function save_atomic(file, d)
% Write to a .part name, then rename, so a killed run never leaves a
% half-written file with the final name.
    out_dir = fileparts(file);
    if ~isfolder(out_dir), mkdir(out_dir); end
    tmp = [file '.part'];
    cmos_all_data = d; %#ok<NASGU>
    save(tmp, 'cmos_all_data', '-v7.3');
    if isfile(file), delete(file); end
    movefile(tmp, file);
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


% =========================================================================
%  Conditioning log: one row per conditioned recording, keyed on the
%  conditioned file, merged into the CSV after every recording.
% =========================================================================
function [names, types] = cond_log_schema()
    s = {'Conditioned_file','string'; 'Raw_folder','string'; 'Converted_file','string'; ...
         'Converted','string'; 'Format','string'; 'Cameras','double'; 'Frames','double'; ...
         'acqFreq_Hz','double'; 'Pacing','double'; 'Polarity','string'; 'QC_fail','double'; ...
         'QC_warn','double'; 'QC','string'; 'Stage_errors','string'; 'Saved','double'; ...
         'Errors','string'; 'Elapsed_s','double'; 'Conditioned_on','string'};
    names = s(:,1);  types = s(:,2);
end

function row = cond_row(job)
    [names, types] = cond_log_schema();
    row = struct();
    for k = 1:numel(names)
        if strcmp(types{k}, 'double'), row.(names{k}) = NaN; else, row.(names{k}) = ""; end
    end
    row.Conditioned_file = string(job.conditioned_file);
    row.Raw_folder       = string(job.raw_folder);
    row.Converted_file   = string(job.converted_file);
    row.Format           = string(job.raw_entry.ext);
    row.Saved            = 0;
end

function append_cond_log(csvfile, row)
    [names, types] = cond_log_schema();
    rows = {};
    if isfile(csvfile)
        try
            old = readtable(csvfile, 'TextType', 'string', 'Delimiter', ',');
            for i = 1:height(old)
                r = struct();
                for k = 1:numel(names)
                    if ismember(names{k}, old.Properties.VariableNames)
                        v = old.(names{k})(i);
                        if strcmp(types{k}, 'double')
                            if isnumeric(v), r.(names{k}) = double(v); else, r.(names{k}) = str2double(string(v)); end
                        else
                            v = string(v); if ismissing(v), v = ""; end
                            r.(names{k}) = v;
                        end
                    elseif strcmp(types{k}, 'double'), r.(names{k}) = NaN;
                    else, r.(names{k}) = "";
                    end
                end
                rows{end+1} = r; %#ok<AGROW>
            end
        catch ME
            warning('cadence_batch_pipeline:condlog', 'Could not read %s (%s); rewriting it.', csvfile, ME.message);
            rows = {};
        end
    end
    keys = cellfun(@(r) char(r.Conditioned_file), rows, 'UniformOutput', false);
    k = find(strcmp(keys, char(row.Conditioned_file)), 1);
    if isempty(k), rows{end+1} = row; else, rows{k} = row; end

    n = numel(rows);
    T = table();
    for k = 1:numel(names)
        if strcmp(types{k}, 'double')
            col = nan(n, 1);
            for i = 1:n, if ~isempty(rows{i}.(names{k})), col(i) = double(rows{i}.(names{k})); end, end
        else
            col = strings(n, 1);
            for i = 1:n, col(i) = string(rows{i}.(names{k})); end
            col(ismissing(col)) = "";
        end
        T.(names{k}) = col;
    end
    out_dir = fileparts(csvfile);
    if ~isempty(out_dir) && ~isfolder(out_dir), mkdir(out_dir); end
    writetable(T, csvfile);
end
