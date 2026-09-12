function T = cadence_recompute_medians(output_root, varargin)
%CADENCE_RECOMPUTE_MEDIANS  Rebuild the medians CSV and workbook from saved metrics files.
%
%   T = cadence_recompute_medians(output_root)
%   T = cadence_recompute_medians(output_root, 'Folders', {'ZT2/R2'}, 'BuildSummary', true)
%
%   Scans <output_root>/**/*-metrics.mat, recomputes cadence_recording_medians
%   for each (no feature extraction, ~10-15 s per file for loading) and
%   rewrites <output_root>/cadence_recording_medians.csv and the Excel summary.
%   Use it after changing cadence_metric_labels / cadence_recording_medians,
%   or to summarise metrics files produced by the app itself.
%
%   Options: 'Folders' (relative folders under output_root), 'MediansFile',
%   'SummaryFile', 'BuildSummary' (default true), 'SaveMapsPDF' (default
%   false: write <name>-maps.pdf next to every metrics file, the cardiac maps
%   of each feature and camera, see cadence_maps_pdf; existing PDFs are
%   rewritten), 'ExcludePattern'.

    p = inputParser;
    p.addParameter('Folders', {}, @(x) iscellstr(x) || isstring(x) || ischar(x));
    p.addParameter('MediansFile', '', @(x) ischar(x) || isstring(x));
    p.addParameter('SummaryFile', '', @(x) ischar(x) || isstring(x));
    p.addParameter('BuildSummary', true, @islogical);
    p.addParameter('SaveMapsPDF', false, @islogical);
    p.addParameter('ExcludePattern', '(?i)a[r]{1,2}[rh]?h?y+t?h?m', @(x) ischar(x) || isstring(x));
    p.parse(varargin{:});
    o = p.Results;
    o.Folders = cellstr(o.Folders); o.Folders = o.Folders(~cellfun(@isempty, o.Folders));
    output_root = char(output_root);
    if isempty(o.MediansFile), o.MediansFile = fullfile(output_root, 'cadence_recording_medians.csv'); end
    if isempty(o.SummaryFile), o.SummaryFile = fullfile(output_root, 'CADENCE_median_summary.xlsx'); end
    if ~exist('compute_lat_50', 'file'), cadence_batch_paths(); end

    files = dir(fullfile(output_root, '**', '*-metrics.mat'));
    files = files(~[files.isdir] & ~startsWith({files.name}, '.'));
    keep = true(numel(files), 1);
    for k = 1:numel(files)
        rel = strrep(files(k).folder, [output_root filesep], ''); rel = strrep(rel, filesep, '/');
        if any(cellfun(@(q) startsWith(q, '.'), strsplit(rel, '/'))), keep(k) = false; end
        if ~isempty(o.Folders)
            keep(k) = keep(k) && any(cellfun(@(f) strcmp(rel, f) || startsWith(rel, [f '/']), strrep(o.Folders, '\', '/')));
        end
    end
    if ~isempty(char(o.ExcludePattern))
        keep = keep & cellfun(@(nm) isempty(regexp(nm, char(o.ExcludePattern), 'once')), {files.name}');
    end
    files = files(keep);
    fprintf('%d metrics file(s) under %s\n', numel(files), output_root);

    L = cadence_metric_labels();
    rows = cell(1, numel(files));
    for k = 1:numel(files)
        f = fullfile(files(k).folder, files(k).name);
        t0 = tic;
        m = cadence_parse_recording_name(f);
        r = struct('ZT', m.ZT, 'ZT_folder', m.zt_folder, 'Experiment', m.experiment, 'Condition', m.condition, ...
                   'CL_ms', m.CL_ms, 'Tag', m.tag, 'Run', m.run, 'Rat', m.rat, 'Date', m.date, 'Clock', m.clock, ...
                   'Stim', m.stim, 'File', string([m.base '.mat']), 'Source', "", 'Metrics_file', string(f));
        for q = 1:size(L,1), r.(L{q,1}) = NaN; end
        r.Features = ""; r.Errors = ""; r.Elapsed_s = NaN; r.Extracted_on = ""; r.Maps_pdf = "";
        try
            S = load(f);
            if isfield(S, 'cmos_all_data'), d = S.cmos_all_data; else, fn = fieldnames(S); d = S.(fn{1}); end
            clear S
            if isfield(d, 'ep_metrics'), r.Features = string(strjoin(fieldnames(d.ep_metrics)', ' ')); end
            v = cadence_recording_medians(d);
            fn = fieldnames(v);
            for q = 1:numel(fn), r.(fn{q}) = v.(fn{q}); end
            r.Extracted_on = string(datestr(files(k).datenum, 'yyyy-mm-dd HH:MM:SS'));
            pdf = regexprep(f, '-metrics\.mat$', '-maps.pdf');
            if o.SaveMapsPDF
                try
                    n = cadence_maps_pdf(d, pdf, 'Name', files(k).name);
                    if n > 0, r.Maps_pdf = string(pdf); end
                catch ME
                    r.Errors = string(['maps PDF: ' ME.message]);
                end
            elseif isfile(pdf)
                r.Maps_pdf = string(pdf);
            end
            clear d
        catch ME
            r.Errors = string(['FATAL: ' ME.message]);
        end
        r.Elapsed_s = round(toc(t0));
        rows{k} = r;
        fprintf('[%d/%d] %-70s %3.0f s\n', k, numel(files), files(k).name, r.Elapsed_s);
    end

    T = struct2table_uniform(rows);
    writetable(T, o.MediansFile);
    fprintf('wrote %s\n', o.MediansFile);
    if o.BuildSummary && ~isempty(rows)
        cadence_build_summary(T, o.SummaryFile);
    end
end

function T = struct2table_uniform(rows)
    if isempty(rows), T = table(); return; end
    names = fieldnames(rows{1});
    T = table();
    for k = 1:numel(names)
        v0 = rows{1}.(names{k});
        if isnumeric(v0)
            col = nan(numel(rows), 1);
            for i = 1:numel(rows), x = rows{i}.(names{k}); if ~isempty(x), col(i) = double(x); end; end
        else
            col = strings(numel(rows), 1);
            for i = 1:numel(rows), col(i) = string(rows{i}.(names{k})); end
            col(ismissing(col)) = "";
        end
        T.(names{k}) = col;
    end
end
