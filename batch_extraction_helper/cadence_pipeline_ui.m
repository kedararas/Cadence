function fig = cadence_pipeline_ui()
%CADENCE_PIPELINE_UI  Window for run_cadence_pipeline: pick three folders and run.
%
%   cadence_pipeline_ui
%   fig = cadence_pipeline_ui
%
%   A small front end for the headless raw -> converted -> conditioned ->
%   metrics pipeline (cadence_batch_pipeline), laid out like the other
%   CADENCE modules:
%     STEP 1  raw data parent folder (SELECT / REFRESH), a list of its
%             sub-folders to restrict the run (select none = whole tree), and
%             the sampling frequency used for .tif recordings
%     STEP 2  processed folder (gets converted/ and conditioned/ sub-trees and
%             the conditioning log)
%     STEP 3  metrics output folder (metrics files, medians CSV, workbook, log)
%     STEP 4  the Signal Conditioning dropdowns with their app defaults, the
%             extraction field of view, and the run options (dry run, resume,
%             re-condition, re-convert, keep converted files, condition the
%             arrhythmia-tagged recordings, save a PDF of the cardiac maps)
%   RUN PIPELINE runs cadence_batch_pipeline with those settings; every
%   console line of the batch is echoed in the text area (newest on top, as
%   in the other modules) and STOP ends the run before the next recording,
%   keeping everything done so far.  The same messages also go to
%   <metrics>/cadence_batch_log.txt.
%
%   Console-free equivalent: run_cadence_pipeline(raw, processed, metrics, ...).
%
%   See also run_cadence_pipeline, cadence_batch_pipeline, cadence_condition_data.

    addpath(fileparts(mfilename('fullpath')));
    cadence_batch_paths();

    black = [0 0 0]; white = [1 1 1]; blue = [0 0.4471 0.7412];

    fig = uifigure('Name', 'Cadence - Batch Pipeline (raw -> converted -> conditioned -> metrics)', ...
                   'Position', [300 80 1240 940], 'Visible', 'off');
    g = uigridlayout(fig, [6 1]);
    g.RowHeight = {250, 78, 78, 228, 48, '1x'};
    g.Padding = [8 8 8 8];  g.RowSpacing = 8;
    h = struct();

    % ---- STEP 1: raw parent folder ----------------------------------------------
    p1 = uipanel(g, 'Title', 'STEP 1: SELECT RAW DATA PARENT FOLDER (sub-folders holding .tif or .gsh/.gsd camera files)', ...
                 'FontWeight', 'bold');
    p1.Layout.Row = 1;
    g1 = uigridlayout(p1, [3 3]);
    g1.RowHeight = {26, '1x', 26};  g1.ColumnWidth = {'1x', 190, 190};  g1.Padding = [6 6 6 6];
    h.RawField = uieditfield(g1, 'text', 'Placeholder', 'raw data parent folder', ...
                             'ValueChangedFcn', @(s, e) refresh_list(fig));
    h.RawField.Layout.Row = 1;  h.RawField.Layout.Column = 1;
    h.RawButton = uibutton(g1, 'push', 'Text', 'SELECT RAW FOLDER', 'BackgroundColor', black, 'FontColor', white, ...
                           'ButtonPushedFcn', @(s, e) select_dir(fig, 'RawField', true));
    h.RawButton.Layout.Row = 1;  h.RawButton.Layout.Column = 2;
    h.RefreshButton = uibutton(g1, 'push', 'Text', 'REFRESH FOLDER LIST', 'BackgroundColor', black, 'FontColor', white, ...
                               'ButtonPushedFcn', @(s, e) refresh_list(fig));
    h.RefreshButton.Layout.Row = 1;  h.RefreshButton.Layout.Column = 3;
    h.FolderList = uilistbox(g1, 'Items', {}, 'Value', {}, 'Multiselect', 'on');
    h.FolderList.Layout.Row = 2;  h.FolderList.Layout.Column = [1 3];
    lbl = uilabel(g1, 'Text', ['Sub-folders to process (each with everything below it). ' ...
                  'Select none to process the whole tree.']);
    lbl.Layout.Row = 3;  lbl.Layout.Column = 1;
    lbl = uilabel(g1, 'Text', 'Sampling freq. (Hz), .tif only', 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
    lbl.Layout.Row = 3;  lbl.Layout.Column = 2;
    h.FsField = uieditfield(g1, 'numeric', 'Value', 1000, 'Limits', [1 Inf], 'FontWeight', 'bold');
    h.FsField.Layout.Row = 3;  h.FsField.Layout.Column = 3;

    % ---- STEP 2 / 3: processed and metrics folders ------------------------------
    [h.ProcField, h.ProcButton] = dir_row(g, 2, ...
        'STEP 2: SELECT PROCESSED FOLDER (converted/ and conditioned/ trees are written here)', ...
        'SELECT PROCESSED FOLDER', 'ProcField', 'processed folder', ...
        'NOTE: if not selected, a folder named <raw folder>_processed next to the raw folder is used.');
    [h.MetField, h.MetButton] = dir_row(g, 3, ...
        'STEP 3: SELECT METRICS OUTPUT FOLDER (metrics files, medians CSV, Excel summary, log)', ...
        'SELECT METRICS FOLDER', 'MetField', 'metrics output folder', ...
        'NOTE: if not selected, a folder named <raw folder>_metrics next to the raw folder is used.');

    % ---- STEP 4: settings -------------------------------------------------------
    p4 = uipanel(g, 'Title', 'STEP 4: SIGNAL CONDITIONING SETTINGS (defaults = recommended), EXTRACTION AND RUN OPTIONS', ...
                 'FontWeight', 'bold');
    p4.Layout.Row = 4;
    g4 = uigridlayout(p4, [5 8]);
    g4.RowHeight = {26, 26, 26, 26, 26};  g4.ColumnWidth = {'1x', 100, '1x', 100, '1x', 100, '1x', 100};
    g4.Padding = [6 6 6 6];
    h.SVD      = dd(g4, 1, 1, 'A) SVD DENOISING', {'YES', 'NO'}, 'YES');
    h.Binning  = dd(g4, 1, 3, 'B) BINNING', {'NO', 'YES'}, 'NO');
    h.BinSize  = dd(g4, 1, 5, 'B) BINNING BOX', {'3 x 3', '5 x 5', '7 x 7', '9 x 9'}, '3 x 3');
    h.Filter   = dd(g4, 1, 7, 'C) TEMPORAL FILTER (Hz)', {'NONE', '[0, 150]', '[0, 100]', '[0, 75]', '[0, 50]'}, '[0, 50]');
    h.Drift    = dd(g4, 2, 1, 'D) DRIFT CORRECTION', {'NO', 'YES'}, 'YES');
    h.Normalize= dd(g4, 2, 3, 'E) NORMALIZATION', {'NO', 'YES'}, 'YES');
    h.Ensemble = dd(g4, 2, 5, 'F) ENSEMBLE AVERAGING', {'NO', 'YES'}, 'YES');
    h.Motion   = dd(g4, 2, 7, 'G) MOTION CORRECTION', {'NO', 'YES'}, 'NO');
    h.DryRun      = cb(g4, 3, [1 2], 'Dry run (list only, nothing processed)', false);
    h.Resume      = cb(g4, 3, [3 4], 'Resume (skip recordings already done)', true);
    h.Recondition = cb(g4, 3, [5 6], 'Re-condition (keep converted files)', false);
    h.Reconvert   = cb(g4, 3, [7 8], 'Re-convert (start over from raw files)', false);
    h.SaveConv    = cb(g4, 4, [1 2], 'Keep converted (raw .mat) files', true);
    h.CondExcl    = cb(g4, 4, [3 4], 'Condition arrhythmia recordings too (not extracted)', true);
    lbl = uilabel(g4, 'Text', 'Field of view for CV (mm)', 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
    lbl.Layout.Row = 4;  lbl.Layout.Column = [5 7];
    h.FOV = uieditfield(g4, 'numeric', 'Value', 20, 'Limits', [0.1 Inf], 'FontWeight', 'bold');
    h.FOV.Layout.Row = 4;  h.FOV.Layout.Column = 8;
    h.MapsPDF     = cb(g4, 5, [1 4], 'Save a PDF of the cardiac maps for each recording (<name>-maps.pdf next to the metrics file)', false);

    % ---- run / stop --------------------------------------------------------------
    g5 = uigridlayout(g, [1 2]);
    g5.Layout.Row = 5;  g5.ColumnWidth = {'1x', 180};  g5.Padding = [0 0 0 0];
    h.RunButton = uibutton(g5, 'push', 'Text', 'RUN PIPELINE: CONVERT -> CONDITION -> EXTRACT FEATURES -> SUMMARY', ...
                           'BackgroundColor', blue, 'FontColor', white, 'FontSize', 16, 'FontWeight', 'bold', ...
                           'ButtonPushedFcn', @(s, e) run_pressed(fig));
    h.RunButton.Layout.Column = 1;
    h.StopButton = uibutton(g5, 'push', 'Text', 'STOP AFTER CURRENT FILE', 'BackgroundColor', [0.65 0.1 0.1], ...
                            'FontColor', white, 'FontWeight', 'bold', 'Enable', 'off', ...
                            'ButtonPushedFcn', @(s, e) stop_pressed(fig));
    h.StopButton.Layout.Column = 2;

    % ---- console -----------------------------------------------------------------
    h.Console = uitextarea(g, 'Editable', 'off', 'FontWeight', 'bold', 'FontColor', [0 0.451 0.7412]);
    h.Console.Layout.Row = 6;
    h.Console.Value = {'Welcome to CADENCE 2026!! '; ''; ...
        'STEP 1) Select the parent folder of the raw recordings. Every sub-folder (any depth) that holds .tif volumes or SciMedia .gsh/.gsd pairs is one recording. Optionally pick sub-folders to restrict the run.'; ''; ...
        'STEP 2) Select the processed folder: converted .mat files go to converted/, conditioned files to conditioned/, with a conditioning log.'; ''; ...
        'STEP 3) Select the metrics output folder: -metrics.mat files, the per-recording medians CSV, the Excel summary and the batch log.'; ''; ...
        'STEP 4) Check the conditioning settings (defaults = recommended: SVD denoising, [0, 50] Hz filter, drift correction, normalization, ensemble averaging) and run. Features are extracted with the adaptive SNR mask on the ensemble beat. The run resumes where it stopped.'};

    fig.UserData = struct('h', h, 'stop', false, 'running', false, 'last_table', table());
    fig.Visible = 'on';
    if nargout == 0, clear fig; end
end


% =========================================================================
%  Layout helpers
% =========================================================================
function [field, button] = dir_row(g, row, title, button_text, field_name, placeholder, note)
    p = uipanel(g, 'Title', title, 'FontWeight', 'bold');
    p.Layout.Row = row;
    gg = uigridlayout(p, [2 2]);
    gg.RowHeight = {26, 18};  gg.ColumnWidth = {'1x', 190};  gg.Padding = [6 4 6 4];  gg.RowSpacing = 2;
    fig = ancestor(g, 'figure');
    field = uieditfield(gg, 'text', 'Placeholder', placeholder);
    field.Layout.Row = 1;  field.Layout.Column = 1;
    button = uibutton(gg, 'push', 'Text', button_text, 'BackgroundColor', [0 0 0], 'FontColor', [1 1 1], ...
                      'ButtonPushedFcn', @(s, e) select_dir(fig, field_name, false));
    button.Layout.Row = 1;  button.Layout.Column = 2;
    lbl = uilabel(gg, 'Text', note);
    lbl.Layout.Row = 2;  lbl.Layout.Column = [1 2];
end

function d = dd(g, row, col, label, items, value)
    lbl = uilabel(g, 'Text', label, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    lbl.Layout.Row = row;  lbl.Layout.Column = col;
    d = uidropdown(g, 'Items', items, 'Value', value, 'FontWeight', 'bold');
    d.Layout.Row = row;  d.Layout.Column = col + 1;
end

function c = cb(g, row, cols, text, value)
    c = uicheckbox(g, 'Text', text, 'Value', value);
    c.Layout.Row = row;  c.Layout.Column = cols;
end


% =========================================================================
%  Callbacks
% =========================================================================
function select_dir(fig, field_name, is_raw)
    h = fig.UserData.h;
    start = h.(field_name).Value;
    if isempty(start) || ~isfolder(start), start = pwd; end
    try
        fig.Visible = 'off';                     % as the other modules do (modal dialog on macOS)
        chosen = uigetdir(start);
        fig.Visible = 'on';
    catch ME
        fig.Visible = 'on';
        console(fig, ['ERROR - select folder: ' ME.message]);
        return;
    end
    if isequal(chosen, 0)
        console(fig, 'Folder not updated.');
        return;
    end
    h.(field_name).Value = chosen;
    console(fig, sprintf('%s selected: %s', strrep(field_name, 'Field', ' folder'), chosen));
    if is_raw, refresh_list(fig); end
end

function refresh_list(fig)
    h = fig.UserData.h;
    raw = h.RawField.Value;
    h.FolderList.Items = {};  h.FolderList.Value = {};
    if isempty(raw) || ~isfolder(raw)
        if ~isempty(raw), console(fig, sprintf('Raw folder not found: %s', raw)); end
        return;
    end
    d = dir(raw);
    d = d([d.isdir]);
    names = {d.name};
    names = names(~startsWith(names, '.'));
    h.FolderList.Items = sort(names);
    console(fig, sprintf('%d sub-folder(s) found in %s', numel(names), raw));
end

function stop_pressed(fig)
    fig.UserData.stop = true;
    fig.UserData.h.StopButton.Enable = 'off';
    console(fig, 'STOP requested: the run ends after the current recording.');
end

function run_pressed(fig)
    if fig.UserData.running, return; end
    h = fig.UserData.h;

    raw = strtrim(h.RawField.Value);
    if isempty(raw) || ~isfolder(raw)
        console(fig, 'Please select an existing raw data parent folder (STEP 1).');
        return;
    end
    raw = regexprep(raw, '[\\/]+$', '');
    [parent, raw_name] = fileparts(raw);
    proc = strtrim(h.ProcField.Value);
    if isempty(proc)
        proc = fullfile(parent, [raw_name '_processed']);
        h.ProcField.Value = proc;
        console(fig, sprintf('Processed folder not selected; using %s', proc));
    end
    met = strtrim(h.MetField.Value);
    if isempty(met)
        met = fullfile(parent, [raw_name '_metrics']);
        h.MetField.Value = met;
        console(fig, sprintf('Metrics folder not selected; using %s', met));
    end

    % conditioning settings -> cadence_condition_data options
    co = struct();
    co.SVD       = strcmp(h.SVD.Value, 'YES');
    co.Binning   = strcmp(h.Binning.Value, 'YES');
    co.BinSize   = str2double(regexp(h.BinSize.Value, '^\d+', 'match', 'once'));
    co.Drift     = strcmp(h.Drift.Value, 'YES');
    co.Normalize = strcmp(h.Normalize.Value, 'YES');
    co.Ensemble  = strcmp(h.Ensemble.Value, 'YES');
    co.Motion    = strcmp(h.Motion.Value, 'YES');
    fz = regexp(h.Filter.Value, '\[0,\s*(\d+)\]', 'tokens', 'once');
    if isempty(fz), co.FilterHz = []; else, co.FilterHz = str2double(fz{1}); end
    eo = struct('FOV_mm', h.FOV.Value);

    folders = cellstr(h.FolderList.Value);
    folders = folders(~cellfun(@isempty, folders));

    args = {'Folders', folders, 'SamplingHz', h.FsField.Value, 'ConditionOpts', co, 'ExtractOpts', eo, ...
            'DryRun', logical(h.DryRun.Value), 'Resume', logical(h.Resume.Value), ...
            'Recondition', logical(h.Recondition.Value), 'Reconvert', logical(h.Reconvert.Value), ...
            'SaveConverted', logical(h.SaveConv.Value), 'ConditionExcluded', logical(h.CondExcl.Value), ...
            'SaveMapsPDF', logical(h.MapsPDF.Value), ...
            'LogFcn', @(s) console(fig, s), 'ShouldStop', @() should_stop(fig)};

    fig.UserData.stop = false;
    fig.UserData.running = true;
    h.RunButton.Enable = 'off';
    h.StopButton.Enable = 'on';
    console(fig, ' ');
    console(fig, sprintf('==== RUN started %s ====', datestr(now, 'yyyy-mm-dd HH:MM:SS')));
    if isempty(folders), console(fig, 'Scope: whole tree');
    else, console(fig, ['Scope: ' strjoin(folders, ', ')]); end
    console(fig, sprintf('Conditioning: SVD %s | binning %s (%d) | filter %s | drift %s | normalize %s | ensemble %s | motion %s | FOV %g mm | maps PDF %s', ...
        h.SVD.Value, h.Binning.Value, co.BinSize, h.Filter.Value, h.Drift.Value, h.Normalize.Value, h.Ensemble.Value, h.Motion.Value, eo.FOV_mm, tern(h.MapsPDF.Value, 'yes', 'no')));
    drawnow;

    T = table();
    try
        T = cadence_batch_pipeline(raw, proc, met, args{:});
        if h.DryRun.Value
            console(fig, 'DRY RUN COMPLETED: nothing was processed.');
        elseif fig.UserData.stop
            console(fig, sprintf('RUN STOPPED: %d recording(s) in the medians table so far.', height(T)));
        else
            console(fig, sprintf('PIPELINE COMPLETED SUCCESSFULLY: %d recording(s) in %s', ...
                height(T), fullfile(met, 'CADENCE_median_summary.xlsx')));
        end
    catch ME
        console(fig, ['ERROR: ' ME.message]);
    end
    fig.UserData.last_table = T;
    fig.UserData.running = false;
    fig.UserData.stop = false;
    h.RunButton.Enable = 'on';
    h.StopButton.Enable = 'off';
end

function v = tern(c, a, b)
    if c, v = a; else, v = b; end
end

function tf = should_stop(fig)
    drawnow limitrate;                           % let a STOP click through
    tf = isvalid(fig) && fig.UserData.stop;
end

function console(fig, msg)
% Newest line on top, as in the other CADENCE modules; capped so a long batch
% does not slow the window down (the full log is in cadence_batch_log.txt).
    if ~isvalid(fig), return; end
    ta = fig.UserData.h.Console;
    existing = ta.Value;
    if ~iscell(existing), existing = cellstr(existing); end
    if numel(existing) > 400, existing = existing(1:400); end
    ta.Value = [{char(msg)}; existing(:)];
    drawnow limitrate;
end
