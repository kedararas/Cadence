function cadence_build_summary(T, out_xlsx, varargin)
%CADENCE_BUILD_SUMMARY  Excel workbook of per-recording metric medians.
%
%   cadence_build_summary(medians_csv, out_xlsx)
%   cadence_build_summary(medians_table, out_xlsx)
%
%   Sheets
%     Recordings     one row per recording, grouped by ZT then experiment
%                    (rat), then condition, then pacing CL (longest first)
%     Per camera     one row per camera of each recording (only written when
%                    'PerCamera' is given and some recording has more than one
%                    camera of the same signal type, e.g. a four-camera
%                    voltage rig); same order as Recordings, then camera
%     By experiment  one row per ZT x experiment x condition x CL: median of
%                    the recordings at that CL (usually one) + N
%     By ZT          one row per ZT x condition x CL: median across the
%                    experiments at that ZT + N experiments
%     Alternans onset      per experiment: onset CL of V / Ca alternans, which
%                          came first, concordant -> discordant -> arrhythmia
%                          transition, functional refractory period, and the
%                          same onsets read from the file-name tags
%     Alternans per rec.   the per-recording classification behind it
%     Circadian cosinor    24 h cosinor fit of every metric (per-experiment
%                          value at RefCL, Baseline) and of the onset measures
%     <metric>...    one sheet per metric: rows = CL, columns = experiments
%                    (ZT<n>_R<n>_<condition>) ordered by ZT, followed by one
%                    "ZT<n> <condition> median" column per ZT
%     README         metric definitions and the settings used
%
%   CSV copies of the tables are written next to the workbook.
%
%   Options: 'RefCL' (default 150) reference CL for the cosinor of steady-state
%   metrics; 'OnsetOpts' struct passed to cadence_alternans_onset;
%   'PerCamera' table (or CSV) from cadence_camera_medians.

    p = inputParser;
    p.addParameter('RefCL', 150);
    p.addParameter('OnsetOpts', struct());
    p.addParameter('PerCamera', []);
    p.parse(varargin{:});
    opt = p.Results;

    if ischar(T) || isstring(T)
        T = readtable(char(T), 'TextType', 'string', 'Delimiter', ',');
    end
    out_xlsx = char(out_xlsx);
    [L, arr_names] = cadence_metric_labels();
    have = ismember(L(:,1), T.Properties.VariableNames);
    L = L(have, :);
    % Arrhythmia recordings (Rhythm column, written by the batch) get their own
    % sheet; every other sheet summarizes the paced recordings and leaves the
    % arrhythmia-only metrics out.  The onset analysis still sees all rows.
    Tall = T;
    if ismember('Rhythm', T.Properties.VariableNames)
        is_arr = lower(string(T.Rhythm)) == "arrhythmia";
        is_arr(ismissing(is_arr)) = false;
        T = T(~is_arr, :);
    else
        is_arr = false(height(T), 1);
    end
    Larr = L(ismember(L(:,1), [arr_names; {'df_v'; 'ri_v'; 'oi_v'; 'pacing_hz'; 'capture_ratio'; 'tissue_frac_v'; 'v_polarity_conf'}]), :);
    L = L(~ismember(L(:,1), arr_names), :);
    mnames = L(:,1);  mlabels = L(:,2);

    % ---- normalise meta columns --------------------------------------------
    T.Condition = fillstr(T, 'Condition');
    T.Tag       = fillstr(T, 'Tag');
    T.Experiment = fillstr(T, 'Experiment');
    expnum = nan(height(T), 1);
    for i = 1:height(T)
        t = regexp(char(T.Experiment(i)), '\d+', 'match', 'once');
        if ~isempty(t), expnum(i) = str2double(t); end
    end
    T.ExpNum = expnum;
    condrank = double(T.Condition ~= "Baseline");     % Baseline first
    [~, ord] = sortrows([nz(T.ZT), nz(T.ExpNum), condrank, -nz(T.CL_ms), nz(T.Run)]);
    T = T(ord, :);

    if isfile(out_xlsx), delete(out_xlsx); end
    [outdir, stem] = fileparts(out_xlsx);

    % ---- Sheet 1: Recordings ----------------------------------------------------
    meta_cols = {'ZT','Experiment','Condition','CL_ms','Tag','Run','Date','File'};
    meta_lbl  = {'ZT','Experiment','Condition','CL (ms)','Tag','Run','Date','File'};
    C = [meta_lbl, mlabels', {'Beats/Errors'}];
    body = cell(height(T), numel(C));
    for i = 1:height(T)
        for k = 1:numel(meta_cols), body{i,k} = cellval(T.(meta_cols{k})(i)); end
        for k = 1:numel(mnames), body{i,numel(meta_cols)+k} = cellval(T.(mnames{k})(i)); end
        e = "";
        if ismember('Errors', T.Properties.VariableNames), e = T.Errors(i); end
        body{i,end} = cellval(e);
    end
    writecell([C; body], out_xlsx, 'Sheet', 'Recordings');
    Trec = T(:, [meta_cols, mnames']);
    writetable(Trec, fullfile(outdir, [stem '_recordings.csv']));

    % ---- Sheet 1b: Per camera (multi-camera rigs) -------------------------------
    n_percam = write_per_camera(opt.PerCamera, out_xlsx, outdir, stem);

    % ---- Sheet 2: By experiment -----------------------------------------------------
    keyE = strcat(string(nz(T.ZT)), "|", T.Experiment, "|", T.Condition, "|", string(nz(T.CL_ms)));
    [gE, firstE] = unique_groups(keyE);
    E = cell(numel(firstE), 5 + numel(mnames));
    for g = 1:numel(firstE)
        rows = find(gE == g);
        i0 = rows(1);
        E(g, 1:5) = {cellval(T.ZT(i0)), cellval(T.Experiment(i0)), cellval(T.Condition(i0)), cellval(T.CL_ms(i0)), numel(rows)};
        for k = 1:numel(mnames)
            E{g, 5+k} = cellval(nanmed(T.(mnames{k})(rows)));
        end
    end
    hdrE = [{'ZT','Experiment','Condition','CL (ms)','N recordings'}, mlabels'];
    writecell([hdrE; E], out_xlsx, 'Sheet', 'By experiment');
    writecell([hdrE; E], fullfile(outdir, [stem '_by_experiment.csv']));

    % ---- Sheet 3: By ZT  (median across experiments) --------------------------------
    % first collapse duplicates within an experiment (E above), then pool.
    Ezt  = string(cellfun(@(x) num2str(x), E(:,1), 'UniformOutput', false));
    Econ = string(E(:,3));
    Ecl  = string(cellfun(@(x) num2str(x), E(:,4), 'UniformOutput', false));
    keyZ = strcat(Ezt, "|", Econ, "|", Ecl);
    [gZ, firstZ] = unique_groups(keyZ);
    Z = cell(numel(firstZ), 4 + numel(mnames));
    for g = 1:numel(firstZ)
        rows = find(gZ == g);
        i0 = rows(1);
        Z(g, 1:4) = {E{i0,1}, E{i0,3}, E{i0,4}, numel(rows)};
        for k = 1:numel(mnames)
            v = cellfun(@(x) tonum(x), E(rows, 5+k));
            Z{g, 4+k} = cellval(nanmed(v));
        end
    end
    hdrZ = [{'ZT','Condition','CL (ms)','N experiments'}, mlabels'];
    writecell([hdrZ; Z], out_xlsx, 'Sheet', 'By ZT');
    writecell([hdrZ; Z], fullfile(outdir, [stem '_by_ZT.csv']));

    % ---- alternans onset / circadian ----------------------------------------------
    O = table(); S = table(); Cz = table();
    try
        oo = struct2nv(opt.OnsetOpts);
        exf = fullfile(outdir, 'cadence_excluded_recordings.csv');
        if isfile(exf) && ~isfield(opt.OnsetOpts, 'Excluded'), oo = [oo, {'Excluded', exf}]; end
        [O, S] = cadence_alternans_onset(Tall, oo{:});
        writetable(O, out_xlsx, 'Sheet', 'Alternans onset');
        writetable(S, out_xlsx, 'Sheet', 'Alternans per recording');
        writetable(O, fullfile(outdir, [stem '_alternans_onset.csv']));
        writetable(S, fullfile(outdir, [stem '_alternans_per_recording.csv']));
    catch ME
        warning('cadence_build_summary:onset', 'Alternans onset sheet skipped: %s', ME.message);
    end
    try
        Cz = cadence_cosinor(T, O, 'RefCL', opt.RefCL);
        writetable(Cz, out_xlsx, 'Sheet', 'Circadian cosinor');
        writetable(Cz, fullfile(outdir, [stem '_circadian_cosinor.csv']));
    catch ME
        warning('cadence_build_summary:cosinor', 'Circadian sheet skipped: %s', ME.message);
    end

    % ---- per-metric sheets: rows = CL, cols = experiments then ZT medians ----------
    expid = strcat("ZT", string(nz(T.ZT)), "_", T.Experiment, "_", T.Condition);
    expkey = [nz(T.ZT), nz(T.ExpNum), condrank];
    [ukey, ia] = unique(expkey, 'rows');           % sorted by ZT, rat, condition
    exps = expid(ia);
    ztid  = strcat("ZT", string(ukey(:,1)), " ", T.Condition(ia), " median");
    [uzt, iz] = unique([ukey(:,1), ukey(:,3)], 'rows');
    ztcols = ztid(iz);
    cls = unique(T.CL_ms(isfinite(T.CL_ms)));
    cls = sort(cls, 'descend');
    used = {};
    for k = 1:numel(mnames)
        M = nan(numel(cls), numel(exps));
        for e = 1:numel(exps)
            sel = expid == exps(e);
            for r = 1:numel(cls)
                v = T.(mnames{k})(sel & T.CL_ms == cls(r));
                M(r, e) = nanmed(v);
            end
        end
        Mz = nan(numel(cls), size(uzt,1));
        for z = 1:size(uzt,1)
            cols = ukey(:,1) == uzt(z,1) & ukey(:,3) == uzt(z,2);
            Mz(:, z) = nanmed_rows(M(:, cols));
        end
        hdr = [{'CL (ms)'}, cellstr(exps'), cellstr(ztcols')];
        body = [num2cell(cls), num2cell_nan([M, Mz])];
        [sname, used] = sheet_name(mlabels{k}, used);
        writecell([hdr; body], out_xlsx, 'Sheet', sname);
    end

    % ---- Arrhythmia: one row per arrhythmia recording, its own metric set -------------
    if any(is_arr)
        A = Tall(is_arr, :);
        meta = intersect({'ZT','Experiment','Condition','CL_ms','Tag','Run','Date','File', ...
                          'Rhythm','Capture','Rhythm_basis','Name_tag'}, A.Properties.VariableNames, 'stable');
        A = A(:, [meta, Larr(:,1)']);
        hdr = [meta, Larr(:,2)'];
        writecell([hdr; table2cell(A)], out_xlsx, 'Sheet', 'Arrhythmia');
        writetable(A, fullfile(outdir, [stem '_arrhythmia.csv']));
    end

    % ---- README ---------------------------------------------------------------------
    R = {'CADENCE median summary', ''; ...
         'Generated', datestr(now, 'yyyy-mm-dd HH:MM:SS'); ...
         'Recordings', num2str(height(T)); ...
         'Value', 'median over finite pixels of the masked map (slot 1), per recording'; ...
         'Mask', 'adaptive SNR mask per camera: SNR >= max(1.5, 0.5 x median tissue SNR), holes filled, specks removed'; ...
         'Analysis window', 'ensemble-averaged beat (CAM<n>_average) for all timing metrics; full record for alternans and DF/RI/OI'; ...
         'Cameras', 'voltage metrics from the first voltage camera, calcium metrics from the first calcium camera, as assigned for the run (camera_roles in each metrics file); files extracted without roles assume CAM1 = voltage, CAM2 = calcium'; ...
         'Per camera', per_camera_note(n_percam); ...
         'Conduction velocity', 'interior pixels only (8 px erosion of the tissue mask), cm/s'; ...
         'Alternans present', 'fraction of tissue pixels with alternans ratio > 0.10 and p < 0.05 is >= 0.10 (SigFrac); discordant when concordance ratio < 0.80'; ...
         'Rhythm', 'each recording is classified before extraction (Rhythm column; Rhythm_basis says why): paced with DF within 5% of the stimulus rate = 1:1 capture -> paced metrics; any other capture, or no stimulus with an irregular rhythm (beat CV > 0.15 or OI < 0.35), or no stimulus, regular, and an arrhythmia file tag -> arrhythmia (voltage only: DF/RI/OI, wavefront and rotor dynamics, Arrhythmia sheet). Overrides: cadence_rhythm_overrides.csv in the metrics folder'; ...
         'Arrhythmia (onset)', 'classified arrhythmia, a file tag containing "arrhythm", or capture ratio outside 1 +/- 0.15 with OI (V) < 0.5'; ...
         'Functional refractory period', 'shortest CL still at 1:1 capture (|DF/pacing - 1| <= 0.15)'; ...
         'Circadian cosinor', sprintf('y = M + A cos(2 pi t/24) + B sin(2 pi t/24) on per-experiment values at CL %g ms (Baseline); p = F-test vs flat mean', opt.RefCL); ...
         '', ''; ...
         'Column', 'Definition'};
    R = [R; L(:, 2:3); Larr(ismember(Larr(:,1), arr_names), 2:3)];
    writecell(R, out_xlsx, 'Sheet', 'README');
    fprintf('Wrote %s (%d recordings, %d metrics)\n', out_xlsx, height(T), numel(mnames));
end


% =========================================================================
function n = write_per_camera(Tc, out_xlsx, outdir, stem)
% 'Per camera' sheet + CSV.  Returns the number of rows written (0 = no sheet).
    n = 0;
    if isempty(Tc), return; end
    if ischar(Tc) || isstring(Tc)
        if ~isfile(char(Tc)), return; end
        Tc = readtable(char(Tc), 'TextType', 'string', 'Delimiter', ',');
    end
    if ~istable(Tc) || height(Tc) == 0 || ~all(ismember({'Camera','Signal'}, Tc.Properties.VariableNames))
        return;
    end
    Tc.Camera = fillstr(Tc, 'Camera');  Tc.Signal = fillstr(Tc, 'Signal');
    key = fillstr(Tc, 'File');
    if ismember('Source', Tc.Properties.VariableNames), key = fillstr(Tc, 'Source'); end
    % only worth a sheet when some recording has two cameras of one signal type
    [~, ~, g] = unique(strcat(key, "|", Tc.Signal));
    if max(accumarray(g, 1)) < 2, return; end

    Tc.Condition  = fillstr(Tc, 'Condition');
    Tc.Experiment = fillstr(Tc, 'Experiment');
    expnum = nan(height(Tc), 1);  camnum = nan(height(Tc), 1);
    for i = 1:height(Tc)
        t = regexp(char(Tc.Experiment(i)), '\d+', 'match', 'once');
        if ~isempty(t), expnum(i) = str2double(t); end
        t = regexp(char(Tc.Camera(i)), '\d+', 'match', 'once');
        if ~isempty(t), camnum(i) = str2double(t); end
    end
    num = @(c) numcol(Tc, c);
    [~, ord] = sortrows([nz(num('ZT')), nz(expnum), double(Tc.Condition ~= "Baseline"), ...
                         -nz(num('CL_ms')), nz(num('Run')), double(Tc.Signal ~= "voltage"), nz(camnum)]);
    Tc = Tc(ord, :);

    L = cadence_metric_labels();
    L = L(ismember(L(:,1), Tc.Properties.VariableNames), :);
    has = cellfun(@(c) any(isfinite(double(Tc.(c)))), L(:,1));
    L = L(has, :);                                  % drop metrics no camera has
    meta_cols = {'ZT','Experiment','Condition','CL_ms','Tag','Run','Date','File','Camera','Signal'};
    meta_lbl  = {'ZT','Experiment','Condition','CL (ms)','Tag','Run','Date','File','Camera','Signal'};
    keep = ismember(meta_cols, Tc.Properties.VariableNames);
    meta_cols = meta_cols(keep);  meta_lbl = meta_lbl(keep);
    C = [meta_lbl, L(:,2)'];
    body = cell(height(Tc), numel(C));
    for i = 1:height(Tc)
        for k = 1:numel(meta_cols), body{i,k} = cellval(Tc.(meta_cols{k})(i)); end
        for k = 1:size(L,1), body{i,numel(meta_cols)+k} = cellval(Tc.(L{k,1})(i)); end
    end
    writecell([C; body], out_xlsx, 'Sheet', 'Per camera');
    writetable(Tc(:, [meta_cols, L(:,1)']), fullfile(outdir, [stem '_per_camera.csv']));
    n = height(Tc);
end

function s = per_camera_note(n)
    if n > 0
        s = sprintf(['%d rows: every camera of every recording, voltage metrics for voltage cameras and ' ...
                     'calcium metrics for calcium cameras, same medians as Recordings; pair metrics ' ...
                     '(V-Ca delay, Ca-AP coupling, risk) and recording-level metrics stay on Recordings'], n);
    else
        s = 'not written: no recording has more than one camera of the same signal type';
    end
end

function v = numcol(T, c)
    if ismember(c, T.Properties.VariableNames), v = double(T.(c)); else, v = nan(height(T), 1); end
end

function nv = struct2nv(s)
    f = fieldnames(s); nv = cell(1, 2*numel(f));
    for k = 1:numel(f), nv{2*k-1} = f{k}; nv{2*k} = s.(f{k}); end
end

function s = fillstr(T, col)
    if ismember(col, T.Properties.VariableNames)
        s = string(T.(col)); s(ismissing(s)) = "";
    else
        s = strings(height(T), 1);
    end
end

function v = nz(v)
    v = double(v); v(isnan(v)) = -1;
end

function m = nanmed(x)
    x = double(x(:)); x = x(isfinite(x));
    if isempty(x), m = NaN; else, m = median(x); end
end

function m = nanmed_rows(M)
    m = nan(size(M,1), 1);
    for r = 1:size(M,1), m(r) = nanmed(M(r,:)); end
end

function c = cellval(v)
    if isstring(v) || ischar(v)
        c = char(v);
    elseif isnumeric(v) || islogical(v)
        if isempty(v) || (isscalar(v) && isnan(v)), c = []; else, c = double(v); end
    else
        c = [];
    end
end

function v = tonum(c)
    if isempty(c), v = NaN; else, v = double(c); end
end

function C = num2cell_nan(M)
    C = num2cell(M);
    C(isnan(M)) = {[]};
end

function [g, first] = unique_groups(keys)
    [~, first, g] = unique(keys, 'stable');
end

function [name, used] = sheet_name(label, used)
    name = regexprep(label, '[\[\]\:\*\?\/\\]', '');
    name = strtrim(regexprep(name, '\s+', ' '));
    if numel(name) > 31, name = name(1:31); end
    base = name; k = 2;
    while any(strcmpi(used, name))
        suffix = sprintf(' %d', k);
        name = [base(1:min(end, 31 - numel(suffix))) suffix]; k = k + 1;
    end
    used{end+1} = name;
end
