function cadence_build_summary(T, out_xlsx, varargin)
%CADENCE_BUILD_SUMMARY  Excel workbook of per-recording metric medians.
%
%   cadence_build_summary(medians_csv, out_xlsx)
%   cadence_build_summary(medians_table, out_xlsx)
%
%   Sheets
%     Recordings     one row per recording, grouped by ZT then experiment
%                    (rat), then condition, then pacing CL (longest first)
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
%   metrics; 'OnsetOpts' struct passed to cadence_alternans_onset.

    p = inputParser;
    p.addParameter('RefCL', 150);
    p.addParameter('OnsetOpts', struct());
    p.parse(varargin{:});
    opt = p.Results;

    if ischar(T) || isstring(T)
        T = readtable(char(T), 'TextType', 'string', 'Delimiter', ',');
    end
    out_xlsx = char(out_xlsx);
    L = cadence_metric_labels();
    have = ismember(L(:,1), T.Properties.VariableNames);
    L = L(have, :);
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
        [O, S] = cadence_alternans_onset(T, oo{:});
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

    % ---- README ---------------------------------------------------------------------
    R = {'CADENCE median summary', ''; ...
         'Generated', datestr(now, 'yyyy-mm-dd HH:MM:SS'); ...
         'Recordings', num2str(height(T)); ...
         'Value', 'median over finite pixels of the masked map (slot 1), per recording'; ...
         'Mask', 'adaptive SNR mask per camera: SNR >= max(1.5, 0.5 x median tissue SNR), holes filled, specks removed'; ...
         'Analysis window', 'ensemble-averaged beat (CAM<n>_average) for all timing metrics; full record for alternans and DF/RI/OI'; ...
         'Cameras', 'CAM1 = voltage, CAM2 = calcium'; ...
         'Conduction velocity', 'interior pixels only (8 px erosion of the tissue mask), cm/s'; ...
         'Alternans present', 'fraction of tissue pixels with alternans ratio > 0.10 and p < 0.05 is >= 0.10 (SigFrac); discordant when concordance ratio < 0.80'; ...
         'Arrhythmia', 'file tag containing "arrhythm", or capture ratio outside 1 +/- 0.15 with OI (V) < 0.5'; ...
         'Functional refractory period', 'shortest CL still at 1:1 capture (|DF/pacing - 1| <= 0.15)'; ...
         'Circadian cosinor', sprintf('y = M + A cos(2 pi t/24) + B sin(2 pi t/24) on per-experiment values at CL %g ms (Baseline); p = F-test vs flat mean', opt.RefCL); ...
         '', ''; ...
         'Column', 'Definition'};
    R = [R; L(:, 2:3)];
    writecell(R, out_xlsx, 'Sheet', 'README');
    fprintf('Wrote %s (%d recordings, %d metrics)\n', out_xlsx, height(T), numel(mnames));
end


% =========================================================================
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
