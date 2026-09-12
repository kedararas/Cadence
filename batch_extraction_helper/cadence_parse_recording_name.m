function meta = cadence_parse_recording_name(filepath)
%CADENCE_PARSE_RECORDING_NAME  Study metadata from a recording's path.
%
%   meta = cadence_parse_recording_name('/.../ZT2_8am/R2/05-Rat-R2-20250511-8AM-Baseline-200ms-3v.mat')
%
%   File name pattern (any piece may be missing; missing -> NaN / ""):
%     <run>-Rat-R<rat>-<yyyymmdd>-<clock>-<Baseline|IP>-<CL>ms-<stim>v[-<tag>]
%   Folder pattern: <root>/ZT<n>[_anything]/R<n>/file
%
%   Fields
%     file, base        file name; base without .mat/-metrics/-conditioned
%     zt_folder, rat_folder   the two parent folder names
%     ZT                ZT number from the ZT folder (NaN if not parsable)
%     experiment        "R<n>" from the rat folder, else from the file name
%     run, rat, date, clock, condition, CL_ms, stim, tag

    filepath = char(filepath);
    [folder, base, ext] = fileparts(filepath);
    meta = struct();
    meta.file = [base ext];
    base = regexprep(base, '(-metrics|-conditioned)$', '');
    meta.base = base;

    parts = strsplit(folder, filesep);
    parts = parts(~cellfun(@isempty, parts));
    meta.rat_folder = ""; meta.zt_folder = "";
    if numel(parts) >= 1, meta.rat_folder = string(parts{end}); end
    if numel(parts) >= 2, meta.zt_folder = string(parts{end-1}); end

    % ZT: nearest ancestor folder named ZT<n>..., else the file name
    meta.ZT = NaN;
    for q = numel(parts):-1:1
        t = regexp(parts{q}, '^ZT(\d+)', 'tokens', 'once');
        if ~isempty(t), meta.ZT = str2double(t{1}); meta.zt_folder = string(parts{q}); break; end
    end
    if isnan(meta.ZT)
        t = regexp(base, '[_-]ZT(\d+)', 'tokens', 'once');
        if ~isempty(t), meta.ZT = str2double(t{1}); end
    end
    % experiment: nearest ancestor folder named R<n> (rat / heart id)
    exp_folder = "";
    for q = numel(parts):-1:1
        if ~isempty(regexp(parts{q}, '^R\d+$', 'once')), exp_folder = string(parts{q}); break; end
    end

    meta.run  = tok(base, '^(\d+)[-_]');
    meta.rat  = tok(base, '-R(\d+)-');
    if isnan(meta.rat), meta.rat = tok(base, 'Rat-R(\d+)'); end
    meta.date = tokstr(base, '-(\d{8})-');
    meta.clock = tokstr(base, '-(\d{1,2}(?::\d{2})?\s?[AaPp][Mm])-');
    if ~isempty(regexp(base, '-IP-', 'once'))
        meta.condition = "IP";
    elseif ~isempty(regexp(base, '-Baseline-', 'once', 'ignorecase'))
        meta.condition = "Baseline";
    else
        meta.condition = "";
    end
    meta.CL_ms = tok(base, '[-_](\d+)ms');
    meta.stim  = tokstr(base, '[-_](\d+(?:\.\d+)?[vV])(?:[-_]|$)');
    t = regexp(base, '[-_]\d+(?:\.\d+)?[vV][-_](.+)$', 'tokens', 'once');
    if isempty(t), meta.tag = ""; else, meta.tag = string(t{1}); end

    if strlength(exp_folder) > 0
        meta.experiment = exp_folder;
    elseif ~isnan(meta.rat)
        meta.experiment = string(sprintf('R%d', meta.rat));
    else
        meta.experiment = meta.rat_folder;
    end
end

function v = tok(s, pat)
    t = regexp(s, pat, 'tokens', 'once');
    if isempty(t), v = NaN; else, v = str2double(t{1}); end
end

function v = tokstr(s, pat)
    t = regexp(s, pat, 'tokens', 'once');
    if isempty(t), v = ""; else, v = string(t{1}); end
end
