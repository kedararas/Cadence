function meta = cadence_parse_recording_name(filepath, lights_on)
%CADENCE_PARSE_RECORDING_NAME  Study metadata from a recording's path.
%
%   meta = cadence_parse_recording_name('/.../ZT2_8am/R2/05-Rat-R2-20250511-8AM-Baseline-200ms-3v.mat')
%   meta = cadence_parse_recording_name(path, lights_on)   % lights-on hour, default 6
%
%   File name pattern (any piece may be missing; missing -> NaN / ""):
%     <run>-Rat-R<rat>-<yyyymmdd>-<clock>-<Baseline|IP>-<CL>ms-<stim>v[-<tag>]
%   The older underscore form is read too:
%     <run>_<m>_<d>_<yyyy>_Rat_WH_[CT<n>_]<clock>_R<rat>_<IP|Baseline>_..._<CL>ms_<stim>v
%   Folder pattern: <root>/ZT<n>[_anything]/R<n>/file
%
%   Fields
%     file, base        file name; base without .mat/-metrics/-conditioned
%     zt_folder, rat_folder   the two parent folder names
%     ZT                Zeitgeber time.  An explicit label wins: a ZT<n> folder,
%                       else ZT<n> in the file name.  Otherwise it is derived
%                       from the clock in the label, ZT = clock hour - lights_on
%                       (mod 24): with lights on at 06:00, 8AM -> ZT2, 12PM ->
%                       ZT6, 4PM -> ZT10, 8PM -> ZT14, 12AM -> ZT18, 4AM -> ZT22.
%                       NaN when the label has neither.
%     zt_source         "folder" | "file name" | "clock (lights on HH:00)" | ""
%     experiment        "R<n>" from the rat folder, else from the file name
%     run, rat, date, clock, condition, CL_ms, stim, tag
%     clock_hour        the labelled clock as hours since midnight (NaN if none)
%     acquired          labelled acquisition date and time, "yyyy-mm-dd HH:MM"
%                       (date only, or "", when the label lacks a piece).  This
%                       is the NOMINAL time in the file label, not the camera's
%                       own time stamp.

    if nargin < 2 || isempty(lights_on), lights_on = 6; end
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
    meta.ZT = NaN;  meta.zt_source = "";
    for q = numel(parts):-1:1
        t = regexp(parts{q}, '^ZT(\d+)', 'tokens', 'once');
        if ~isempty(t)
            meta.ZT = str2double(t{1}); meta.zt_folder = string(parts{q}); meta.zt_source = "folder"; break;
        end
    end
    if isnan(meta.ZT)
        t = regexp(base, '[_-](ZT|CT)(\d+)', 'tokens', 'once');      % CT<n>: circadian time, same scale
        if ~isempty(t), meta.ZT = str2double(t{2}); meta.zt_source = "file name"; end
    end
    % experiment: nearest ancestor folder named R<n> (rat / heart id)
    exp_folder = "";
    for q = numel(parts):-1:1
        if ~isempty(regexp(parts{q}, '^R\d+$', 'once')), exp_folder = string(parts{q}); break; end
    end

    meta.run  = tok(base, '^(\d+)[-_]');
    meta.rat  = tok(base, '[-_]R(\d+)(?:[-_]|$)');
    if isnan(meta.rat), meta.rat = tok(base, 'Rat-R(\d+)'); end
    meta.date = tokstr(base, '-(\d{8})-');
    if strlength(meta.date) == 0                      % older form: <run>_<m>_<d>_<yyyy>_
        t = regexp(base, '^\d+_(\d{1,2})_(\d{1,2})_(\d{4})_', 'tokens', 'once');
        if ~isempty(t), meta.date = string(sprintf('%s%02d%02d', t{3}, str2double(t{1}), str2double(t{2}))); end
    end
    meta.clock = tokstr(base, '[-_](\d{1,2}(?::\d{2})?\s?[AaPp][Mm])(?:[-_]|$)');
    % clock as hours since midnight, the labelled acquisition time, and ZT from
    % the clock when the label carries no explicit ZT
    meta.clock_hour = NaN;
    ck = char(meta.clock);
    hh = tok(ck, '^(\d{1,2})');  mi = tok(ck, ':(\d{2})');  ap = tokstr(ck, '([AaPp])[Mm]$');
    if ~isnan(hh) && strlength(ap) > 0
        if isnan(mi), mi = 0; end
        meta.clock_hour = mod(hh, 12) + 12 * strcmpi(ap, "p") + mi / 60;
    end
    meta.acquired = "";
    dd = char(meta.date);
    if numel(dd) == 8
        meta.acquired = string(sprintf('%s-%s-%s', dd(1:4), dd(5:6), dd(7:8)));
        if ~isnan(meta.clock_hour)
            meta.acquired = meta.acquired + sprintf(' %02d:%02d', floor(meta.clock_hour), round(60 * mod(meta.clock_hour, 1)));
        end
    end
    if isnan(meta.ZT) && ~isnan(meta.clock_hour)
        meta.ZT = mod(meta.clock_hour - lights_on, 24);
        meta.zt_source = string(sprintf('clock (lights on %02d:00)', lights_on));
    end
    if ~isempty(regexp(base, '[-_]IP(?:[-_]|$)', 'once'))
        meta.condition = "IP";
    elseif ~isempty(regexp(base, '[-_]Baseline(?:[-_]|$)', 'once', 'ignorecase'))
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
