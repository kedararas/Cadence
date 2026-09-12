function [O, S] = cadence_alternans_onset(T, varargin)
%CADENCE_ALTERNANS_ONSET  Alternans onset, concordance and arrhythmia per experiment.
%
%   [O, S] = cadence_alternans_onset(T)             T = medians table / CSV
%   [O, S] = cadence_alternans_onset(T, 'SigFrac', 0.10, 'ConcordanceMin', 0.80, ...)
%
%   Walks each experiment (ZT x Experiment x Condition) down its pacing
%   protocol (longest CL first) and classifies every recording, then reports
%   where things started:
%
%     Onset_CL_V / Onset_CL_Ca   longest CL at which voltage / calcium
%                                alternans is present (and still present at
%                                the next shorter CL, unless it is the last)
%     First                      "Ca first" | "V first" | "simultaneous" |
%                                "V only" | "Ca only" | "none"
%     Onset_CL_discordant        longest CL with V alternans whose
%                                concordance ratio < ConcordanceMin
%     Arrhythmia_CL              CL of the FIRST (lowest run number) recording
%                                flagged arrhythmic: a recording
%                                listed in the excluded-recordings table
%                                ('Excluded', the arrhythmia runs that were not
%                                extracted), a file tag containing "arrhythm",
%                                or capture ratio / organization index outside
%                                the paced range
%     FRP_CL                     functional refractory period: shortest CL
%                                still at ~1:1 capture (|DF/pacing - 1| <= CaptureTol)
%     Loss_of_capture_CL         longest CL below FRP_CL with lost capture
%     Transition                 e.g. "Con@85 -> Dis@70 -> Arrhythmia@65"
%     Sequence                   every CL with its state, e.g.
%                                "200:- 175:- ... 85:V+Ca(con) 80:V+Ca(dis) 65:ARR"
%     Annotated_*_CL             the same onsets read from the file-name tags
%                                (Alternans / Dis / Arrhythmia), for comparison
%
%   S is the per-recording classification behind O.
%
%   Criteria (options)
%     'SigFrac'         alternans present when the fraction of tissue pixels
%                       with ratio > 0.10 and p < 0.05 is >= this (default 0.10)
%     'ConcordanceMin'  discordant when concordance ratio < this (default 0.80)
%     'CaptureTol'      1:1 capture when |capture ratio - 1| <= this (default 0.15)
%     'OIMin'           arrhythmic by metric when OI (V) < this AND capture is
%                       lost (default 0.50).  Independently, a recording whose
%                       DF exceeds the pacing rate by more than CaptureTol
%                       (capture ratio > 1 + CaptureTol) is flagged arrhythmic:
%                       the tissue is running faster than the stimulus.  Ratios
%                       below 1 are block (2:1 etc.), not arrhythmia.
%     'RequireConsecutive'  metric-based onsets (V / Ca alternans, concordant,
%                       discordant) must also hold at the next shorter CL,
%                       unless it is the shortest CL recorded (default true).
%                       Suppresses single borderline recordings.  Arrhythmic
%                       recordings are removed from the sequence first, so an
%                       arrhythmia between two alternating recordings does not
%                       break the test.  Arrhythmia itself and the file-name
%                       annotations are taken as given.

    p = inputParser;
    p.addParameter('SigFrac', 0.10);
    p.addParameter('ConcordanceMin', 0.80);
    p.addParameter('CaptureTol', 0.15);
    p.addParameter('OIMin', 0.50);
    p.addParameter('RequireConsecutive', true);
    p.addParameter('Excluded', []);      % table/CSV of recordings not extracted (arrhythmia runs)
    p.parse(varargin{:});
    o = p.Results;

    if ischar(T) || isstring(T), T = readtable(char(T), 'TextType', 'string', 'Delimiter', ','); end
    X = o.Excluded;
    if ischar(X) || isstring(X)
        if isfile(char(X)), X = readtable(char(X), 'TextType', 'string', 'Delimiter', ','); else, X = []; end
    end
    n = height(T);
    col = @(name) getcol(T, name, n);

    ZT = col('ZT'); Exp = strcol(T, 'Experiment', n); Cond = strcol(T, 'Condition', n);
    CL = col('CL_ms'); Run = col('Run'); Tag = lower(strcol(T, 'Tag', n));
    sigV = col('alt_sig_frac_v'); sigCa = col('ca_alt_sig_frac');
    conc = col('concordance_v'); cap = col('capture_ratio'); oi = col('oi_v');

    % ---- per-recording classification -----------------------------------------------
    V_alt  = sigV  >= o.SigFrac;
    Ca_alt = sigCa >= o.SigFrac;
    disc   = V_alt & conc < o.ConcordanceMin;
    % Tag matching tolerates the spelling variants present in the corpus
    % ("Arhhythmia", "alternan"), and must not confuse the two with each other.
    tag_alt = ~cellfun(@isempty, regexp(cellstr(Tag), 'alternan', 'once'));
    tag_dis = contains(Tag, 'dis');
    tag_arr = ~cellfun(@isempty, regexp(cellstr(Tag), 'a[r]{1,2}[rh]?h?y+t?h?m', 'once'));
    capt_ok = isfinite(cap) & abs(cap - 1) <= o.CaptureTol;
    capt_lost = isfinite(cap) & ~capt_ok;
    df_exceeds = isfinite(cap) & cap > 1 + o.CaptureTol;      % tissue faster than the stimulus: not driven
    arr_oi     = capt_lost & oi < o.OIMin;                     % disorganized and not captured
    arr = tag_arr | df_exceeds | arr_oi;
    arr_source = strings(n, 1);
    for i = 1:n
        src = strings(0, 1);
        if tag_arr(i),    src(end+1) = "tag"; end %#ok<AGROW>
        if df_exceeds(i), src(end+1) = "DF>pacing"; end %#ok<AGROW>
        if arr_oi(i),     src(end+1) = "OI low & capture lost"; end %#ok<AGROW>
        arr_source(i) = strjoin(src, "; ");
    end

    state = strings(n, 1);
    for i = 1:n
        if arr(i)
            state(i) = "ARR";
        elseif V_alt(i) && Ca_alt(i)
            state(i) = "V+Ca";
        elseif V_alt(i)
            state(i) = "V";
        elseif Ca_alt(i)
            state(i) = "Ca";
        else
            state(i) = "-";
        end
        if V_alt(i) && ~arr(i)
            if disc(i), state(i) = state(i) + "(dis)"; else, state(i) = state(i) + "(con)"; end
        end
        if capt_lost(i) && ~arr(i), state(i) = state(i) + "!capture"; end
    end

    S = table(ZT, Exp, Cond, CL, Run, sigV, sigCa, V_alt, Ca_alt, conc, disc, cap, capt_ok, df_exceeds, oi, arr, arr_source, ...
              tag_alt, tag_dis, tag_arr, state, ...
              'VariableNames', {'ZT','Experiment','Condition','CL_ms','Run','V_sig_fraction','Ca_sig_fraction', ...
              'V_alternans','Ca_alternans','Concordance','Discordant','Capture_ratio','Capture_1to1','DF_exceeds_pacing','OI_V', ...
              'Arrhythmia','Arrhythmia_source','Tag_alternans','Tag_discordant','Tag_arrhythmia','State'});
    [~, ord] = sortrows([nz(ZT), expnum(Exp), double(Cond ~= "Baseline"), -nz(CL), nz(Run)]);
    S = S(ord, :);

    % ---- per-experiment summary -------------------------------------------------------
    key = strcat(string(nz(S.ZT)), "|", S.Experiment, "|", S.Condition);
    ukey = unique(key, 'stable');
    O = table();
    for g = 1:numel(ukey)
        R = S(key == ukey(g), :);
        R = sortrows(R, 'CL_ms', 'descend');
        cls = R.CL_ms;
        row = struct();
        row.ZT = R.ZT(1); row.Experiment = R.Experiment(1); row.Condition = R.Condition(1);
        row.N_recordings = height(R);
        row.CL_longest = max(cls); row.CL_shortest = min(cls);
        % Arrhythmic recordings are REMOVED from the sequence, not marked as
        % "no alternans": AND-ing ~Arrhythmia into the flag made an arrhythmia
        % between two alternating recordings fail the consecutive test below,
        % pushing the onset down a rung (G1 ZT22 R23: Ca onset read 65 ms
        % where the sequence V+Ca(con) ARR V+Ca(con) V+Ca(con) gives 75 ms).
        keep = ~R.Arrhythmia;
        kcl  = cls(keep);
        row.Onset_CL_V  = onset(kcl, R.V_alternans(keep),  o.RequireConsecutive);
        row.Onset_CL_Ca = onset(kcl, R.Ca_alternans(keep), o.RequireConsecutive);
        row.First = first_of(row.Onset_CL_V, row.Onset_CL_Ca);
        row.Onset_CL_discordant = onset(kcl, R.Discordant(keep), o.RequireConsecutive);
        row.Onset_CL_concordant = onset(kcl, R.V_alternans(keep) & ~R.Discordant(keep), o.RequireConsecutive);
        % Arrhythmia: the FIRST arrhythmia recording in run order (the protocol
        % descends in CL, so the induction CL is the earliest run; later
        % arrhythmia runs are often recorded back at a long CL while the
        % episode continues and must not be read as the onset).
        arr_cl = R.CL_ms(R.Arrhythmia); arr_run = R.Run(R.Arrhythmia);
        row.Sequence_excluded = "";
        if istable(X) && ~isempty(X)
            xs = X(nz(X.ZT) == nz(row.ZT) & string(X.Experiment) == row.Experiment & ...
                   string(X.Condition) == row.Condition, :);
            xs = xs(contains(lower(string(xs.Tag)), 'arrhythm') & isfinite(xs.CL_ms), :);
            if ~isempty(xs)
                arr_cl = [arr_cl; xs.CL_ms]; arr_run = [arr_run; xs.Run];
                [~, ord] = sort(nz(xs.Run));
                row.Sequence_excluded = strjoin(strcat("run", string(xs.Run(ord)), ":", string(xs.CL_ms(ord)), "ms ARR(excluded)")', " ");
            end
        end
        row.Arrhythmia_CL = first_by_run(arr_cl, arr_run);
        row.Peak_V_sig_fraction  = maxnan(R.V_sig_fraction(~R.Arrhythmia));
        row.Peak_Ca_sig_fraction = maxnan(R.Ca_sig_fraction(~R.Arrhythmia));
        % functional refractory period: shortest CL still captured 1:1
        c_ok = R.Capture_1to1 & ~R.Arrhythmia;
        if any(c_ok), row.FRP_CL = min(cls(c_ok)); else, row.FRP_CL = NaN; end
        lost = isfinite(R.Capture_ratio) & ~R.Capture_1to1 & cls < row.FRP_CL;
        if any(lost), row.Loss_of_capture_CL = max(cls(lost)); else, row.Loss_of_capture_CL = NaN; end
        row.Transition = transition(row.Onset_CL_concordant, row.Onset_CL_discordant, row.Arrhythmia_CL);
        row.Sequence = strjoin(strcat(string(cls), ":", R.State)', " ");
        % annotations from the file names
        row.Annotated_alternans_CL  = onset(cls, R.Tag_alternans, false);      % deliberate labels: take as given
        row.Annotated_discordant_CL = onset(cls, R.Tag_discordant, false);
        tcl = R.CL_ms(R.Tag_arrhythmia); trun = R.Run(R.Tag_arrhythmia);
        if istable(X) && ~isempty(X) && exist('xs', 'var') && ~isempty(xs)
            tcl = [tcl; xs.CL_ms]; trun = [trun; xs.Run];
        end
        row.Annotated_arrhythmia_CL = first_by_run(tcl, trun);
        O = [O; struct2table(row, 'AsArray', true)]; %#ok<AGROW>
    end
    if ~isempty(O)
        O.First = string(O.First); O.Transition = string(O.Transition); O.Sequence = string(O.Sequence);
    end
end


% -------------------------------------------------------------------------
function cl = onset(cls, flag, consecutive)
% longest CL where flag holds (and, if consecutive, also at the next shorter CL)
    cl = NaN;
    flag = logical(flag(:)); cls = cls(:);
    for i = 1:numel(cls)
        if ~flag(i), continue; end
        if consecutive && i < numel(cls) && ~flag(i+1), continue; end
        cl = cls(i); return;
    end
end

function cl = first_by_run(cls, runs)
% CL of the earliest recording (smallest run number); shortest CL if runs unknown
    cl = NaN;
    ok = isfinite(cls); cls = cls(ok); runs = runs(ok);
    if isempty(cls), return; end
    if any(isfinite(runs))
        runs(~isfinite(runs)) = Inf;
        [~, i] = min(runs); cl = cls(i);
    else
        cl = min(cls);
    end
end

function s = first_of(v, ca)
    if isnan(v) && isnan(ca), s = "none";
    elseif isnan(ca), s = "V only";
    elseif isnan(v), s = "Ca only";
    elseif ca > v, s = "Ca first";
    elseif v > ca, s = "V first";
    else, s = "simultaneous";
    end
end

function s = transition(con, dis, arr)
% events ordered as the protocol reaches them (longest CL first)
    ev = [con, dis, arr]; lab = ["Con", "Dis", "Arrhythmia"];
    keep = ~isnan(ev); ev = ev(keep); lab = lab(keep);
    [ev, k] = sort(ev, 'descend'); lab = lab(k);
    if isempty(ev), s = "none"; return; end
    parts = strings(1, numel(ev));
    for i = 1:numel(ev), parts(i) = sprintf("%s@%g", lab(i), ev(i)); end
    s = strjoin(parts, " -> ");
end

function m = maxnan(x)
    x = x(isfinite(x)); if isempty(x), m = NaN; else, m = max(x); end
end

function v = getcol(T, name, n)
    if ismember(name, T.Properties.VariableNames)
        v = T.(name); if ~isnumeric(v), v = str2double(string(v)); end
        v = double(v(:));
    else
        v = nan(n, 1);
    end
end

function s = strcol(T, name, n)
    if ismember(name, T.Properties.VariableNames), s = string(T.(name)); s(ismissing(s)) = "";
    else, s = strings(n, 1); end
    s = s(:);
end

function v = nz(v)
    v = double(v); v(isnan(v)) = -1;
end

function k = expnum(e)
    k = nan(numel(e), 1);
    for i = 1:numel(e)
        t = regexp(char(e(i)), '\d+', 'match', 'once');
        if ~isempty(t), k(i) = str2double(t); end
    end
    k(isnan(k)) = -1;
end
