function recs = validate_schema(d, stage, recording)
%VALIDATE_SCHEMA  Structural contract check for cmos_all_data at a boundary.
%
%   recs = validate_schema(cmos_all_data, stage)
%   recs = validate_schema(cmos_all_data, stage, recording)
%
%   Runs at LOAD time (before any stage consumes the struct).  It checks only
%   STRUCTURE — required fields present and of the right type — not values
%   (that is validate_conditioned / validate_metrics).  Catches the
%   field-name / type mismatches that otherwise surface as cryptic
%   "Unrecognized field" errors or silent NaNs deep inside a stage, e.g.:
%     * 'frame_rate' expected but only 'frequency'/'acqFreq' present
%     * 'df_map' present but a cell instead of a numeric matrix
%     * 'pixel_size' absent (thresholds silently fall back to geometry)
%
%   stage      which contract to enforce: "raw" | "conditioned" | "metrics".
%   recording  optional identifier for the records (default "").
%
%   Returns qc_check records with stage label "schema".  A FAIL means a
%   downstream stage will error or misbehave — fail fast at the boundary.

    if nargin < 3, recording = ""; end
    name = string(recording);

    recs = repmat(qc_check("","","" ,"PASS",NaN,""), 0, 1);
    function add(check, status, message)
        recs(end+1, 1) = qc_check(name, "schema", check, status, NaN, message);
    end

    specs = schema_for(string(stage));
    if isempty(specs)
        add("stage", "FAIL", sprintf("unknown schema stage '%s'", stage));
        return;
    end

    for s = 1:numel(specs)
        sp = specs(s);
        present = isfield(d, sp.name);

        if present
            if check_type(d.(sp.name), sp.type)
                add(sp.name, "PASS", "");
            else
                add(sp.name, "FAIL", sprintf("'%s' must be %s, got %s", ...
                    sp.name, sp.type, class(d.(sp.name))));
            end
            continue;
        end

        % Field absent — is a likely-misnamed alias present?
        alias_hit = sp.aliases(arrayfun(@(a) isfield(d, a), sp.aliases));
        if ~isempty(alias_hit)
            add(sp.name, "FAIL", sprintf( ...
                "required '%s' missing; found alias '%s' — rename at the producing stage", ...
                sp.name, alias_hit(1)));
        elseif sp.required
            add(sp.name, "FAIL", sprintf("required field '%s' missing", sp.name));
        else
            add(sp.name, "WARN", sprintf("optional field '%s' absent", sp.name));
        end
    end
end


% ======================================================================
function specs = schema_for(stage)
%SCHEMA_FOR  Declarative field contract per pipeline stage.
%   spec(name, type, required, aliases)

    raw = [ ...
        spec("CAM1",     "numeric", true,  []), ...
        spec("acqFreq",  "numeric", true,  ["frequency","frame_rate","Fs"]), ...
        spec("analog1",  "numeric", false, []) ];

    conditioned = [ raw, ...
        spec("CAM1_SNR",     "numeric", true,  []), ...
        spec("CAM1_average", "numeric", false, []) ];

    % The 'metrics' struct is what Feature Extraction saves via compile_data
    % (the *-metrics.mat files) and what Signal Analysis / Conduction Velocity
    % consume.  It carries the conditioned stacks plus the extracted EP
    % features under ep_metrics — NOT df_map / ps_dynamics (those belong to the
    % arrhythmia-dynamics path).  Keep this list in sync with compile_data.
    metrics = [ ...
        spec("CAM1",       "numeric", true,  []), ...
        spec("acqFreq",    "numeric", true,  ["frequency","frame_rate","Fs"]), ...
        spec("ep_metrics", "struct",  true,  ["metrics"]), ...
        spec("num_files",  "numeric", true,  []), ...
        spec("window",     "numeric", false, []), ...
        spec("data_masks", "cell",    false, []), ...
        spec("analog1",    "numeric", false, []) ];

    switch stage
        case "raw",         specs = raw;
        case "conditioned", specs = conditioned;
        case "metrics",     specs = metrics;
        % Unknown stage -> empty struct array of the spec shape (the caller
        % checks isempty(specs) and reports an "unknown stage" FAIL). NB: spec
        % is a local function, not a class, so spec.empty is invalid.
        otherwise,          specs = struct('name', {}, 'type', {}, ...
                                            'required', {}, 'aliases', {});
    end
end


function s = spec(name, type, required, aliases)
    if isempty(aliases), aliases = strings(1,0); end
    s = struct('name', string(name), 'type', string(type), ...
               'required', logical(required), 'aliases', string(aliases));
end


function tf = check_type(v, t)
    switch t
        case "numeric", tf = isnumeric(v);
        case "cell",    tf = iscell(v);
        case "char",    tf = ischar(v) || isstring(v);
        case "logical", tf = islogical(v);
        case "struct",  tf = isstruct(v);
        case "object",  tf = isobject(v) || isstruct(v);   % handle classes or struct
        otherwise,      tf = true;                          % "any"
    end
end
