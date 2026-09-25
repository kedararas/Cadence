function [Tcam, names] = cadence_camera_medians(d, varargin)
%CADENCE_CAMERA_MEDIANS  Per-camera medians: one row per camera of a recording.
%
%   Tcam = cadence_camera_medians(metrics)
%   Tcam = cadence_camera_medians(metrics, 'Meta', struct('File', "x.mat", ...))
%   [Tcam, names] = cadence_camera_medians([])     % empty table + metric columns
%
%   cadence_recording_medians gives one value per metric per recording, read
%   from the FIRST voltage and FIRST calcium camera.  On a rig with several
%   voltage cameras (e.g. the four-camera pig wedge) the other cameras never
%   reach the summary.  This function runs the same median code once per
%   camera, so every camera gets its own row:
%     voltage camera k : cadence_recording_medians(d, 'VoltageCam', k, 'CalciumCam', 0)
%                        -> voltage metrics (APD, activation, CV, DF/RI/OI,
%                           AP alternans, Vm-only substrate, voltage QC)
%     calcium camera k : cadence_recording_medians(d, 'VoltageCam', 0, 'CalciumCam', k)
%                        -> calcium metrics (CaTD, decay, tau, Ca alternans, Ca QC)
%   Metrics that need a voltage-calcium pair (V-Ca delay, Ca-AP coupling,
%   combined substrate risk) or that are properties of the recording (beats,
%   pacing rate) stay in the recording-level table and are not repeated here.
%
%   Camera roles come from metrics.camera_roles (written by
%   cadence_extract_features).  Files extracted before roles were recorded
%   fall back to CAM1 = voltage, CAM2 = calcium.  Only cameras present in the
%   file are listed.
%
%   Output: a table with the 'Meta' fields first (in the order given), then
%   Camera ("CAM3"), Signal ("voltage" / "calcium"), then one column per
%   metric in the voltage and calcium families (NaN for the other family).
%   names is the list of those metric columns.

    p = inputParser;
    p.addParameter('Meta', struct(), @isstruct);
    p.addParameter('CVEdgeMargin', 8);
    p.addParameter('Polyfit', true);
    p.parse(varargin{:});
    o = p.Results;

    L = cadence_metric_labels();
    all_names = L(:,1);
    is_v  = metric_family(all_names, 'voltage');
    is_ca = metric_family(all_names, 'calcium');
    names = all_names(is_v | is_ca);

    Tcam = empty_table(o.Meta, names);
    if isempty(d), return; end
    if ischar(d) || isstring(d)
        S = load(char(d));
        if isfield(S, 'cmos_all_data'), d = S.cmos_all_data;
        else, fn = fieldnames(S); d = S.(fn{1}); end
        clear S
    end

    [vcams, cacams] = roles_of(d);
    jobs = [vcams(:), ones(numel(vcams), 1); cacams(:), 2 * ones(numel(cacams), 1)];
    rows = cell(size(jobs, 1), 1);
    for j = 1:size(jobs, 1)
        cam = jobs(j, 1);
        if jobs(j, 2) == 1
            vals = cadence_recording_medians(d, 'VoltageCam', cam, 'CalciumCam', 0, ...
                       'CVEdgeMargin', o.CVEdgeMargin, 'Polyfit', o.Polyfit);
            keep = is_v;  sig = "voltage";
        else
            vals = cadence_recording_medians(d, 'VoltageCam', 0, 'CalciumCam', cam, ...
                       'CVEdgeMargin', o.CVEdgeMargin, 'Polyfit', o.Polyfit);
            keep = is_ca; sig = "calcium";
        end
        r = o.Meta;
        r.Camera = string(sprintf('CAM%d', cam));
        r.Signal = sig;
        for k = 1:numel(all_names)
            if ~(is_v(k) || is_ca(k)), continue; end
            if keep(k), r.(all_names{k}) = double(vals.(all_names{k}));
            else,       r.(all_names{k}) = NaN; end
        end
        rows{j} = r;
    end
    if ~isempty(rows)
        Tcam = struct2table([rows{:}], 'AsArray', true);
        Tcam = coerce(Tcam, o.Meta, names);
    end
end


% =========================================================================
function [vc, cc] = roles_of(d)
    present = @(k) isfield(d, sprintf('CAM%d', k));
    if isfield(d, 'camera_roles') && isstruct(d.camera_roles)
        vc = double(d.camera_roles.voltage(:)');
        cc = double(d.camera_roles.calcium(:)');
    else
        vc = 1;  cc = 2;
    end
    vc = vc(arrayfun(present, vc));
    cc = cc(arrayfun(present, cc));
end

function tf = metric_family(n, fam)
% Voltage: *_v and v_* (per-camera QC) plus capture ratio.  Calcium: ca_*,
% catd*, *_ca.  Pair and recording-level metrics belong to neither.
    pair = {'vc_delay', 'ca_ap_coupling', 'ca_ap_inphase', 'risk_map', 'risk_global', ...
            'n_beats', 'pacing_hz', 'v_ca_identical'};
    n = n(:);
    switch fam
        case 'voltage'
            tf = endsWith(n, '_v') | startsWith(n, 'v_') | strcmp(n, 'capture_ratio');
        case 'calcium'
            tf = startsWith(n, 'ca_') | startsWith(n, 'catd') | endsWith(n, '_ca');
    end
    tf = tf & ~ismember(n, pair);
end

function T = empty_table(meta, names)
    r = meta;  r.Camera = "";  r.Signal = "";
    for k = 1:numel(names), r.(names{k}) = NaN; end
    T = coerce(struct2table(r, 'AsArray', true), meta, names);
    T = T([], :);
end

function T = coerce(T, meta, names)
% Strings stay strings (missing -> ""), metric columns are double.
    mf = fieldnames(meta);
    for k = 1:numel(mf)
        if ~(isnumeric(meta.(mf{k})) || islogical(meta.(mf{k})))
            c = string(T.(mf{k}));  c(ismissing(c)) = "";  T.(mf{k}) = c;
        end
    end
    T.Camera = string(T.Camera);  T.Signal = string(T.Signal);
    for k = 1:numel(names), T.(names{k}) = double(T.(names{k})); end
end
