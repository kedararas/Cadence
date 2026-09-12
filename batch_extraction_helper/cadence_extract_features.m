function [d, status] = cadence_extract_features(d, opts)
%CADENCE_EXTRACT_FEATURES  Headless port of the Feature Extraction app's EXTRACT FEATURES.
%
%   [metrics, status] = cadence_extract_features(cmos_all_data)
%   [metrics, status] = cadence_extract_features(cmos_all_data, opts)
%
%   Runs, without the UI, the same per-recording pipeline that
%   Cadence_Feature_Extraction.mlapp runs when EXTRACT FEATURES / BATCH
%   PROCESSING is pressed, and returns the 'metrics' struct that the app's
%   SAVE DATA would write (cmos_all_data + ep_metrics + window/data_masks/FOV).
%   Signal Analysis and Conduction Velocity load the result unchanged.
%
%   The numerical work is done by the SAME helper functions the app calls
%   (compute_lat_50, conduction_velocity, analyzeAPAlternans,
%   analyzeCaTransientAlternans, assess_arrhythmia_substrate,
%   cardiacSpectralMetrics, create_snr_mask, extract_beat_frames,
%   normalize_data).  The app-private glue (extract_apd, extract_rise_time,
%   extract_tau, despeckle_act_map, extract_combo_masks, extract_v_c_delay)
%   is reproduced here line for line, with the app object replaced by a
%   parameter struct.  If those app methods change, this file must follow.
%
%   Defaults reproduce the app's batch behaviour on a dual-camera rig:
%     * no signal window drawn  -> every metric except alternans and
%       DF/RI/OI is measured on the ensemble-averaged beat (CAM<n>_average)
%     * no user ROI / SNR toggle -> per-camera ADAPTIVE SNR mask
%       (create_snr_mask(snr, [], true))
%     * CAM1 = voltage (LAT, APD80, repolarization 80%, rise 20-90%,
%       local CV at FOV mm), CAM2 = calcium (V-Ca delay, CaTD80, decay 80%,
%       rise 20-90%, tau 90%->20%), AP alternans on CAM1, Ca alternans on
%       CAM2 with the combined arrhythmia substrate, DF/RI/OI on CAM1.
%     * representative pixel = image centre (what the app uses when the
%       user has not right-clicked).
%
%   opts fields (all optional)
%     FOV_mm          field of view, mm (default 20 -- the value stored in the
%                     lab's existing *-metrics files)
%     VoltageCams     cameras treated as voltage   (default 1)
%     CalciumCams     cameras treated as calcium   (default 2, if present)
%     APD_pct         APD level, %                 (default 80)
%     Rep_pct         repolarization level, %      (default 80)
%     APRise          [lo hi] AP rise levels, %    (default [20 90])
%     CaTD_pct        Ca transient duration, %     (default 80)
%     CaDecay_pct     Ca decay time level, %       (default 80)
%     CaRise          [lo hi] Ca rise levels, %    (default [20 90])
%     TauStart/TauEnd tau fit window as amplitude fractions (default 0.9 / 0.2,
%                     the app's dropdown defaults)
%     DoActTime, DoAPD, DoRep, DoAPRise, DoCV, DoVCDelay, DoCaTD, DoCaDecay,
%     DoCaRise, DoTau, DoAlternans, DoComplexity   logical switches (all true)
%     ComplexityCams  cameras for DF/RI/OI          (default = VoltageCams)
%     ExtraAPDLevels  additional duration levels, %, stored as ep_metrics
%                     fields apd<L>_data (voltage cams) and ca<L>_data (calcium
%                     cams) in the standard 6-slot layout (default 50, i.e.
%                     APD50 and CaTD50).  These fields are an ADDITION to what
%                     the app writes; pass [] to omit them.
%     Log             function handle for messages (default @(s) fprintf('%s\n',s))
%
%   status
%     .features   cellstr of ep_metrics fields produced
%     .errors     cellstr of stage errors (empty = clean run)
%     .timing     struct of seconds per stage
%     .settings   the resolved opts
%
%   See also cadence_batch_extract, cadence_recording_medians.

    if nargin < 2 || isempty(opts), opts = struct(); end
    P = resolve_options(d, opts);
    logf = P.Log;

    status = struct('features', {{}}, 'errors', {{}}, 'timing', struct(), 'settings', P);

    % The app starts each recording with an empty ep_metrics and the
    % window/mask state cleared; window stays empty (ensemble-average path).
    d.ep_metrics = struct();
    P.num_files = double(d.num_files);
    P.acqFreq   = double(d.acqFreq);

    % representative pixel: image centre (get_pixel_location with app.pixel empty)
    [r, c] = size(d.CAM1(:,:,1));
    P.row = round(r/2);  P.col = round(c/2);

    % ---- 1B) combo masks (adaptive SNR mask per camera) -------------------
    t0 = tic;
    P.combo_masks = cell(4,1);
    for i = 1:P.num_files
        snr = get_field(d, sprintf('CAM%d_SNR', i));
        if isempty(snr)
            logf(sprintf('CAM%d: no SNR map, no mask applied.', i));
            continue;
        end
        [snr_mask, snr_info] = create_snr_mask(snr, [], true);
        snr_mask = double(snr_mask);
        snr_mask(~snr_mask) = NaN;
        P.combo_masks{i,1} = snr_mask;
        logf(sprintf('CAM%d: adaptive SNR mask, threshold %.2f (%.0f%% tissue).', ...
            i, snr_info.threshold, 100*snr_info.tissue_fraction));
    end
    status.timing.masks = toc(t0);

    % ---- stages, in the app's order ---------------------------------------
    stages = { ...
        'act_times',     P.DoActTime,    @() stage_act_time(d, P); ...
        'apd_data',      P.DoAPD,        @() stage_v_apd(d, P); ...
        'rep_data',      P.DoRep,        @() stage_v_rep(d, P); ...
        'ap_rise_times', P.DoAPRise,     @() stage_v_rise(d, P); ...
        'local_cv',      P.DoCV,         @() stage_cv(d, P); ...
        'vc_delay',      P.DoVCDelay,    @() stage_v_c_delay(d, P); ...
        'ca_data',       P.DoCaTD,       @() stage_c_apd(d, P); ...
        'ca_rep_data',   P.DoCaDecay,    @() stage_c_rep(d, P); ...
        'ca_rise_times', P.DoCaRise,     @() stage_c_rise(d, P); ...
        'ca_tau',        P.DoTau,        @() stage_tau(d, P); ...
        'extra_apd',     ~isempty(P.ExtraAPDLevels), @() stage_extra_apd(d, P); ...
        'complexity_data', P.DoComplexity, @() stage_complexity(d, P); ...
        'alternans_data',  P.DoAlternans,  @() stage_alternans(d, P)};

    for s = 1:size(stages, 1)
        name = stages{s,1};
        if ~stages{s,2}, continue; end
        t0 = tic;
        try
            em = stages{s,3}();               % struct of ep_metrics fields to merge
            fn = fieldnames(em);
            for k = 1:numel(fn)
                if isempty(em.(fn{k})), continue; end     % stage produced nothing for this field
                d.ep_metrics.(fn{k}) = merge_cells(get_field(d.ep_metrics, fn{k}), em.(fn{k}));
                if ~any(strcmp(status.features, fn{k}))
                    status.features{end+1} = fn{k}; %#ok<AGROW>
                end
            end
            logf(sprintf('%-16s done in %.1f s', name, toc(t0)));
        catch ME
            msg = sprintf('%s FAILED: %s', name, ME.message);
            logf(msg);
            status.errors{end+1} = msg; %#ok<AGROW>
        end
        status.timing.(name) = toc(t0);
    end

    % ---- compile_data: what SAVE DATA writes --------------------------------
    d.window     = [];
    d.data_masks = cell(4,2);
    d.FOV        = P.FOV_mm;
end


% =========================================================================
%                              OPTIONS
% =========================================================================
function P = resolve_options(d, opts)
    def = struct( ...
        'FOV_mm', 20, ...
        'VoltageCams', 1, ...
        'CalciumCams', [], ...
        'APD_pct', 80, 'Rep_pct', 80, 'APRise', [20 90], ...
        'CaTD_pct', 80, 'CaDecay_pct', 80, 'CaRise', [20 90], ...
        'TauStart', 0.9, 'TauEnd', 0.2, ...
        'DoActTime', true, 'DoAPD', true, 'DoRep', true, 'DoAPRise', true, ...
        'DoCV', true, 'DoVCDelay', true, 'DoCaTD', true, 'DoCaDecay', true, ...
        'DoCaRise', true, 'DoTau', true, 'DoAlternans', true, 'DoComplexity', true, ...
        'ComplexityCams', [], ...
        'ExtraAPDLevels', 50, ...
        'Log', @(s) fprintf('%s\n', s));
    P = def;
    fn = fieldnames(opts);
    for k = 1:numel(fn)
        if ~isfield(def, fn{k})
            error('cadence_extract_features:badOption', 'Unknown option ''%s''.', fn{k});
        end
        P.(fn{k}) = opts.(fn{k});
    end
    nf = double(d.num_files);
    if isempty(P.CalciumCams) && nf >= 2 && ~isfield(opts, 'CalciumCams')
        P.CalciumCams = 2;
    end
    P.VoltageCams = P.VoltageCams(P.VoltageCams <= nf);
    P.CalciumCams = P.CalciumCams(P.CalciumCams <= nf);

    % num_files is a header field and can overstate what the file holds: a
    % voltage-only recording saved with num_files = 2 has no CAM2 array.  Keep
    % only cameras that are actually there, or every stage touching the missing
    % camera throws -- and stage_alternans then loses the voltage alternans it
    % had already computed along with the failed calcium branch.
    present = @(i) isfield(d, sprintf('CAM%d', i)) && ~isempty(d.(sprintf('CAM%d', i)));
    dropped = [P.VoltageCams(~arrayfun(present, P.VoltageCams)), ...
               P.CalciumCams(~arrayfun(present, P.CalciumCams))];
    P.VoltageCams = P.VoltageCams(arrayfun(present, P.VoltageCams));
    P.CalciumCams = P.CalciumCams(arrayfun(present, P.CalciumCams));
    if ~isempty(dropped)
        P.Log(sprintf('CAM%s declared by num_files but not in the file; those stages are skipped.', ...
                      strjoin(arrayfun(@(i) num2str(i), dropped, 'UniformOutput', false), ', CAM')));
    end

    if isempty(P.ComplexityCams), P.ComplexityCams = P.VoltageCams; end
    P.ComplexityCams = P.ComplexityCams(arrayfun(present, P.ComplexityCams));
end


function v = get_field(s, f)
    if isstruct(s) && isfield(s, f), v = s.(f); else, v = []; end
end


function out = merge_cells(existing, new)
% Merge per-camera rows of a {camera x slot} cell into the existing field.
    if isempty(existing), out = new; return; end
    out = existing;
    for i = 1:size(new, 1)
        for j = 1:size(new, 2)
            if i <= size(out,1) && j <= size(out,2) && isempty(new{i,j}) && ~isempty(out{i,j})
                continue;   % keep what an earlier stage wrote for this slot
            end
            out{i,j} = new{i,j};
        end
    end
end


function rep = pct_to_threshold(val)
    rep = max(0, 1 - val / 100);
end


function [cmos_data, avg_cmos_data] = get_cmos_data(d, cam)
    cmos_data = d.(sprintf('CAM%d', cam));
    avg_cmos_data = get_field(d, sprintf('CAM%d_average', cam));
end


function win = analysis_window(d, cam)
% ensemble || isempty(oap_window) -> the averaged beat (no window is ever drawn here)
    [~, win] = get_cmos_data(d, cam);
    if isempty(win)
        error('cadence_extract_features:noAverage', ...
              'CAM%d has no ensemble average (CAM%d_average); run Signal Conditioning with ensemble averaging.', cam, cam);
    end
end


function bg = get_bg_image(d, cam)
    bg = d.(sprintf('CAM%d_image', cam));
end


% =========================================================================
%                        ACTIVATION TIME  (extract_act_time)
% =========================================================================
function em = stage_act_time(d, P)
    em = struct(); em.act_times = cell(0, 6);
    for i = P.VoltageCams
        win_data = analysis_window(d, i);
        act_times = compute_lat_50(win_data);
        act_times = act_times * (1000 / P.acqFreq);            % frames -> ms
        mask_act_times = act_times;
        if ~isempty(P.combo_masks{i,1})
            mask_act_times = act_times .* P.combo_masks{i,1};
        end
        [mask_act_times, n_spk] = despeckle_act_map(mask_act_times);
        if n_spk > 0
            P.Log(sprintf('CAM%d activation map: despeckled %d isolated early/late pixel(s).', i, n_spk));
        end
        em.act_times{i,1} = mask_act_times;
        em.act_times{i,2} = get_bg_image(d, i);
        em.act_times{i,3} = squeeze(normalize_data(win_data(P.row, P.col, :)));
        em.act_times{i,4} = mask_act_times(P.row, P.col);
        em.act_times{i,5} = num2str(i);
        em.act_times{i,6} = act_times;
    end
end


function [clean_map, n_speckle] = despeckle_act_map(act_map)
% Verbatim from the app: isolated pixels > 5 ms from their 3x3 median are
% replaced by that median.
    speckle_ms = 5;
    [nr, nc] = size(act_map);
    a = double(act_map);  a(a == 0) = NaN;
    clean_map = act_map;
    [rr, cc] = find(isfinite(a));
    n_speckle = 0;
    for k = 1:numel(rr)
        r = rr(k);  c = cc(k);
        r0 = max(1,r-1); r1 = min(nr,r+1);
        c0 = max(1,c-1); c1 = min(nc,c+1);
        w  = a(r0:r1, c0:c1);
        mw = median(w(isfinite(w)), 'omitnan');
        if isfinite(mw) && abs(a(r,c) - mw) > speckle_ms
            n_speckle = n_speckle + 1;
            clean_map(r,c) = mw;
        end
    end
end


% =========================================================================
%                     APD / REPOLARIZATION  (extract_apd)
% =========================================================================
function em = stage_v_apd(d, P)
    em = struct(); em.apd_data = cell(0,6); em.rep_data = cell(0,6);
    for i = P.VoltageCams
        a = cadence_extract_apd(d, P, i, pct_to_threshold(P.APD_pct), sprintf('APD%d', P.APD_pct));
        em.apd_data(i,1:6) = {a{1,1}, a{1,2}, a{1,3}, a{1,4}, a{1,5}, a{1,6}};
        em.rep_data(i,1:6) = {a{1,7}, a{1,2}, a{1,3}, a{1,8}, a{1,5}, a{1,9}};
    end
end


function em = stage_v_rep(d, P)
    em = struct(); em.rep_data = cell(0,6);
    if P.DoAPD && P.Rep_pct == P.APD_pct
        return;   % identical to what stage_v_apd already wrote; the app recomputes it
    end
    for i = P.VoltageCams
        a = cadence_extract_apd(d, P, i, pct_to_threshold(P.Rep_pct), sprintf('Rep%d', P.Rep_pct));
        em.rep_data(i,1:6) = {a{1,7}, a{1,2}, a{1,3}, a{1,8}, a{1,5}, a{1,9}};
    end
end


function em = stage_c_apd(d, P)
    em = struct(); em.ca_data = cell(0,6); em.ca_rep_data = cell(0,6);
    for i = P.CalciumCams
        a = cadence_extract_apd(d, P, i, pct_to_threshold(P.CaTD_pct), sprintf('CaTD%d', P.CaTD_pct));
        em.ca_data(i,1:6)     = {a{1,1}, a{1,2}, a{1,3}, a{1,4}, a{1,5}, a{1,6}};
        em.ca_rep_data(i,1:6) = {a{1,7}, a{1,2}, a{1,3}, a{1,8}, a{1,5}, a{1,9}};
    end
end


function em = stage_c_rep(d, P)
    em = struct(); em.ca_rep_data = cell(0,6);
    if P.DoCaTD && P.CaDecay_pct == P.CaTD_pct
        return;
    end
    for i = P.CalciumCams
        a = cadence_extract_apd(d, P, i, pct_to_threshold(P.CaDecay_pct), sprintf('CaDecay%d', P.CaDecay_pct));
        em.ca_rep_data(i,1:6) = {a{1,7}, a{1,2}, a{1,3}, a{1,8}, a{1,5}, a{1,9}};
    end
end


function em = stage_extra_apd(d, P)
% Additional duration levels (e.g. APD50 / CaTD50), same engine, extra fields.
    em = struct();
    for L = P.ExtraAPDLevels(:)'
        fv = sprintf('apd%d_data', round(L));  fc = sprintf('ca%d_data', round(L));
        em.(fv) = cell(0,6);  em.(fc) = cell(0,6);
        for i = P.VoltageCams
            a = cadence_extract_apd(d, P, i, pct_to_threshold(L), sprintf('APD%d', round(L)));
            em.(fv)(i,1:6) = {a{1,1}, a{1,2}, a{1,3}, a{1,4}, a{1,5}, a{1,6}};
        end
        for i = P.CalciumCams
            a = cadence_extract_apd(d, P, i, pct_to_threshold(L), sprintf('CaTD%d', round(L)));
            em.(fc)(i,1:6) = {a{1,1}, a{1,2}, a{1,3}, a{1,4}, a{1,5}, a{1,6}};
        end
    end
end


% =========================================================================
%                         RISE TIME  (extract_rise_time)
% =========================================================================
function em = stage_v_rise(d, P)
    em = struct(); em.ap_rise_times = cell(0,6);
    for i = P.VoltageCams
        r = extract_rise_time(d, P, i, P.APRise(1)/100, P.APRise(2)/100);
        em.ap_rise_times(i,1:6) = r;
    end
end


function em = stage_c_rise(d, P)
    em = struct(); em.ca_rise_times = cell(0,6);
    for i = P.CalciumCams
        r = extract_rise_time(d, P, i, P.CaRise(1)/100, P.CaRise(2)/100);
        em.ca_rise_times(i,1:6) = r;
    end
end


function data = extract_rise_time(d, P, camera, rise_min, rise_max)
% Port of Cadence_Feature_Extraction/extract_rise_time (ensemble path).
    row = P.row; col = P.col;
    win_data = analysis_window(d, camera);
    [num_rows, num_cols, T] = size(win_data);

    Y  = double(reshape(win_data, num_rows * num_cols, T));
    N  = size(Y, 1);
    tt = 1:T;

    b0  = prctile(Y,  5, 2);
    p0  = prctile(Y, 99, 2);
    mid = b0 + 0.5 * (p0 - b0);

    hi_mid    = Y >= mid;
    sustained = hi_mid(:, 1:end-1) & hi_mid(:, 2:end);
    [has_anchor, anchor] = max(sustained, [], 2);
    anchor(~has_anchor)  = T;

    pre  = Y;  pre(tt >= anchor)  = NaN;
    base = median(pre, 2, 'omitnan');
    post = Y;  post(tt <  anchor) = NaN;
    pk   = max(post, [], 2);
    amp  = pk - base;

    lo_thr = base + rise_min * amp;
    hi_thr = base + rise_max * amp;

    [has_end,   i_end] = max((Y >= hi_thr) & (tt >= anchor), [], 2);
    [has_start, j_rev] = max(fliplr((Y <= lo_thr) & (tt <= anchor)), [], 2);
    i_start = T - j_rev + 1;

    idx = (1:N)';

    ya    = Y(sub2ind([N T], idx, max(i_end - 1, 1)));
    yb    = Y(sub2ind([N T], idx, i_end));
    den   = yb - ya;
    f_end = (hi_thr - ya) ./ den;
    f_end(den <= 0) = 0;
    t_end = (i_end - 1) + min(max(f_end, 0), 1);

    yc      = Y(sub2ind([N T], idx, i_start));
    yd      = Y(sub2ind([N T], idx, min(i_start + 1, T)));
    den2    = yd - yc;
    f_start = (lo_thr - yc) ./ den2;
    f_start(den2 <= 0) = 0;
    t_start = i_start + min(max(f_start, 0), 1);

    rt = t_end - t_start;
    ok = has_anchor & has_end & has_start & isfinite(base) & ...
         (amp > 0) & (i_end > 1) & (i_start < i_end) & (rt > 0);

    rise_frames     = nan(N, 1);
    rise_frames(ok) = rt(ok);
    rise_time       = reshape(rise_frames, num_rows, num_cols) * (1000 / P.acqFreq);

    p = sub2ind([num_rows num_cols], row, col);
    if isfinite(base(p)) && amp(p) > 0
        oap = (Y(p, :)' - base(p)) / amp(p);
    else
        oap = nan(T, 1);
    end
    rise_plot = nan(T, 1);
    if ok(p)
        rs = max(1, floor(t_start(p)));
        re = min(T, ceil(t_end(p)));
        rise_plot(rs:re) = oap(rs:re);
    end

    data = cell(1, 6);
    data{1,1} = rise_time;
    if ~isempty(P.combo_masks{camera,1})
        data{1,1} = rise_time .* P.combo_masks{camera,1};
    end
    data{1,2} = get_bg_image(d, camera);
    data{1,3} = oap;
    data{1,4} = rise_plot;
    data{1,5} = num2str(camera);
    data{1,6} = rise_time;
end


% =========================================================================
%                    LOCAL CONDUCTION VELOCITY  (extract_cv)
% =========================================================================
function em = stage_cv(d, P)
    em = struct(); em.local_cv = cell(0,4);
    for i = P.VoltageCams
        win_data = analysis_window(d, i);
        act_times = compute_lat_50(win_data);
        mask_act_times = act_times;
        if ~isempty(P.combo_masks{i,1})
            mask_act_times = act_times .* P.combo_masks{i,1};
        end
        num_pixels = size(act_times, 1);
        pixel_resolution = P.FOV_mm / num_pixels;                 % mm
        mask_act_times = mask_act_times * (1000 / P.acqFreq);     % frames -> ms
        mask_act_times = despeckle_act_map(mask_act_times);
        [Vmag, Vx, Vy, Gx, Gy, meta] = conduction_velocity(mask_act_times, pixel_resolution, pixel_resolution, []);
        Vmag(meta.smooth_support < 0.75) = NaN;

        cv_data = cell(1, 7);
        cv_data{1,1} = Vmag*100;   % cm/s
        cv_data{1,2} = Vx;  cv_data{1,3} = Vy;
        cv_data{1,4} = Gx;  cv_data{1,5} = Gy;
        cv_data{1,6} = meta;
        cv_data{1,7} = mask_act_times;
        em.local_cv{i,1} = cv_data;
        em.local_cv{i,2} = get_bg_image(d, i);
        em.local_cv{i,3} = mask_act_times;
        em.local_cv{i,4} = num2str(i);
    end
end


% =========================================================================
%                  VOLTAGE-CALCIUM DELAY  (extract_v_c_delay)
% =========================================================================
function em = stage_v_c_delay(d, P)
    em = struct(); em.vc_delay = cell(0,6);
    % Gate on the cameras actually resolved, not on num_files: that header field
    % can claim a second camera the file does not contain.
    if isempty(P.VoltageCams) || isempty(P.CalciumCams)
        P.Log(sprintf(['Voltage-Calcium Delay skipped: needs a voltage and a calcium camera ' ...
            '(this recording has %d declared camera(s); voltage %s, calcium %s).'], ...
            P.num_files, mat2str(P.VoltageCams), mat2str(P.CalciumCams)));
        return;
    end
    v_camera = P.VoltageCams(1);  ca_camera = P.CalciumCams(1);
    if ~isfield(d, 'analog1') || isempty(d.analog1)
        P.Log(['Voltage-Calcium Delay skipped: ensemble averaging with no analog1 pacing trace. ' ...
               'Each camera was aligned to its own upstroke, which removes the delay being measured.']);
        return;
    end
    v_win_data  = analysis_window(d, v_camera);
    ca_win_data = analysis_window(d, ca_camera);

    v_act_times  = compute_lat_50(v_win_data)  * (1000 / P.acqFreq);
    ca_act_times = compute_lat_50(ca_win_data) * (1000 / P.acqFreq);

    mask_v_act_times  = v_act_times;
    mask_ca_act_times = ca_act_times;
    if ~isempty(P.combo_masks{v_camera,1})
        mask_v_act_times = v_act_times .* P.combo_masks{v_camera,1};
    end
    if ~isempty(P.combo_masks{ca_camera,1})
        mask_ca_act_times = ca_act_times .* P.combo_masks{ca_camera,1};
    end

    v_c_oap = cell(2,1);
    v_c_oap{1,1} = squeeze(normalize_data(v_win_data(P.row, P.col, :)));
    v_c_oap{2,1} = squeeze(normalize_data(ca_win_data(P.row, P.col, :)));

    em.vc_delay{ca_camera,1} = mask_ca_act_times - mask_v_act_times;   % SIGNED: Ca minus Vm
    em.vc_delay{ca_camera,2} = get_bg_image(d, ca_camera);
    em.vc_delay{ca_camera,3} = v_c_oap;
    em.vc_delay{ca_camera,4} = mask_ca_act_times(P.row, P.col);
    em.vc_delay{ca_camera,5} = 'VCDelay';
    em.vc_delay{ca_camera,6} = ca_act_times - v_act_times;

    dl = em.vc_delay{ca_camera,1};
    dl = dl(isfinite(dl));
    if ~isempty(dl) && median(dl) < 0
        P.Log(sprintf(['WARNING: median V-Ca delay is %.1f ms (negative): calcium appears to LEAD ' ...
            'voltage. Check that CAM1 is voltage and CAM2 calcium.'], median(dl)));
    end
end


% =========================================================================
%                     CALCIUM DECAY CONSTANT  (extract_tau)
% =========================================================================
function em = stage_tau(d, P)
    em = struct(); em.ca_tau = cell(0,6);
    acqFreq = P.acqFreq;
    dt      = 1 / acqFreq;
    row = P.row; col = P.col;
    tau_start = P.TauStart;  tau_end = P.TauEnd;

    for i = P.CalciumCams
        win_data = analysis_window(d, i);
        [num_rows, num_cols, num_frames] = size(win_data);

        data_mask = ones(num_rows, num_cols);
        if ~isempty(P.combo_masks{i,1})
            data_mask = data_mask .* P.combo_masks{i,1};
        end

        Y  = double(reshape(win_data, num_rows * num_cols, num_frames));
        N  = size(Y, 1);
        T  = num_frames;
        tt = 1:T;
        idx = (1:N)';

        ws = max(3, round(0.005 * acqFreq));
        wc = max(5, 2*floor(0.015 * acqFreq / 2) + 1);
        Ys = movmedian(Y, ws, 2);
        Yc = movmedian(Y, wc, 2);

        b0  = prctile(Ys,  5, 2);
        p0  = prctile(Ys, 99, 2);
        rg  = p0 - b0;
        mid = b0 + 0.5 * rg;
        low = b0 + 0.1 * rg;

        hi_mid    = Ys >= mid;
        sustained = hi_mid(:, 1:end-1) & hi_mid(:, 2:end);
        [has_anchor, anchor] = max(sustained, [], 2);
        anchor(~has_anchor)  = T;

        cand = (Ys <= low) & (tt < anchor);
        [has_foot, j_rev] = max(fliplr(cand), [], 2);
        foot = T - j_rev + 1;
        foot(~has_foot) = max(anchor(~has_foot) - 1, 1);

        pre     = Y;  pre(tt > foot) = NaN;
        base    = median(pre, 2, 'omitnan');
        n_pre   = sum(~isnan(pre), 2);
        min_pre = max(3, round(0.005 * acqFreq));
        base(n_pre < min_pre) = NaN;
        clear pre

        sigma = 1.4826 * median(abs(diff(Y, 1, 2)), 2, 'omitnan') / sqrt(2);

        post         = Ys;  post(tt < anchor) = NaN;
        [pk_val, pk] = max(post, [], 2);
        amp          = pk_val - base;
        pk_amp       = Y(sub2ind([N T], idx, pk)) - base;
        clear post Ys

        sigma_f  = sigma * sqrt(pi / (2 * ws));
        amp_gate = (sqrt(2 * log(max(T / ws, 3))) + 3) * sigma_f;
        ok_amp   = isfinite(base) & isfinite(sigma) & (sigma > 0) & ...
                   (amp > 0) & (amp > amp_gate) & (pk_amp > 0) & ...
                   ~isnan(data_mask(:));

        amp_safe = amp;     amp_safe(~ok_amp) = NaN;
        pk_safe  = pk_amp;  pk_safe(~ok_amp)  = NaN;
        Yn       = (Y - base) ./ amp_safe;

        % slot 6: model-free 1/e
        lev  = 1 / exp(1);
        labs = base + lev .* pk_safe;
        [hit1, ix] = max(cumsum((Yc <= labs) & (tt > pk), 2) == 1, [], 2);
        ix   = max(ix, 2);
        y1   = Yc(sub2ind([N T], idx, ix - 1));
        y2   = Yc(sub2ind([N T], idx, ix));
        den1 = y1 - y2;
        frac = (y1 - labs) ./ den1;
        frac(~(den1 > 0)) = 0;
        frac = min(max(frac, 0), 1);
        t1e  = ((ix - 1 + frac) - pk) * dt * 1000;
        tau_mode_constant = nan(N, 1);
        g1e  = hit1 & ok_amp & isfinite(t1e) & (t1e > 0);
        tau_mode_constant(g1e) = t1e(g1e);
        clear Yc

        % slot 1: decay fit with a baseline term
        if tau_start >= 1
            rep_start = pk;
            has_start = ok_amp;
        else
            [has_start, rep_start] = max(cumsum((Yn <= tau_start) & (tt >= pk), 2) == 1, [], 2);
            rep_start(~has_start)  = T;
        end
        [has_end, rep_end] = max(cumsum((Yn <= tau_end) & (tt >= rep_start), 2) == 1, [], 2);
        rep_end(~has_end) = T;
        rep_end = max(rep_end, rep_start);

        v_s = Yn(sub2ind([N T], idx, rep_start));
        v_e = Yn(sub2ind([N T], idx, rep_end));
        ok_span = (v_s > 0) & (v_e <= 0.5 * v_s);

        D     = Y - base;
        inwin = (tt >= rep_start) & (tt <= rep_end);
        Trel  = ((tt - 1) - (rep_start - 1)) * dt;

        n_win = sum(inwin, 2);
        [hf, i_first] = max(inwin, [], 2);
        [~,  j2]      = max(fliplr(inwin), [], 2);
        span_ms = (T - j2 + 1 - i_first) * dt * 1000;
        min_pts     = max(4, round(0.005 * acqFreq));
        min_span_ms = 5;

        valid = inwin & (D > 0) & isfinite(D);
        Z = zeros(N, T);  Z(valid) = log(D(valid));
        w = zeros(N, T);  w(valid) = D(valid);
        S0  = sum(w, 2);
        S1  = sum(w .* Trel, 2);
        S2  = sum(w .* Trel.^2, 2);
        Sy  = sum(w .* Z, 2);
        Sty = sum(w .* Trel .* Z, 2);
        den = S0 .* S2 - S1.^2;
        den(~(den > 0)) = NaN;
        tau_seed_s = -1 ./ ((S0 .* Sty - S1 .* Sy) ./ den);
        clear Z w valid

        Dm = D;     Dm(~inwin) = 0;
        Tm = Trel;  Tm(~inwin) = 0;
        tau_s = tau_seed_s;
        tau_s(~isfinite(tau_s) | tau_s <= 0) = 20 * dt;
        tau_s = min(max(tau_s, dt), T * dt);
        for it = 1:8
            E    = exp(-Tm ./ tau_s) .* inwin;
            Aa   = sum(Dm .* E, 2) ./ max(sum(E.^2, 2), eps);
            Rr   = (Dm - Aa .* E) .* inwin;
            G    = Aa .* E .* Tm ./ tau_s.^2;
            step = sum(Rr .* G, 2) ./ max(sum(G.^2, 2), eps);
            step = min(max(step, -0.5 * tau_s), 0.5 * tau_s);
            tau_s = min(max(tau_s + step, dt), T * dt);
        end
        clear Dm Tm E Rr G

        ok_fit = ok_amp & ok_span & has_start & hf & ...
                 (n_win >= min_pts) & (span_ms >= min_span_ms) & ...
                 isfinite(tau_s) & (tau_s > dt) & (tau_s < T * dt);
        tau_constant = nan(N, 1);
        tau_constant(ok_fit) = tau_s(ok_fit) * 1000;

        p = sub2ind([num_rows num_cols], row, col);
        if ok_amp(p)
            oap = (Y(p, :)' - base(p)) / amp(p);
        else
            oap = nan(T, 1);
        end
        fitted_oap = nan(T, 1);
        if ok_fit(p)
            fitted_oap(rep_start(p):rep_end(p)) = oap(rep_start(p):rep_end(p));
        end
        clear Y Yn D inwin Trel

        tau_constant      = reshape(tau_constant,      num_rows, num_cols);
        tau_mode_constant = reshape(tau_mode_constant, num_rows, num_cols);

        em.ca_tau{i,1} = tau_constant;
        em.ca_tau{i,6} = tau_mode_constant;
        if ~isempty(P.combo_masks{i,1})
            em.ca_tau{i,1} = tau_constant      .* data_mask;
            em.ca_tau{i,6} = tau_mode_constant .* data_mask;
        end
        em.ca_tau{i,2} = get_bg_image(d, i);
        em.ca_tau{i,3} = oap;
        em.ca_tau{i,4} = fitted_oap;
        em.ca_tau{i,5} = num2str(i);

        elig  = ~isnan(data_mask);
        n_all = nnz(elig);
        n_bad = nnz(elig & isnan(em.ca_tau{i,1}));
        n_bd6 = nnz(elig & isnan(em.ca_tau{i,6}));
        P.Log(sprintf('Tau CAM%d: %d of %d masked pixels (%.1f%%) censored in the fit, %.1f%% in the 1/e estimate.', ...
            i, n_bad, n_all, 100 * n_bad / max(n_all,1), 100 * n_bd6 / max(n_all,1)));
    end
end


% =========================================================================
%                  DF / RI / OI  (extract_arr_complexity)
% =========================================================================
function em = stage_complexity(d, P)
    em = struct(); em.complexity_data = cell(0,5);
    for i = P.ComplexityCams
        if ~isempty(P.combo_masks{i,1})
            mask = P.combo_masks{i,1};
        else
            mask = extract_image_mask(get_bg_image(d, i), 50);
        end
        cmos_data = d.(sprintf('CAM%d', i));            % full raw stack
        bg_idx   = find(isnan(mask));
        n_pixels = size(cmos_data, 1) * size(cmos_data, 2);
        bg_3d  = bg_idx + (0:size(cmos_data,3)-1) * n_pixels;
        cmos_data(bg_3d) = NaN;
        [DF, RI, OI, ~] = cardiacSpectralMetrics(cmos_data, P.acqFreq);
        if isempty(DF), continue; end
        em.complexity_data{i,1} = {DF, RI, OI};
        em.complexity_data{i,2} = get_bg_image(d, i);
        em.complexity_data{i,3} = squeeze(normalize_data(cmos_data(P.row, P.col, :)));
        em.complexity_data{i,4} = [DF(P.row, P.col), RI(P.row, P.col), OI(P.row, P.col)];
        em.complexity_data{i,5} = num2str(i);
    end
end


% =========================================================================
%                         ALTERNANS  (extract_alternans)
% =========================================================================
function em = stage_alternans(d, P)
    em = struct(); em.alternans_data = cell(0,6);
    cams = [P.VoltageCams(:)', P.CalciumCams(:)'];
    if isempty(cams), return; end

    if isfield(d, 'analog1') && ~isempty(d.analog1)
        pacing = d.analog1;
        [sf, ef] = extract_beat_frames(pacing, 'FixedLength', true);
    else
        pacing = auto_detect_peaks(d.CAM1);
        [sf, ef] = extract_beat_frames(pacing, 'PreWindow', 10, 'FixedLength', true);
    end
    beat_frames = [sf(:), ef(:)];
    if size(beat_frames, 1) < 2
        error('cadence_extract_features:noBeats', 'fewer than 2 beats segmented; alternans needs beat pairs.');
    end
    P.Log(sprintf('Alternans: %d beats of %d frames.', size(beat_frames,1), beat_frames(1,2)-beat_frames(1,1)+1));

    voltage_alternans = 0;  calcium_alternans = 0;
    for i = P.VoltageCams
        a = extract_v_alternans(d, P, beat_frames, i);
        em.alternans_data(i,1:6) = a;
        voltage_alternans = i;
    end
    for i = P.CalciumCams
        a = extract_c_alternans(d, P, beat_frames, i);
        em.alternans_data(i,1:6) = a;
        calcium_alternans = i;
    end

    if voltage_alternans > 0 && calcium_alternans > 0
        v_alternans  = em.alternans_data{voltage_alternans,1};
        ca_alternans = em.alternans_data{calcium_alternans,1};
        substrate = assess_arrhythmia_substrate(v_alternans, 'CaResult', ca_alternans, ...
            'CaData', get_cmos_data(d, calcium_alternans), 'BeatFrames', beat_frames);
        em.alternans_data{calcium_alternans,5} = substrate;
    end
end


function data = extract_v_alternans(d, P, pacing, camera)
    mask = P.combo_masks{camera,1};
    cmos_data = get_cmos_data(d, camera);
    time = (0:size(cmos_data,3)-1) * (1000 / P.acqFreq);        % real ms
    alternans = analyzeAPAlternans(time, cmos_data, 'BeatFrames', pacing, 'Mask', mask);
    substrate = assess_arrhythmia_substrate(alternans, 'BeatFrames', pacing);
    data = cell(1, 6);
    data{1,1} = alternans;
    data{1,2} = get_bg_image(d, camera);
    data{1,3} = squeeze(normalize_data(cmos_data(P.row, P.col, :)));
    data{1,4} = num2str(camera);
    data{1,5} = substrate;
    data{1,6} = 'AP';
end


function data = extract_c_alternans(d, P, pacing, camera)
    mask = P.combo_masks{camera,1};
    cmos_data = get_cmos_data(d, camera);
    time = (0:size(cmos_data,3)-1) * (1000 / P.acqFreq);
    alternans = analyzeCaTransientAlternans(time, cmos_data, 'BeatFrames', pacing, 'Mask', mask);
    data = cell(1, 6);
    data{1,1} = alternans;
    data{1,2} = get_bg_image(d, camera);
    data{1,3} = squeeze(normalize_data(cmos_data(P.row, P.col, :)));
    data{1,4} = num2str(camera);
    data{1,5} = [];
    data{1,6} = 'Ca';
end
