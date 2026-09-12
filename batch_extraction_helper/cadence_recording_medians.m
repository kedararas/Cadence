function [vals, labels] = cadence_recording_medians(d, varargin)
%CADENCE_RECORDING_MEDIANS  Per-recording median of every summary metric.
%
%   vals = cadence_recording_medians(metrics)
%   vals = cadence_recording_medians('path/to/x-metrics.mat')
%   [vals, labels] = cadence_recording_medians(..., 'VoltageCam', 1, 'CalciumCam', 2)
%
%   metrics is the cmos_all_data struct with ep_metrics (what Feature
%   Extraction saves, or what cadence_extract_features returns).  Returns a
%   struct with one field per metric in cadence_metric_labels (NaN where a
%   metric is absent) and, optionally, the {name, label, desc} list.
%
%   Median = median over finite pixels of the masked slot-1 map.  Conduction
%   velocity is pooled over the interior of the tissue only (8 px erosion),
%   the same gate the Signal Analysis and Conduction Velocity modules use,
%   because one-sided smoothing biases CV high along the mask border.

    p = inputParser;
    p.addParameter('VoltageCam', 1);
    p.addParameter('CalciumCam', 2);
    p.addParameter('CVEdgeMargin', 8);
    p.addParameter('Polyfit', true);
    p.parse(varargin{:});
    o = p.Results;

    if ischar(d) || isstring(d)
        S = load(char(d));
        if isfield(S, 'cmos_all_data'), d = S.cmos_all_data;
        else, fn = fieldnames(S); d = S.(fn{1}); end
    end

    labels = cadence_metric_labels();
    vals = struct();
    for k = 1:size(labels, 1), vals.(labels{k,1}) = NaN; end

    if ~isfield(d, 'ep_metrics') || ~isstruct(d.ep_metrics), return; end
    em = d.ep_metrics;
    v  = o.VoltageCam;  ca = o.CalciumCam;

    % ---- voltage maps ------------------------------------------------------
    act = slot(em, 'act_times', v, 1);
    act(act == 0) = NaN;                       % legacy zero background
    vals.act_time_v = nanmed(act);
    apd = slot(em, 'apd_data', v, 1);
    vals.apd_v      = nanmed(apd);
    vals.rep_time_v = nanmed(slot(em, 'rep_data', v, 1));
    vals.apd50_v    = nanmed(level_map(d, em, 'apd50_data', v, 0.5, 'APD50'));
    vals.ap_rise_v  = nanmed(slot(em, 'ap_rise_times', v, 1));
    if ~isempty(act)
        vals.tissue_frac_v = nnz(isfinite(act)) / numel(act);
        if ~isempty(apd)
            vals.apd_valid_frac_v = nnz(isfinite(apd)) / max(nnz(isfinite(act)), 1);
        end
    end

    % ---- conduction velocity ----------------------------------------------
    cvc = slot(em, 'local_cv', v, 1);
    if iscell(cvc) && numel(cvc) >= 1 && isnumeric(cvc{1,1})
        vals.cv_v = nanmed(interior_only(cvc{1,1}, o.CVEdgeMargin));
        if o.Polyfit && numel(cvc) >= 7 && isnumeric(cvc{1,7}) && ~isempty(cvc{1,7})
            try
                act_map = cvc{1,7};
                fov = 20; if isfield(d, 'FOV') && ~isempty(d.FOV), fov = double(d.FOV); end
                px  = fov / size(act_map, 1);
                Vmag = conduction_velocity(act_map, px, px, struct('method', 'polyfit'));
                vals.cv_polyfit_v = nanmed(interior_only(Vmag * 100, o.CVEdgeMargin));
            catch
            end
        end
    end

    % ---- calcium maps -----------------------------------------------------
    vals.vc_delay    = nanmed(slot(em, 'vc_delay', ca, 1));
    catd             = slot(em, 'ca_data', ca, 1);
    vals.catd_ca     = nanmed(catd);
    vals.ca_decay_ca = nanmed(slot(em, 'ca_rep_data', ca, 1));
    vals.catd50_ca   = nanmed(level_map(d, em, 'ca50_data', ca, 0.5, 'CaTD50'));
    vals.ca_rise_ca  = nanmed(slot(em, 'ca_rise_times', ca, 1));
    vals.ca_tau_ca   = nanmed(slot(em, 'ca_tau', ca, 1));
    vals.ca_tau_1e_ca = nanmed(slot(em, 'ca_tau', ca, 6));
    ca_mask_n = NaN;
    if isfield(d, sprintf('CAM%d_SNR', ca))
        try
            m = create_snr_mask(d.(sprintf('CAM%d_SNR', ca)), [], true);
            ca_mask_n = nnz(m);
            vals.tissue_frac_ca = ca_mask_n / numel(m);
        catch
        end
    end
    if ~isempty(catd) && isfinite(ca_mask_n) && ca_mask_n > 0
        vals.catd_valid_frac_ca = nnz(isfinite(catd)) / ca_mask_n;
    end

    % ---- QC: are the voltage and calcium results really from different cameras? ----
    vals.v_ca_identical = 0;
    act_ca = slot(em, 'act_times', ca, 1);
    if ~isempty(apd) && ~isempty(catd) && isequaln(apd, catd), vals.v_ca_identical = 1; end
    if ~isempty(act) && ~isempty(act_ca) && isequaln(act, act_ca), vals.v_ca_identical = 1; end
    if v ~= ca && isfield(d, sprintf('CAM%d', v)) && isfield(d, sprintf('CAM%d', ca))
        try
            a1 = d.(sprintf('CAM%d', v))(:,:,1);  a2 = d.(sprintf('CAM%d', ca))(:,:,1);
            if isequaln(a1, a2), vals.v_ca_identical = 1; end
        catch
        end
    end

    % ---- QC: polarity, judged against the pacing stimulus -----------------
    % Conditioning decides polarity per recording from signal shape, and that
    % test fails at fast rates (the transient fills the cycle, so the signal
    % spends most of its time elevated and the histogram criterion votes
    % "inverted" while the slope criterion collapses to zero).  The stimulus is
    % immune to that: a correctly oriented signal deflects UP after a capture.
    if isfield(d, 'analog1') && ~isempty(d.analog1) && isfield(d, 'acqFreq') && ...
            exist('check_polarity_paced', 'file') == 2
        for cc = unique([v ca])
            cf = sprintf('CAM%d', cc);  sf = sprintf('CAM%d_SNR', cc);
            if ~isfield(d, cf) || isempty(d.(cf)), continue; end
            try
                mk = [];
                if isfield(d, sf) && ~isempty(d.(sf)), mk = create_snr_mask(d.(sf), 3, true); end
                [~, cnf] = check_polarity_paced(d.(cf), d.analog1, d.acqFreq, mk);
                if cc == v,  vals.v_polarity_conf  = cnf; end
                if cc == ca, vals.ca_polarity_conf = cnf; end
            catch
            end
        end
    end

    % ---- QC: pacing rate from the stimulus channel ----------------------
    if isfield(d, 'analog1') && ~isempty(d.analog1) && isfield(d, 'acqFreq')
        try
            pac = double(d.analog1(:));
            % same detector as extract_beat_frames (pulse > half max), but with a
            % short refractory so fast pacing (50 ms) is not merged
            [~, onsets] = findpeaks(pac, 'MinPeakHeight', 0.5 * max(pac), ...
                                    'MinPeakDistance', max(3, round(0.02 * double(d.acqFreq))));
            if numel(onsets) >= 2
                vals.pacing_hz = double(d.acqFreq) / median(diff(onsets));
            end
        catch
        end
    end

    % ---- voltage spectral -------------------------------------------------
    cx = slot(em, 'complexity_data', v, 1);
    if iscell(cx) && numel(cx) >= 3
        vals.df_v = nanmed(cx{1}); vals.ri_v = nanmed(cx{2}); vals.oi_v = nanmed(cx{3});
        if isfinite(vals.pacing_hz) && vals.pacing_hz > 0
            vals.capture_ratio = vals.df_v / vals.pacing_hz;
        end
    end

    % ---- alternans ----------------------------------------------------------
    AP  = slot(em, 'alternans_data', v, 1);
    CA  = slot(em, 'alternans_data', ca, 1);
    sub = slot(em, 'alternans_data', ca, 5);          % combined Vm+Ca substrate
    if ~isstruct(sub), sub = slot(em, 'alternans_data', v, 5); end
    if isstruct(AP) && ~isstruct(CA) && isstruct(slot(em, 'alternans_data', v, 1))
        % single-camera file: only the AP struct exists
    end
    if isstruct(AP)
        vals.apd80_alt_ratio_v = nanmed(fld(AP, 'APD80_ratio_map'));
        vals.apd50_alt_ratio_v = nanmed(fld(AP, 'APD50_ratio_map'));
        vals.amp_alt_ratio_v   = nanmed(fld(AP, 'amp_ratio_map'));
        vals.dvdt_alt_ratio_v  = nanmed(fld(AP, 'dvdt_ratio_map'));
        vals.sai_v             = nanmed(fld(AP, 'spectral_map'));
        if isfield(AP, 'num_beats'), vals.n_beats = double(AP.num_beats); end
        r = fld(AP, 'APD80_ratio_map'); pv = fld(AP, 'pval_map');
        if ~isempty(r) && ~isempty(pv) && isequal(size(r), size(pv))
            ok = isfinite(r);
            if any(ok(:)), vals.alt_sig_frac_v = nnz(ok & r > 0.10 & pv < 0.05) / nnz(ok); end
        end
    end
    if isstruct(CA)
        vals.ca_amp_alt_ratio = nanmed(fld(CA, 'amp_ratio_map'));
        vals.ca_d50_alt_ratio = nanmed(fld(CA, 'D50_ratio_map'));
        vals.ca_d80_alt_ratio = nanmed(fld(CA, 'D80_ratio_map'));
        vals.ca_release_alt   = nanmed(fld(CA, 'release_alt_map'));
        vals.ca_load_alt      = nanmed(fld(CA, 'load_alt_map'));
        vals.ca_sai           = nanmed(fld(CA, 'spectral_map'));
        vals.ca_alt_phase     = nanmed(fld(CA, 'phase_map'));
        r = fld(CA, 'amp_ratio_map'); pv = fld(CA, 'pval_map');
        if ~isempty(r) && ~isempty(pv) && isequal(size(r), size(pv))
            ok = isfinite(r);
            if any(ok(:)), vals.ca_alt_sig_frac = nnz(ok & r > 0.10 & pv < 0.05) / nnz(ok); end
        end
        if isnan(vals.n_beats) && isfield(CA, 'num_beats'), vals.n_beats = double(CA.num_beats); end
    end
    % ---- beat-windowed metrics: per-pixel median over beats, then over pixels ----
    if isstruct(AP)
        lev = fld(AP, 'APD_levels');
        vals.apd50_bw_v = nanmed(beat_median(fld(AP, 'APD_beat'), lev, 50));
        vals.apd80_bw_v = nanmed(beat_median(fld(AP, 'APD_beat'), lev, 80));
    end
    if isstruct(CA)
        lev = fld(CA, 'DecayLevels');
        ttp = beat_stack(fld(CA, 'TTP_beat'), [], []);
        d50 = beat_stack(fld(CA, 'CaTD_beat'), lev, 50);
        d80 = beat_stack(fld(CA, 'CaTD_beat'), lev, 80);
        vals.ca_ttp_bw  = nanmed(med3(ttp));
        vals.ca_d50_bw  = nanmed(med3(d50));
        vals.ca_d80_bw  = nanmed(med3(d80));
        if ~isempty(ttp) && isequal(size(ttp), size(d50)), vals.catd50_bw = nanmed(med3(ttp + d50)); end
        if ~isempty(ttp) && isequal(size(ttp), size(d80)), vals.catd80_bw = nanmed(med3(ttp + d80)); end
        vals.ca_diast_bw = nanmed(med3(beat_stack(fld(CA, 'diastolic_beat'), [], [])));
    end

    % ---- ensemble-window validity flags ------------------------------------
    % Same rule the app guard applies (check_analysis_window: wrapped window OR
    % incomplete relaxation), plus outcome checks the guard cannot make because
    % it never sees the extracted values.
    ok = isfinite(vals.catd_ca) && vals.catd_ca >= 15 && ...
         (~isfinite(vals.vc_delay) || vals.vc_delay >= 0) && ...
         (~isfinite(vals.catd_valid_frac_ca) || vals.catd_valid_frac_ca >= 0.5);
    if ok && isfinite(vals.catd80_bw) && vals.catd80_bw > 0
        ok = abs(vals.catd_ca - vals.catd80_bw) <= 0.5 * vals.catd80_bw;
    end
    ok = ok && window_guard_ok(d, ca) && ...
         ~(isfinite(vals.ca_polarity_conf) && vals.ca_polarity_conf < 0);
    if isfinite(vals.catd_ca), vals.ca_ensemble_valid = double(ok); end
    ok = isfinite(vals.apd_v) && vals.apd_v >= 10 && ...
         (~isfinite(vals.apd_valid_frac_v) || vals.apd_valid_frac_v >= 0.5);
    if ok && isfinite(vals.apd80_bw_v) && vals.apd80_bw_v > 0
        ok = abs(vals.apd_v - vals.apd80_bw_v) <= 0.5 * vals.apd80_bw_v;
    end
    ok = ok && window_guard_ok(d, v) && ...
         ~(isfinite(vals.v_polarity_conf) && vals.v_polarity_conf < 0);
    if isfinite(vals.apd_v), vals.v_ensemble_valid = double(ok); end

    if isstruct(sub)
        vals.apd_alt_ratio_v = nanmed(fld(sub, 'alt_ratio_map'));
        vals.alt_phase_v     = nanmed(fld(sub, 'phase_map'));
        vals.apd_gradient_v  = nanmed(fld(sub, 'gradient_mag'));
        vals.restitution_v   = nanmed(fld(sub, 'restitution_map'));
        vals.concordance_v   = scalar(fld(sub, 'concordance_ratio'));
        vals.ca_ap_coupling  = nanmed(fld(sub, 'coupling_map'));
        vals.ca_ap_inphase   = scalar(fld(sub, 'coupling_fraction'));
        vals.risk_map        = nanmed(fld(sub, 'risk_map'));
        vals.risk_global     = scalar(fld(sub, 'risk_global'));
    end
end


% -------------------------------------------------------------------------
function S = beat_stack(C, levels, level)
% R x C x nBeats stack of per-beat maps from a {nBeats x nLevels} (or {nBeats x 1}) cell.
    S = [];
    if ~iscell(C) || isempty(C), return; end
    if isempty(level)
        col = 1;
    else
        col = find(abs(double(levels(:)') - level) < 1e-6, 1);
        if isempty(col) || col > size(C, 2), return; end
    end
    maps = C(:, col);
    keep = cellfun(@(x) isnumeric(x) && ismatrix(x) && ~isempty(x), maps);
    if ~any(keep), return; end
    S = cat(3, maps{keep});
end

function m = med3(S)
    if isempty(S), m = []; else, m = median(double(S), 3, 'omitnan'); end
end

function m = beat_median(C, levels, level)
    m = med3(beat_stack(C, levels, level));
end

function ok = window_guard_ok(d, cam)
% The app-side guard, applied to the stored ensemble beat.  Missing or
% unreadable data leaves the verdict untouched (true), exactly as in the app.
    ok = true;
    f_avg = sprintf('CAM%d_average', cam);
    if ~isfield(d, f_avg) || isempty(d.(f_avg)) || exist('check_analysis_window', 'file') ~= 2
        return;
    end
    try
        mk = [];
        f_snr = sprintf('CAM%d_SNR', cam);
        if isfield(d, f_snr) && ~isempty(d.(f_snr))
            mk = create_snr_mask(d.(f_snr), [], true);
        end
        ok = check_analysis_window(d.(f_avg), mk, cam, 'ensemble metrics');
    catch
        ok = true;
    end
end

function m = level_map(d, em, field, cam, rep_val, label)
% Masked duration map at an extra level: from the stored field when present,
% otherwise computed on the fly with the same engine (older metrics files).
    m = slot(em, field, cam, 1);
    if ~isempty(m), return; end
    f_avg = sprintf('CAM%d_average', cam);  f_snr = sprintf('CAM%d_SNR', cam);
    if ~isfield(d, f_avg) || isempty(d.(f_avg)) || ~isfield(d, 'acqFreq'), return; end
    try
        P = struct('acqFreq', double(d.acqFreq), 'combo_masks', {cell(4,1)}, 'Log', @(s) []);
        [r, c] = size(d.(f_avg)(:,:,1));  P.row = round(r/2);  P.col = round(c/2);
        if isfield(d, f_snr) && ~isempty(d.(f_snr))
            mk = double(create_snr_mask(d.(f_snr), [], true)); mk(~mk) = NaN;
            P.combo_masks{cam,1} = mk;
        end
        a = cadence_extract_apd(d, P, cam, rep_val, label);
        m = a{1,1};
    catch
        m = [];
    end
end

function x = slot(em, field, cam, k)
    x = [];
    if ~isfield(em, field), return; end
    c = em.(field);
    if ~iscell(c) || size(c,1) < cam || size(c,2) < k, return; end
    x = c{cam, k};
end

function x = fld(s, f)
    if isstruct(s) && isfield(s, f), x = s.(f); else, x = []; end
end

function m = nanmed(x)
    m = NaN;
    if isempty(x) || ~isnumeric(x) && ~islogical(x), return; end
    x = double(x(:));
    x = x(isfinite(x));
    if ~isempty(x), m = median(x); end
end

function s = scalar(x)
    s = NaN;
    if isnumeric(x) && isscalar(x), s = double(x); end
end

function m = interior_only(m, margin)
% Same gate as Signal Analysis get_cv_curated_data / CV interior_valid_mask.
    if isempty(m) || margin <= 0, return; end
    try
        tissue   = imfill(isfinite(m), 'holes');
        interior = imerode(tissue, true(2*margin + 1));
        m(~interior) = NaN;
    catch
    end
end
