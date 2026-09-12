function data = cadence_extract_apd(d, P, camera, rep_val, label)
%CADENCE_EXTRACT_APD  Per-pixel duration / repolarization at one level (port of the app's extract_apd).
%
%   data = cadence_extract_apd(d, P, camera, rep_val, label)
%
%   d        cmos_all_data (needs CAM<camera>_average and CAM<camera>_image)
%   P        struct with fields acqFreq, row, col, combo_masks (cell per camera), Log
%   camera   camera index
%   rep_val  fraction of amplitude REMAINING at the level (APD80 -> 0.2)
%   label    text for log lines
%
%   Returns the app's 9-slot cell: {1} masked duration map (ms), {2} background
%   image, {3} pixel trace, {4} fitted trace, {5} camera label, {6} unmasked
%   duration, {7} masked repolarization time (ms), {8} pixel repolarization,
%   {9} unmasked repolarization.  Ensemble-average path only.
% Port of Cadence_Feature_Extraction/extract_apd (ensemble path).
    row = P.row; col = P.col;
    win_data = analysis_window(d, camera);

    [~, ~, T] = size(win_data);
    win_norm  = normalize_data(win_data, [0 100]);
    t_idx     = reshape(1:T, 1, 1, T);

    act_times = compute_lat_50(win_data);     % frames, sub-frame, NaN where invalid
    act3      = reshape(act_times, size(act_times,1), size(act_times,2), 1);

    pre_act                = win_norm;
    pre_act(t_idx >= act3) = NaN;
    base                   = median(pre_act, 3, 'omitnan');  % diastolic level

    n_pre   = sum(~isnan(pre_act), 3);
    min_pre = max(3, round(0.005 * P.acqFreq));
    base(n_pre < min_pre) = NaN;

    post_act               = win_norm;
    post_act(t_idx < act3) = NaN;
    [pk_val, pk]           = max(post_act, [], 3);
    amp                    = pk_val - base;

    % incomplete-recovery check (ensemble window)
    tail_n   = max(3, round(0.1 * T));
    tail_lvl = median(win_norm(:,:,end-tail_n+1:end), 3);
    base_err = (base - tail_lvl) ./ amp;
    base_err(isnan(act_times)) = NaN;
    med_err  = median(abs(base_err(:)), 'omitnan');
    if med_err > 0.05
        P.Log(sprintf(['WARNING CAM%d %s: pre-activation baseline differs from end-window diastole ' ...
            'by %.1f%% of amplitude (median); APD may be biased short.'], camera, label, 100 * med_err));
    end

    thr = base + rep_val .* amp;

    [hit, rep_frame] = max(cumsum((win_norm <= thr) & (t_idx > pk), 3) == 1, [], 3);
    rep              = double(rep_frame);

    [Rn, Cn]   = size(rep_frame);
    [rr, cc]   = ndgrid(1:Rn, 1:Cn);
    prevf      = max(rep_frame - 1, 1);
    v_hi       = win_norm(sub2ind([Rn Cn T], rr, cc, rep_frame));
    v_lo       = win_norm(sub2ind([Rn Cn T], rr, cc, prevf));
    den        = v_lo - v_hi;
    ok         = hit & (rep_frame > 1) & (den > 0);
    frac       = zeros(Rn, Cn);
    frac(ok)   = (v_lo(ok) - thr(ok)) ./ den(ok);
    rep(ok)    = (rep_frame(ok) - 1) + frac(ok);

    rep(~hit | ~(amp > 0) | isnan(base)) = NaN;

    apd         = rep - act_times;
    bad         = apd <= 0 | isnan(apd);
    apd(bad)    = NaN;
    rep(bad)    = NaN;

    oap        = squeeze(win_norm(row, col, :));
    act_loc    = round(act_times(row, col));
    rep_loc    = round(rep(row, col));
    fitted_oap = nan(T, 1);
    if ~isnan(act_loc) && ~isnan(rep_loc)
        fitted_oap(act_loc:rep_loc) = oap(act_loc:rep_loc);
    end

    rep = rep * (1000 / P.acqFreq);
    apd = apd * (1000 / P.acqFreq);

    data = cell(1, 9);
    data{1,1} = apd;  data{1,7} = rep;
    if ~isempty(P.combo_masks{camera,1})
        m = ~isnan(act_times .* P.combo_masks{camera,1});
        data{1,1}(~m) = NaN;
        data{1,7}(~m) = NaN;
    end
    data{1,2} = get_bg_image(d, camera);
    data{1,3} = oap;
    data{1,4} = fitted_oap;
    data{1,5} = num2str(camera);
    data{1,6} = apd;
    data{1,8} = rep(row, col);
    data{1,9} = rep;

    elig = ~isnan(act_times);
    if ~isempty(P.combo_masks{camera,1})
        elig = elig & ~isnan(P.combo_masks{camera,1});
    end
    n_all = nnz(elig);
    n_bad = nnz(elig & isnan(data{1,1}));
    P.Log(sprintf('CAM%d %s: %d of %d activated pixels (%.1f%%) censored (no crossing or duration<=0).', ...
        camera, label, n_bad, n_all, 100 * n_bad / max(n_all, 1)));
end




function win = analysis_window(d, cam)
    f = sprintf('CAM%d_average', cam);
    if ~isfield(d, f) || isempty(d.(f))
        error('cadence_extract_apd:noAverage', 'CAM%d has no ensemble average (CAM%d_average).', cam, cam);
    end
    win = d.(f);
end

function bg = get_bg_image(d, cam)
    bg = d.(sprintf('CAM%d_image', cam));
end
