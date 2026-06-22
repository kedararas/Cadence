function substrate = assess_arrhythmia_substrate(ap_result, varargin)
%ASSESS_ARRHYTHMIA_SUBSTRATE  Compute arrhythmia susceptibility metrics from
%   optical-mapping alternans analysis results.
%
%   substrate = assess_arrhythmia_substrate(ap_result)
%   substrate = assess_arrhythmia_substrate(ap_result, Name, Value, ...)
%
%   Inputs
%     ap_result   - struct from analyzeAPAlternans(...) in 3-D mode
%                   (must contain APD_beat, act_beat, APD_levels, dt, num_beats)
%
%   Optional Name-Value Pairs
%     'CaResult'    struct from analyzeCaTransientAlternans in 3-D mode.
%                   Required for Ca-AP phase relationship map.
%     'CaData'      [rows x cols x frames] raw Ca fluorescence array.
%                   Required for diastolic Ca elevation map.
%     'BeatFrames'  [num_beats x 2] [start_frame, end_frame] per beat.
%                   Required for restitution slope and diastolic Ca maps.
%     'APD_level'   Which APD level to use as primary metric (default: 80).
%     'MinAlt'      Minimum APD alternans fraction to flag a pixel as
%                   significant (default: 0.05 = 5% of mean APD).
%
%   Output struct fields
%
%   Phase / concordance
%     .phase_map          [R x C] AP alternans phase: +1 odd>even, -1 odd<even, 0 invalid
%     .phase_angle_map    [R x C] continuous AP phase angle (radians, −π to +π);
%                                 angle of FFT coeff at 0.5 cyc/beat — reveals
%                                 gradual phase gradients and travelling alternans waves
%     .is_discordant      logical scalar: true when >20% of pixels are in minority phase
%     .concordance_ratio  fraction of valid pixels sharing the majority phase [0, 1]
%     .nodal_lines        [R x C] logical: pixels at phase boundaries (highest risk sites)
%
%   Alternans ratio
%     .alt_ratio_map      [R x C] |APD_alt| / mean_APD per pixel (dimensionless)
%                                 normalised alternans magnitude independent of local APD
%
%   Spatial dispersion of repolarization
%     .apd_mean_map       [R x C] mean APD across all beats (ms)
%     .apd_disp_global    scalar: SD of APD map across valid pixels (ms)
%     .gradient_mag       [R x C] spatial gradient magnitude of APD map (ms/pixel)
%     .gradient_dir       [R x C] gradient direction (degrees)
%     .max_gradient       scalar: peak gradient magnitude
%
%   Ca-AP coupling (requires CaResult)
%     .coupling_map       [R x C] Pearson r: +1 = in-phase (Ca drives AP),
%                                            -1 = out-of-phase (V drives Ca)
%     .inphase_mask       [R x C] logical: true where Ca and AP alternate in phase
%     .coupling_fraction  fraction of valid pixels that are in-phase
%     .ca_phase_map       [R x C] Ca alternans phase: +1/-1/0 (binary)
%     .ca_phase_angle_map [R x C] Ca continuous phase angle (radians, −π to +π)
%     .ca_alt_ratio_map   [R x C] |Ca_alt| / mean_Ca per pixel
%
%   Diastolic Ca elevation (requires CaData + BeatFrames)
%     .diast_ca_map       [R x C] mean normalised diastolic Ca level
%     .diast_ca_alt_map   [R x C] diastolic Ca alternans (odd - even diastole)
%
%   APD restitution slope (requires BeatFrames + >= 3 beats)
%     .restitution_map    [R x C] local APD restitution slope per pixel
%     .slope_gt1_mask     [R x C] logical: pixels with slope > 1 (dynamically unstable)
%     .mean_slope         scalar: mean slope across valid pixels
%
%   Composite arrhythmia risk
%     .risk_map           [R x C] composite risk score [0, 1] per pixel
%     .risk_global        scalar: mean risk score across significant pixels
%     .risk_components    struct with individual normalised component maps

    %% ---- Parse inputs ---------------------------------------------------
    p = inputParser;
    addRequired(p,  'ap_result',  @isstruct);
    addParameter(p, 'CaResult',   [],   @(x) isempty(x) || isstruct(x));
    addParameter(p, 'CaData',     [],   @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'BeatFrames', [],   @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'APD_level',  'auto', @(x) isnumeric(x) || ischar(x) || isstring(x));
    addParameter(p, 'MinAlt',     0.15, @isnumeric);
    parse(p, ap_result, varargin{:});
    opt = p.Results;

    ca_result   = opt.CaResult;
    ca_data     = opt.CaData;
    beat_frames = opt.BeatFrames;
    apd_level   = opt.APD_level;
    min_alt     = opt.MinAlt;

    % Validate ap_result is from 3-D mode
    if ~isfield(ap_result, 'APD_beat') || ~strcmp(ap_result.mode, '3D')
        error('assess_arrhythmia_substrate:invalid', ...
            'ap_result must be the output of analyzeAPAlternans in 3-D mode.');
    end

    % Resolve 'auto' APD level (default).  Pick the deepest repolarisation
    % level whose spatial coverage is still close to the best-covered level.
    % At fast pacing the AP may not repolarise to 80% within the beat window,
    % leaving APD80 mostly NaN; auto then falls back to APD50/APD30 so the
    % substrate maps are not dominated by holes.  Pass a numeric APD_level to
    % override.
    if ischar(apd_level) || isstring(apd_level)
        apd_level = select_apd_level(ap_result);
    end

    % Get dimensions from first APD map
    apd80_col = find(ap_result.APD_levels == apd_level, 1);
    if isempty(apd80_col)
        warning('assess_arrhythmia_substrate:noLevel', ...
            'APD_level=%d not found in ap_result.APD_levels. Using last level.', apd_level);
        apd80_col = size(ap_result.APD_beat, 2);
        apd_level = ap_result.APD_levels(apd80_col);
    end

    ref_map = ap_result.APD_beat{1, apd80_col};
    [R, C]  = size(ref_map);
    nBeats  = ap_result.num_beats;
    dt      = ap_result.dt;   % ms per frame

    % Stack per-beat APD maps: R x C x nBeats
    apd_stack    = cat(3, ap_result.APD_beat{:, apd80_col});  % ms
    apd_mean_map = mean(apd_stack, 3, 'omitnan');
    apd_alt      = ap_result.(sprintf('APD%d_alt', apd_level));

    % Valid pixel mask: non-NaN alternans and non-NaN mean APD
    valid = ~isnan(apd_alt) & ~isnan(apd_mean_map) & apd_mean_map > 0;

    % Alternans ratio map: |APD_alternans| / mean_APD per pixel.
    % Prefer precomputed map from ap_result to avoid duplicate work.
    apd_ratio_field = sprintf('APD%d_ratio_map', apd_level);
    if isfield(ap_result, apd_ratio_field)
        alt_ratio_map = ap_result.(apd_ratio_field);
    else
        alt_ratio_map = abs(apd_alt) ./ (apd_mean_map + eps);
    end
    alt_ratio_map(~valid) = nan;


    % =====================================================================
    %  1. PHASE MAP & CONCORDANCE / DISCORDANT ALTERNANS
    % =====================================================================
    % Binary phase map  (+1 / -1 / 0)
    %   +1 : odd beats have longer APD than even beats (positive phase)
    %   -1 : odd beats have shorter APD (negative / inverted phase)
    %    0 : invalid pixel
    % Use precomputed map from ap_result if available (avoids recomputing
    % the FFT), otherwise derive from the alternans sign.
    %
    % Continuous phase angle map  (radians, −π to +π)
    %   Computed as angle of the FFT coefficient at 0.5 cycles/beat.
    %   Gradual spatial gradients = travelling alternans wave.
    %   Pixels separated by ~π radians are in anti-phase (discordant).
    if isfield(ap_result, 'phase_map')
        phase_map = ap_result.phase_map;
    else
        phase_map          = zeros(R, C);
        phase_map(valid)   = sign(apd_alt(valid));
    end

    if isfield(ap_result, 'phase_angle_map')
        phase_angle_map = ap_result.phase_angle_map;
    else
        % Compute from per-beat APD stack
        phase_angle_map = nan(R, C);
        apd80_col_loc   = apd80_col;
        if nBeats >= 4
            all_apd = cat(3, ap_result.APD_beat{:, apd80_col_loc});
            all_apd(isnan(all_apd)) = 0;
            all_apd = all_apd - mean(all_apd, 3);
            Y_ap    = fft(all_apd, [], 3);
            alt_bin = floor(nBeats/2) + 1;
            if alt_bin <= size(Y_ap, 3)
                phase_angle_map = angle(Y_ap(:,:, alt_bin));
            end
        end
    end

    % Fraction significant (|alt| > min_alt * mean_APD)
    sig_mask = valid & abs(apd_alt) > min_alt * apd_mean_map;

    n_sig = sum(sig_mask(:));
    if n_sig > 0
        n_pos = sum(phase_map(sig_mask) > 0);
        n_neg = sum(phase_map(sig_mask) < 0);
        majority_n       = max(n_pos, n_neg);
        concordance_ratio = majority_n / n_sig;
        is_discordant    = concordance_ratio < 0.80;  % >20% in minority phase
    else
        concordance_ratio = 1;
        is_discordant    = false;
    end

    % Nodal lines: pixels bordering both positive and negative phase regions
    se         = strel('disk', 1);
    pos_region = phase_map > 0;
    neg_region = phase_map < 0;
    nodal_lines = imdilate(pos_region, se) & imdilate(neg_region, se);


    % =====================================================================
    %  2. SPATIAL DISPERSION OF REPOLARISATION
    % =====================================================================
    % Large SD of APD map = large spatial heterogeneity = reentry substrate.
    % High gradient magnitude = steep repolarisation border = block site.

    apd_disp_global = std(apd_mean_map(valid), 'omitnan');

    % Gradient of mean APD map
    apd_clean = apd_mean_map;
    apd_clean(~valid) = 0;                     % zero invalid for gradient calc
    [gradient_mag, gradient_dir] = imgradient(apd_clean, 'sobel');
    gradient_mag(~valid) = nan;                % restore invalid mask
    gradient_dir(~valid) = nan;
    max_gradient = max(gradient_mag(:), [], 'omitnan');


    % =====================================================================
    %  3. Ca-AP COUPLING MAP  (requires CaResult)
    % =====================================================================
    % Pearson correlation between the per-beat Ca amplitude series and the
    % per-beat APD series at each pixel.
    %
    %   r > 0  →  in-phase: large Ca transient → long APD
    %             Ca drives AP via NCX / ICaL inactivation
    %   r < 0  →  out-of-phase: large Ca → short APD
    %             Voltage drives Ca via ICaL gating
    %
    % Out-of-phase regions are less mechanistically driven by Ca and may
    % represent a more complex arrhythmia substrate.

    coupling_map      = nan(R, C);
    inphase_mask      = false(R, C);
    coupling_fraction = nan;
    ca_phase_map       = nan(R, C);     % binary  +1/-1/0
    ca_phase_angle_map = nan(R, C);     % continuous  [-pi, pi]
    ca_alt_ratio_map   = nan(R, C);     % |Ca_alt| / mean_Ca per pixel

    if ~isempty(ca_result) && isfield(ca_result, 'amp_beat') && ...
            ca_result.num_beats == nBeats

        ca_series   = cat(3, ca_result.amp_beat{:});     % R x C x nBeats
        odd_idx     = 1:2:nBeats;
        even_idx    = 2:2:nBeats;
        n_pairs_ca  = min(numel(odd_idx), numel(even_idx));
        ap_series  = apd_stack;                           % R x C x nBeats

        % Per-pixel zero-mean normalisation
        ca_dm = ca_series - mean(ca_series, 3, 'omitnan');
        ap_dm = ap_series - mean(ap_series, 3, 'omitnan');

        ca_std = std(ca_series, 0, 3, 'omitnan');
        ap_std = std(ap_series, 0, 3, 'omitnan');

        coupling_map = mean(ca_dm .* ap_dm, 3, 'omitnan') ./ ...
                       (ca_std .* ap_std + eps);

        % Invalidate pixels with negligible signal variance
        no_signal = ca_std < eps('single') | ap_std < eps('single') | ~valid;
        coupling_map(no_signal) = nan;

        inphase_mask = coupling_map > 0 & ~no_signal;

        n_valid_coup = sum(~isnan(coupling_map(:)));
        if n_valid_coup > 0
            coupling_fraction = sum(inphase_mask(:)) / n_valid_coup;
        end

        % Ca phase maps — wire in from ca_result if precomputed
        if isfield(ca_result, 'phase_map')
            ca_phase_map = ca_result.phase_map;
        else
            ca_phase_map = sign(mean(ca_series, 3, 'omitnan') - ...
                               mean(ca_series(:,:,2:2:end), 3, 'omitnan'));
        end
        if isfield(ca_result, 'phase_angle_map')
            ca_phase_angle_map = ca_result.phase_angle_map;
        end

        % Ca alternans ratio map
        if isfield(ca_result, 'amp_ratio_map')
            ca_alt_ratio_map = ca_result.amp_ratio_map;
        else
            ca_amp_alt  = mean(cat(3, ca_result.amp_beat{odd_idx(1:n_pairs_ca)}) - ...
                               cat(3, ca_result.amp_beat{even_idx(1:n_pairs_ca)}), ...
                               3, 'omitnan');
            ca_mean_map = mean(cat(3, ca_result.amp_beat{:}), 3, 'omitnan');
            ca_alt_ratio_map = abs(ca_amp_alt) ./ (ca_mean_map + eps);
        end

        fprintf('Ca-AP coupling: %.1f%% of pixels in-phase ', coupling_fraction*100);
        if coupling_fraction > 0.6
            fprintf('(Ca-driven alternans likely)\n');
        elseif coupling_fraction < 0.4
            fprintf('(Voltage-driven alternans likely)\n');
        else
            fprintf('(Mixed / transitional)\n');
        end
    else
        if isempty(ca_result)
            fprintf('Ca-AP coupling: skipped (no CaResult supplied)\n');
        else
            fprintf('Ca-AP coupling: skipped (beat count mismatch: AP=%d, Ca=%d)\n', ...
                nBeats, ca_result.num_beats);
        end
    end


    % =====================================================================
    %  4. DIASTOLIC Ca ELEVATION MAP  (requires CaData + BeatFrames)
    % =====================================================================
    % Elevated diastolic Ca between beats → delayed afterdepolarisations
    % (DADs) → triggered activity independent of alternans mechanism.
    % Diastolic Ca alternans → beat-to-beat variation in SR loading → SR
    % instability → amplitude alternans origin.

    diast_ca_map     = nan(R, C);
    diast_ca_alt_map = nan(R, C);

    if ~isempty(ca_data) && ~isempty(beat_frames)
        nBF      = size(beat_frames, 1);
        ca_d     = double(ca_data);

        % Per-pixel [0,1] normalise once so diastolic levels are comparable
        ca_lo  = min(ca_d, [], 3);
        ca_hi  = max(ca_d, [], 3);
        ca_rng = max(ca_hi - ca_lo, eps('single'));
        ca_norm = (ca_d - ca_lo) ./ ca_rng;   % R x C x T

        % Diastolic window: last 10% of each inter-beat interval
        diast_ca = nan(R, C, nBF - 1);
        for j = 1:nBF - 1
            diast_end   = beat_frames(j+1, 1) - 1;
            beat_dur    = beat_frames(j, 2) - beat_frames(j, 1) + 1;
            diast_win   = max(2, round(0.10 * beat_dur));
            diast_start = max(beat_frames(j, 2) + 1, diast_end - diast_win + 1);
            if diast_start <= diast_end && diast_end <= size(ca_norm, 3)
                diast_ca(:,:,j) = mean(ca_norm(:,:, diast_start:diast_end), 3);
            end
        end

        diast_ca_map = mean(diast_ca, 3, 'omitnan');

        % Diastolic Ca alternans: odd vs even inter-beat intervals
        odd_d  = diast_ca(:,:, 1:2:end);
        even_d = diast_ca(:,:, 2:2:end);
        n_dp   = min(size(odd_d,3), size(even_d,3));
        if n_dp >= 1
            diast_ca_alt_map = mean(odd_d(:,:,1:n_dp) - even_d(:,:,1:n_dp), ...
                                    3, 'omitnan');
        end

        fprintf('Diastolic Ca: mean=%.3f  range=[%.3f, %.3f]\n', ...
            mean(diast_ca_map(:),'omitnan'), ...
            min(diast_ca_map(:),[],'omitnan'), ...
            max(diast_ca_map(:),[],'omitnan'));
    else
        fprintf('Diastolic Ca: skipped (no CaData or BeatFrames supplied)\n');
    end


    % =====================================================================
    %  5. APD RESTITUTION SLOPE MAP  (requires BeatFrames + >= 3 beats)
    % =====================================================================
    % Local restitution slope: dAPD/dDI estimated from beat-to-beat data.
    %   Slope > 1  →  dynamically unstable → alternans / VF susceptible
    %   Slope < 1  →  stable
    %
    % With alternans present, odd/even beats naturally sample two different
    % (APD, DI) pairs per pixel, giving at least 2 points for slope estimation.
    % A full restitution curve requires an S1-S2 protocol; this provides
    % a local slope estimate at the operating diastolic interval.

    restitution_map = nan(R, C);
    slope_gt1_mask  = false(R, C);
    mean_slope      = nan;

    if ~isempty(beat_frames) && nBeats >= 3
        % Global activation frame for each beat (per pixel)
        % act_beat{j} = activation time in ms from window start
        act_stack = nan(R, C, nBeats);
        for j = 1:nBeats
            act_ms = ap_result.act_beat{j};          % R x C, ms from window start
            act_stack(:,:,j) = beat_frames(j,1) + act_ms/dt - 1;  % global frame
        end

        % Global repolarisation frame = activation frame + APD/dt
        rep_stack = act_stack + apd_stack / dt;      % R x C x nBeats (global frames)

        % DI(j) = act_global(j+1) - rep_global(j)   [frames → ms]
        DI_stack  = (act_stack(:,:,2:end) - rep_stack(:,:,1:end-1)) * dt;
        DI_stack  = max(DI_stack, 0);                % DI cannot be negative
        APD_next  = apd_stack(:,:,2:end);            % APD of beat j+1

        % Vectorised OLS slope: dAPD/dDI per pixel
        % slope = cov(DI, APD) / var(DI)
        DI_mean  = mean(DI_stack,  3, 'omitnan');
        APD_mean = mean(APD_next,  3, 'omitnan');
        DI_dm    = DI_stack  - DI_mean;
        APD_dm   = APD_next  - APD_mean;
        cov_da   = mean(DI_dm .* APD_dm, 3, 'omitnan');
        var_di   = mean(DI_dm .^ 2,      3, 'omitnan');

        restitution_map = cov_da ./ (var_di + eps);
        % Invalidate pixels with insufficient DI variance (constant pacing)
        restitution_map(var_di < 1e-6 | ~valid) = nan;

        slope_gt1_mask = restitution_map > 1 & ~isnan(restitution_map);
        mean_slope     = mean(restitution_map(:), 'omitnan');

        pct_unstable = 100 * sum(slope_gt1_mask(:)) / max(sum(~isnan(restitution_map(:))),1);
        fprintf('Restitution slope: mean=%.3f  %.1f%% of pixels > 1 (unstable)\n', ...
            mean_slope, pct_unstable);
    else
        if isempty(beat_frames)
            fprintf('Restitution slope: skipped (no BeatFrames supplied)\n');
        else
            fprintf('Restitution slope: skipped (need >= 3 beats, have %d)\n', nBeats);
        end
    end


    % =====================================================================
    %  6. COMPOSITE ARRHYTHMIA RISK MAP
    % =====================================================================
    % Each component is normalised to [0, 1] and weighted:
    %
    %   Component             Weight  Rationale
    %   ─────────────────     ──────  ─────────────────────────────────────
    %   Alternans fraction    1.0     Primary substrate for dispersion
    %   APD gradient          1.0     Steepness of repolarisation border
    %   Discordant / nodal    1.0     Directly predicts block site
    %   Restitution slope>1   1.0     Dynamic instability predictor
    %   Ca-AP out-of-phase    0.5     Mechanistic amplifier of risk
    %   High diastolic Ca     0.5     Independent triggered-activity risk
    %
    %   Maximum possible score: 5.0  (or 4.0 without Ca data)
    %   risk_map normalised to [0, 1] by dividing by maximum.

    max_score = 4.0;   % without Ca data
    risk_map  = zeros(R, C);

    % Component 1: alternans ratio  (saturates at 20%)
    % Use the already-computed alt_ratio_map (avoids recomputing the division)
    c1 = min(alt_ratio_map / 0.20, 1);
    c1(isnan(c1) | ~valid) = 0;

    % Component 2: APD gradient (normalised to 95th percentile)
    grad_p95  = prctile(gradient_mag(~isnan(gradient_mag)), 95);
    c2        = min(gradient_mag / (grad_p95 + eps), 1);
    c2(isnan(c2)) = 0;

    % Component 3: nodal lines (binary — highest risk pixels)
    c3 = double(nodal_lines);

    % Component 4: restitution slope > 1 (binary)
    c4 = double(slope_gt1_mask);

    risk_map = c1 + c2 + c3 + c4;

    % Component 5: Ca-AP out-of-phase (adds 0.5 weight)
    c5 = zeros(R, C);
    if ~all(isnan(coupling_map(:)))
        c5 = max(-coupling_map, 0);   % 0 where in-phase, up to 1 where fully out-of-phase
        c5(isnan(coupling_map)) = 0;
        risk_map = risk_map + 0.5 * c5;
        max_score = max_score + 0.5;
    end

    % Component 6: elevated diastolic Ca (adds 0.5 weight)
    c6 = zeros(R, C);
    if ~all(isnan(diast_ca_map(:)))
        diast_p95 = prctile(diast_ca_map(~isnan(diast_ca_map)), 95);
        c6 = min(diast_ca_map / (diast_p95 + eps), 1);
        c6(isnan(diast_ca_map)) = 0;
        risk_map = risk_map + 0.5 * c6;
        max_score = max_score + 0.5;
    end

    % Normalise to [0, 1]
    risk_map(~valid) = nan;
    risk_map(valid)  = risk_map(valid) / max_score;
    risk_global      = mean(risk_map(sig_mask), 'omitnan');


    % =====================================================================
    %  Assemble output
    % =====================================================================
    % ── Phase / concordance ──────────────────────────────────────────────
    substrate.phase_map          = phase_map;         % AP: binary +1/-1/0
    substrate.phase_angle_map    = phase_angle_map;   % AP: continuous [-pi,pi]
    substrate.is_discordant      = is_discordant;
    substrate.concordance_ratio  = concordance_ratio;
    substrate.nodal_lines        = nodal_lines;

    % ── Alternans ratio map ───────────────────────────────────────────────
    substrate.alt_ratio_map      = alt_ratio_map;     % |APD_alt| / mean_APD

    % ── Spatial dispersion ───────────────────────────────────────────────
    substrate.apd_mean_map       = apd_mean_map;
    substrate.apd_disp_global    = apd_disp_global;
    substrate.gradient_mag       = gradient_mag;
    substrate.gradient_dir       = gradient_dir;
    substrate.max_gradient       = max_gradient;

    % ── Ca-AP coupling ───────────────────────────────────────────────────
    substrate.coupling_map       = coupling_map;
    substrate.inphase_mask       = inphase_mask;
    substrate.coupling_fraction  = coupling_fraction;
    substrate.ca_phase_map       = ca_phase_map;       % Ca: binary +1/-1/0
    substrate.ca_phase_angle_map = ca_phase_angle_map; % Ca: continuous [-pi,pi]
    substrate.ca_alt_ratio_map   = ca_alt_ratio_map;   % |Ca_alt| / mean_Ca

    % ── Diastolic Ca ─────────────────────────────────────────────────────
    substrate.diast_ca_map       = diast_ca_map;
    substrate.diast_ca_alt_map   = diast_ca_alt_map;

    % ── Restitution slope ────────────────────────────────────────────────
    substrate.restitution_map    = restitution_map;
    substrate.slope_gt1_mask     = slope_gt1_mask;
    substrate.mean_slope         = mean_slope;

    % ── Composite risk ───────────────────────────────────────────────────
    substrate.risk_map           = risk_map;
    substrate.risk_global        = risk_global;
    substrate.risk_components    = struct( ...
        'alternans_fraction', c1, ...
        'apd_gradient',       c2, ...
        'nodal_proximity',    c3, ...
        'restitution_slope',  c4, ...
        'ca_outofphase',      c5, ...
        'diastolic_ca',       c6);

    % ── Print summary ─────────────────────────────────────────────────────
    fprintf('\n======== Arrhythmia Substrate Assessment ========\n');
    fprintf('Significant pixels    : %d / %d (%.1f%%)\n', ...
        n_sig, sum(valid(:)), 100*n_sig/max(sum(valid(:)),1));
    if is_discordant
        fprintf('Concordance           : DISCORDANT (%.1f%% majority phase) *** HIGH RISK\n', ...
            concordance_ratio*100);
    else
        fprintf('Concordance           : concordant (%.1f%% majority phase)\n', ...
            concordance_ratio*100);
    end
    fprintf('APD spatial dispersion: %.2f ms\n', apd_disp_global);
    fprintf('Max APD gradient      : %.2f ms/pixel\n', max_gradient);
    if ~isnan(mean_slope)
        fprintf('Mean restitution slope: %.3f  (%s)\n', mean_slope, ...
            ternary(mean_slope > 1, 'UNSTABLE > 1', 'stable < 1'));
    end
    fprintf('Global risk score     : %.3f / 1.0\n', risk_global);
    fprintf('=================================================\n\n');
end


% =========================================================================
%  Local helpers
% =========================================================================

function out = ternary(cond, a, b)
%TERNARY  Inline conditional: return a if cond is true, else b.
    if cond, out = a; else, out = b; end
end


function lvl = select_apd_level(ap_result)
%SELECT_APD_LEVEL  Deepest APD level with adequate spatial coverage.
%   Coverage = fraction of tissue pixels with a finite APD-alternans value at
%   that level.  Returns the deepest (highest %) level whose coverage is at
%   least 75% of the best-covered level's, so normal recordings still resolve
%   to APD80 while fast-rate recordings (APD80 mostly NaN) fall back to a
%   level that is actually reachable.
    levels = ap_result.APD_levels(:)';
    if isfield(ap_result, 'amp_alt_map')
        tissue_ref = isfinite(ap_result.amp_alt_map);   % full tissue extent
    else
        tissue_ref = false(size(ap_result.APD_beat{1,1}));
        for L = levels
            tissue_ref = tissue_ref | isfinite(ap_result.(sprintf('APD%d_alt', L)));
        end
    end
    nTissue = max(nnz(tissue_ref), 1);
    cov = zeros(size(levels));
    for k = 1:numel(levels)
        am     = ap_result.(sprintf('APD%d_alt', levels(k)));
        cov(k) = nnz(isfinite(am) & tissue_ref) / nTissue;
    end
    cand = levels(cov >= 0.75 * max(cov));   % adequately-covered levels
    lvl  = max(cand);                        % deepest of those
    fprintf(['assess_arrhythmia_substrate: APD_level=auto -> APD%d ' ...
             '(coverage %.0f%%; deepest APD%d only %.0f%%)\n'], ...
             lvl, 100*cov(levels==lvl), max(levels), 100*cov(levels==max(levels)));
end
