function [d, status] = cadence_condition_data(d, opts)
%CADENCE_CONDITION_DATA  Headless port of the Signal Conditioning app's EXECUTE.
%
%   [cmos_all_data, status] = cadence_condition_data(cmos_all_data)
%   [cmos_all_data, status] = cadence_condition_data(cmos_all_data, opts)
%
%   Runs, without the UI, the same per-recording stage chain that
%   Cadence_Signal_Conditioning.mlapp runs when EXECUTE SIGNAL CONDITIONING is
%   pressed, in the same order, with the same helpers and the same fixed
%   parameters, and returns the struct the app would save as
%   <name>-conditioned.mat (conditioned stamp and qc log included).
%
%   Stage order (the app's), each on every CAM<n>:
%     1  SNR map            extract_snr_mask            -> CAM<n>_SNR   (always)
%     2  drift correction   remove_Drift                               [Drift]
%     3  SVD denoising      denoise_svd(., SVDRank); a camera the SVD
%                           cannot be applied to gets 5 x 5 binning     [SVD]
%                           -> CAM<n>_svd_rank (rank kept; 0 = binning fallback)
%     4  spatial binning    binning(., BinSize)                        [Binning]
%     5  temporal filter    filter_data(., acqFreq, FilterHz)          [FilterHz]
%     6  motion correction  trackCardiacMotion to the frame before
%                           the first beat                              [Motion]
%     7  normalization      normalize_data                             [Normalize]
%     8  polarity check     check_polarity_paced over an SNR>=3 mask,
%                           falling back to check_signal_inversion when
%                           the paced test is not confident (|conf|<0.2);
%                           inverted cameras are flipped 1 - x           (always)
%     9  ensemble average   ensembleAverageFull(., analog1 or auto-
%                           detected peaks)             -> CAM<n>_average [Ensemble]
%    10  SNR mask           create_snr_mask(CAM<n>_SNR, MaskFloor, true);
%                           pixels outside are set to NaN                (always)
%    11  provenance         .conditioned = true, validate_conditioned,
%                           qc_attach                                    (always)
%
%   opts fields (all optional; defaults are the app's dropdown defaults, the
%   "recommended settings" of its welcome text)
%     Drift        true      D) DRIFT CORRECTION
%     SVD          true      A) SPATIAL DENOISING (SVD, rank SVDRank = 8)
%     SVDRank      8
%     Binning      false     B) SPATIAL BINNING (BinSize 3 | 5 | 7 | 9)
%     BinSize      3
%     FilterHz     50        C) TEMPORAL FILTERING [0, FilterHz]; [] or 0 = NONE
%     LowBandHz    []        HYBRID FILTER (batch only). When set (e.g. 50 with
%                            FilterHz 100), a second ensemble average is made
%                            from the stack low-passed at [0, LowBandHz] and
%                            saved as CAM<n>_average_lowband. It is taken before
%                            the final SNR mask, from the same normalized,
%                            polarity-corrected stack and beat windows as
%                            CAM<n>_average, so it equals the average a
%                            FilterHz = LowBandHz run would produce (up to the
%                            cascade of the two low-passes). cadence_extract_features
%                            then measures durations/activation/CV/alternans in
%                            the low band and rise times/DF/RI/OI in the
%                            FilterHz band. Requires Ensemble and LowBandHz < FilterHz.
%     Motion       false     G) MOTION ARTIFACT CORRECTION
%     Normalize    true      E) DATA NORMALIZATION
%     Ensemble     true      F) ENSEMBLE AVERAGING
%     MaskFloor    2         SNR floor of the final tissue mask (app: fixed 2)
%     VoltageCams  []        camera roles, for the Vm-Ca polarity tie-breaker
%     CalciumCams  []        (batch only; the app has no equivalent).  When a
%                            voltage camera's polarity came from the SHAPE check
%                            (no usable stimulus: arrhythmia or unpaced) and a
%                            calcium camera is at least as decisive, the voltage
%                            orientation is set so that calcium FOLLOWS voltage
%                            by -2 to +min(15 ms, 0.3 cycle): the tissue-mean
%                            traces are cross-correlated, and the orientation
%                            that wins by >= 0.2 in correlation is kept.  An
%                            inverted voltage puts the best lag half a cycle
%                            away (Fig 5: -16 / +24 ms as stored, +5 ms flipped).
%                            The source is then 'vm-ca lag'.
%     PolarityPrior []       session polarity prior (batch only): struct array
%                            with fields cam, inverted, n, agree -- the
%                            STIMULUS-based decisions for that camera in the
%                            other recordings of the same experiment (see
%                            cadence_batch_pipeline).  Polarity is set by the dye
%                            and optics, not the rhythm, so a camera whose own
%                            check fell back to signal SHAPE (no usable stimulus)
%                            takes the prior's orientation instead; the source is
%                            then 'session'.  Applied before the Vm-Ca tie-breaker,
%                            which then only sees cameras without a prior.
%     Name         ''        recording name for QC records and messages
%     Log          @(s) fprintf('%s\n', s)
%
%   Like the app, a stage that throws is logged and the data passes through it
%   UNCHANGED; the error is recorded in status.errors so a batch caller can
%   refuse to save (cadence_batch_pipeline does, by default).
%
%   status
%     .messages   cellstr, one line per stage (the app's console lines)
%     .errors     cellstr of stage errors (empty = clean run)
%     .polarity   struct array per camera: cam, inverted, confidence, method
%                 (a Vm-Ca tie-breaker flip is recorded on that camera)
%     .qc         validate_conditioned records (also attached to d.qc)
%     .pacing     true when analog1 drove ensemble averaging
%     .timing     struct of seconds per stage
%     .settings   the resolved opts
%
%   See also cadence_batch_pipeline, cadence_convert_raw, cadence_extract_features.

    if nargin < 2 || isempty(opts), opts = struct(); end
    P = resolve_options(opts);
    logf = P.Log;
    name = char(P.Name);

    status = struct('messages', {{}}, 'errors', {{}}, 'polarity', struct('cam', {}, ...
        'inverted', {}, 'confidence', {}, 'method', {}), 'qc', [], 'pacing', false, ...
        'timing', struct(), 'settings', P);

    if ~isfield(d, 'num_files') || isempty(d.num_files)
        d.num_files = nnz(cellfun(@(f) ~isempty(regexp(f, '^CAM\d+$', 'once')), fieldnames(d)));
    end
    if isfield(d, 'conditioned') && isequal(d.conditioned, true)
        % The app refuses these ("re-conditioning over-smooths"); so do we.
        error('cadence_condition_data:alreadyConditioned', ...
            '%s is already conditioned; start from the raw/converted .mat', name);
    end

    % ---- 1) SNR maps (always; later stages depend on them) ----------------
    t0 = tic;
    try
        for i = 1:d.num_files
            d.(sprintf('CAM%d_SNR', i)) = extract_snr_mask(d.(sprintf('CAM%d', i)));
        end
    catch ME
        fail('snr', ME);
    end
    status.timing.snr = toc(t0);

    % ---- 2) drift correction ---------------------------------------------
    t0 = tic;
    try
        if P.Drift
            for i = 1:d.num_files
                f = sprintf('CAM%d', i);
                d.(f) = remove_Drift(d.(f));
            end
            say('Drift correction done to remove baseline wandering.');
        else
            say('Drift correction option was NOT selected for signal processing');
        end
    catch ME
        fail('drift', ME);
    end
    status.timing.drift = toc(t0);

    % ---- 3) SVD denoising (with binning fallback) -------------------------
    t0 = tic;
    try
        if P.SVD
            for i = 1:d.num_files
                f = sprintf('CAM%d', i);
                [den, info] = denoise_svd(d.(f), P.SVDRank);
                % Provenance: rank kept, or 0 when the camera fell back to
                % 5 x 5 binning. Field absent = SVD not selected.
                if info.applied
                    d.(f) = den;
                    d.(sprintf('CAM%d_svd_rank', i)) = info.K;
                    say(sprintf('SVD Denoising of data completed for: %s (rank %d)', f, info.K));
                else
                    d.(f) = binning(d.(f), 5);
                    d.(sprintf('CAM%d_svd_rank', i)) = 0;
                    say(['Spatial binning done using 5 x 5 box filter for:', f]);
                end
            end
            say('SVD Denoising done to attenuate noise.');
        else
            say('SVD Denoising option was NOT selected for signal processing');
        end
    catch ME
        fail('svd', ME);
    end
    status.timing.svd = toc(t0);

    % ---- 4) spatial binning ----------------------------------------------
    t0 = tic;
    try
        if P.Binning
            for i = 1:d.num_files
                f = sprintf('CAM%d', i);
                d.(f) = binning(d.(f), P.BinSize);
            end
            say(sprintf('Spatial binning done using %d x %d box filter', P.BinSize, P.BinSize));
        else
            say('Spatial binning was NOT selected for signal processing');
        end
    catch ME
        fail('binning', ME);
    end
    status.timing.binning = toc(t0);

    % ---- 5) temporal filtering -------------------------------------------
    t0 = tic;
    try
        if ~isempty(P.FilterHz) && P.FilterHz > 0
            fs = d.acqFreq;
            for i = 1:d.num_files
                f = sprintf('CAM%d', i);
                d.(f) = filter_data(d.(f), fs, P.FilterHz);
            end
            say(sprintf('Temporal filtering done using butterworth IIR filter [0, %d]', P.FilterHz));
        else
            say('Temporal filtering option was NOT selected for signal processing');
        end
    catch ME
        fail('filter', ME);
    end
    status.timing.filter = toc(t0);

    % ---- 6) motion artifact correction -----------------------------------
    t0 = tic;
    try
        if P.Motion
            if isfield(d, 'analog1') && ~isempty(d.analog1)
                [~, locs] = findpeaks(d.analog1);
            else
                % auto_detect_peaks returns a binary onset vector (one output)
                locs = find(auto_detect_peaks(d.CAM1));
            end
            if isempty(locs)
                say('Motion correction skipped: no beat found to choose a reference frame.');
            else
                ref_frame = max(1, locs(1) - 1);
                dm = d;       % commit only if every camera succeeds
                for i = 1:d.num_files
                    f = sprintf('CAM%d', i);
                    dm.(f) = trackCardiacMotion(dm.(f), 'RefFrame', ref_frame);
                end
                d = dm;
                say('Motion correction done to remove motion artifacts.');
            end
        else
            say('Motion correction option was NOT selected for signal processing');
        end
    catch ME
        fail('motion', ME);
    end
    status.timing.motion = toc(t0);

    % ---- 7) normalization ------------------------------------------------
    t0 = tic;
    try
        if P.Normalize
            for i = 1:d.num_files
                f = sprintf('CAM%d', i);
                d.(f) = normalize_data(d.(f));
            end
            say('Data Normalization completed.');
        else
            say('Data Normalization option was NOT selected for signal processing');
        end
    catch ME
        fail('normalize', ME);
    end
    status.timing.normalize = toc(t0);

    % ---- 8) signal inversion check (always) -------------------------------
    t0 = tic;
    try
        pacing = [];
        if isfield(d, 'analog1'), pacing = d.analog1; end
        pol_masks = cell(1, d.num_files);
        for i = 1:d.num_files
            cam_f = sprintf('CAM%d', i);
            avg_f = sprintf('CAM%d_average', i);
            snr_f = sprintf('CAM%d_SNR', i);
            snr_mask = create_snr_mask(d.(snr_f), 3, true);
            pol_masks{i} = snr_mask;

            % Stimulus-anchored first: after a stimulus a correctly oriented
            % signal must deflect UP.  Immune to duty cycle.
            [inv, conf] = check_polarity_paced(d.(cam_f), pacing, d.acqFreq, snr_mask);
            if isnan(conf) || abs(conf) < 0.20
                [inv, conf] = check_signal_inversion(d.(cam_f), d.acqFreq, 'Mask', snr_mask);
                % check_signal_inversion returns an UNSIGNED agreement in
                % [0, 1]; sign it with the decision so the saved value means
                % the same as check_polarity_paced's: + upright, - inverted,
                % magnitude = how decisive.
                conf = conf * (1 - 2*double(inv));
                method = 'signal shape';
                source = 'shape';
                why = sprintf('signal shape, confidence %+.2f (no usable pacing)', conf);
            else
                method = 'post-stimulus deflection';
                source = 'stimulus';
                why = sprintf('post-stimulus deflection, confidence %+.2f', conf);
            end
            % Session prior: the same camera's stimulus-based orientation in
            % this experiment outranks a shape-only decision.
            pr = prior_for(P.PolarityPrior, i);
            if strcmp(source, 'shape') && ~isempty(pr)
                if logical(pr.inverted) ~= logical(inv)
                    say(sprintf('Session prior OVERRULES the shape check for %s (shape said %s, %+.2f).', ...
                        cam_f, tern_str(inv, 'inverted', 'upright'), conf));
                end
                inv  = logical(pr.inverted);
                conf = (1 - 2*double(inv)) * pr.agree;
                method = sprintf('session prior (%d stimulus recordings)', pr.n);
                source = 'session';
                why = sprintf('session prior: %d paced recording(s) of this experiment, %.0f%% %s', ...
                    pr.n, 100*pr.agree, tern_str(inv, 'inverted', 'upright'));
            end
            % Provenance, saved with the file. Signed as measured on the data
            % BEFORE any flip: negative = the recording arrived inverted.
            d.(sprintf('CAM%d_polarity_confidence', i)) = conf;
            d.(sprintf('CAM%d_polarity_source', i))     = source;   % 'stimulus' | 'shape'

            if inv
                d.(cam_f) = 1 - d.(cam_f);
                if isfield(d, avg_f) && ~isempty(d.(avg_f))
                    d.(avg_f) = 1 - d.(avg_f);      % only when averaging already ran
                end
                say(sprintf('Inverted signals for: %s  [%s]', cam_f, why));
            else
                say(sprintf('Polarity OK for: %s  [%s]', cam_f, why));
            end
            status.polarity(end+1) = struct('cam', i, 'inverted', logical(inv), ...
                'confidence', conf, 'method', method); %#ok<AGROW>
        end
        % Vm-Ca timing tie-breaker for voltage cameras decided by shape alone.
        [d, status] = vm_ca_polarity(d, P, pol_masks, status, @say);
        say('Checked for signal inversion... Done');
    catch ME
        fail('inversion', ME);
    end
    status.timing.inversion = toc(t0);

    % ---- 9) ensemble averaging -------------------------------------------
    t0 = tic;
    try
        if P.Ensemble
            pacing = [];
            if isfield(d, 'analog1') && ~isempty(d.analog1)
                pacing = d.analog1;
            end
            for i = 1:d.num_files
                d.(sprintf('CAM%d_average', i)) = ensembleAverageFull(d.(sprintf('CAM%d', i)), pacing);
            end
            status.pacing = ~isempty(pacing);
            if ~isempty(P.LowBandHz) && P.LowBandHz > 0
                % hybrid filter: a low-band average from the same stack and the
                % same beat detection, BEFORE the final mask blanks any pixel
                for i = 1:d.num_files
                    f = sprintf('CAM%d', i);
                    d.([f '_average_lowband']) = ensembleAverageFull(filter_data(d.(f), d.acqFreq, P.LowBandHz), pacing);
                end
                say(sprintf('Low-band ensemble average [0, %g Hz] saved for the hybrid filter (CAM<n>_average_lowband).', P.LowBandHz));
            end
            if isempty(pacing)
                say('Ensemble averaging completed (peaks auto-detected from data).');
            else
                say('Ensemble averaging completed using pacing stimulus.');
            end
        else
            say('Ensemble averaging option was NOT selected for signal processing');
        end
    catch ME
        fail('ensemble', ME);
    end
    status.timing.ensemble = toc(t0);

    % Which temporal filter produced this file (read by cadence_extract_features).
    hz = P.FilterHz; if isempty(hz), hz = 0; end
    lb = P.LowBandHz; if isempty(lb) || ~P.Ensemble, lb = []; end
    d.temporal_filter = struct('conditioned_hz', hz, 'lowband_hz', lb);

    % ---- 10) SNR mask (always) --------------------------------------------
    t0 = tic;
    try
        for i = 1:d.num_files
            f = sprintf('CAM%d', i);
            snr_f = sprintf('CAM%d_SNR', i);
            snr_mask = double(create_snr_mask(d.(snr_f), P.MaskFloor, true));
            snr_mask(~snr_mask) = NaN;
            d.(f) = d.(f) .* snr_mask;
            say(['Data masking completed for: ', f]);
        end
        say('Data masking done using SNR mask');
    catch ME
        fail('mask', ME);
    end
    status.timing.mask = toc(t0);

    % ---- 11) provenance + QC (what the app does just before save) ---------
    t0 = tic;
    d.conditioned = true;   % provenance stamp -- guards against re-conditioning
    try
        vr = validate_conditioned(d, name);
        for r = vr(arrayfun(@(x) x.status ~= "PASS", vr))'
            say(sprintf('COND %s [%s]: %s', r.status, r.check, r.message));
        end
        if isfield(d, 'qc'), prior_qc = d.qc; else, prior_qc = []; end
        d.qc = qc_attach(prior_qc, vr);
        status.qc = vr;
    catch ME
        fail('qc', ME);
    end
    status.timing.qc = toc(t0);

    % ---- nested helpers ----------------------------------------------------
    function say(s)
        status.messages{end+1} = s;
        logf(s);
    end
    function fail(stage, ME)
        s = sprintf('ERROR - %s: %s', stage, ME.message);
        status.errors{end+1} = s;
        logf(s);
    end
end


function P = resolve_options(opts)
    P = struct('Drift', true, 'SVD', true, 'SVDRank', 8, 'Binning', false, 'BinSize', 3, ...
               'FilterHz', 50, 'LowBandHz', [], 'Motion', false, 'Normalize', true, 'Ensemble', true, ...
               'MaskFloor', 2, 'VoltageCams', [], 'CalciumCams', [], 'PolarityPrior', [], ...
               'Name', '', 'Log', @(s) fprintf('%s\n', s));
    fn = fieldnames(opts);
    for k = 1:numel(fn)
        if ~isfield(P, fn{k})
            error('cadence_condition_data:option', 'Unknown option ''%s''', fn{k});
        end
        P.(fn{k}) = opts.(fn{k});
    end
    if ~ismember(P.BinSize, [3 5 7 9])
        error('cadence_condition_data:option', 'BinSize must be 3, 5, 7 or 9 (got %g)', P.BinSize);
    end
    if isempty(P.MaskFloor), P.MaskFloor = 2; end
    if ~isempty(P.LowBandHz) && P.LowBandHz > 0
        if ~P.Ensemble
            error('cadence_condition_data:option', 'LowBandHz needs Ensemble = true (it stores a low-band ensemble average).');
        end
        if ~isempty(P.FilterHz) && P.FilterHz > 0 && P.LowBandHz >= P.FilterHz
            error('cadence_condition_data:option', 'LowBandHz (%g) must be below FilterHz (%g).', P.LowBandHz, P.FilterHz);
        end
    end
end


function [d, status] = vm_ca_polarity(d, P, masks, status, say)
% See VoltageCams / CalciumCams in the header.  Runs after the per-camera
% checks (so calcium is already oriented) and before ensemble averaging.
    ca = P.CalciumCams(P.CalciumCams <= d.num_files);
    vs = P.VoltageCams(P.VoltageCams <= d.num_files);
    if isempty(ca) || isempty(vs), return; end
    ca = ca(1);
    ca_conf = abs(get_num(d, sprintf('CAM%d_polarity_confidence', ca)));
    ca_src  = get_str(d, sprintf('CAM%d_polarity_source', ca));
    for v = vs(:)'
        if ~strcmp(get_str(d, sprintf('CAM%d_polarity_source', v)), 'shape'), continue; end
        v_conf = abs(get_num(d, sprintf('CAM%d_polarity_confidence', v)));
        if ~strcmp(ca_src, 'stimulus') && ca_conf < v_conf
            say(sprintf('Vm-Ca polarity check skipped for CAM%d: calcium (CAM%d, %s %.2f) is less decisive than voltage (%.2f).', ...
                v, ca, ca_src, ca_conf, v_conf));
            continue;
        end
        [keep_up, flip, info] = lag_test(d.(sprintf('CAM%d', v)), masks{v}, ...
                                         d.(sprintf('CAM%d', ca)), masks{ca}, d.acqFreq);
        % saved confidence keeps the convention: negative = arrived inverted
        was_inv = get_num(d, sprintf('CAM%d_polarity_confidence', v)) < 0;
        margin  = abs(info.score_keep - info.score_flip);
        if flip
            d.(sprintf('CAM%d', v)) = 1 - d.(sprintf('CAM%d', v));
            avg_f = sprintf('CAM%d_average', v);
            if isfield(d, avg_f) && ~isempty(d.(avg_f)), d.(avg_f) = 1 - d.(avg_f); end
            d.(sprintf('CAM%d_polarity_source', v)) = 'vm-ca lag';
            d.(sprintf('CAM%d_polarity_confidence', v)) = (1 - 2*double(~was_inv)) * margin;
            k = find([status.polarity.cam] == v, 1, 'last');
            if ~isempty(k)
                status.polarity(k).inverted = ~status.polarity(k).inverted;
                status.polarity(k).method = 'Vm-Ca lag';
            end
            say(sprintf(['Vm-Ca lag check FLIPPED CAM%d: calcium follows the inverted voltage at %+.1f ms ' ...
                '(r %.2f) vs %+.1f ms as it was (r %.2f).'], v, info.lag_flip_ms, info.score_flip, info.lag_keep_ms, info.score_keep));
        elseif keep_up
            d.(sprintf('CAM%d_polarity_source', v)) = 'vm-ca lag';
            d.(sprintf('CAM%d_polarity_confidence', v)) = (1 - 2*double(was_inv)) * margin;
            say(sprintf('Vm-Ca lag check confirms CAM%d: calcium follows at %+.1f ms (r %.2f; flipped %.2f).', ...
                v, info.lag_keep_ms, info.score_keep, info.score_flip));
        else
            say(sprintf('Vm-Ca lag check inconclusive for CAM%d (%s: r %.2f as is, %.2f flipped); shape decision kept.', ...
                v, info.mode, info.score_keep, info.score_flip));
        end
    end
end

function [keep_up, flip, info] = lag_test(V, mv, C, mc, fs)
% Correlate voltage with calcium shifted by lags in [-2 ms, +min(15 ms, 0.3
% cycle)], for the voltage as is and inverted.  Calcium trails voltage by a
% few ms; an inverted voltage puts the best lag half a cycle away.
%   per-pixel (cameras on one pixel grid): median over pixels in both masks of
%     each pixel's correlation -- robust to a rotating wave, whose tissue mean
%     cancels.  Decisive when the winner is >= 0.3 and leads by >= 0.2.
%   tissue mean (grids differ): weaker; decisive only at >= 0.5 and a 0.2 lead.
    keep_up = false;  flip = false;
    info = struct('mode', '', 'score_keep', NaN, 'score_flip', NaN, 'lag_keep_ms', NaN, 'lag_flip_ms', NaN);
    keepv = mv(:) > 0 & isfinite(mv(:));
    keepc = mc(:) > 0 & isfinite(mc(:));
    T = size(V, 3);
    if isequal(size(V), size(C))
        both = find(keepv & keepc);
        if numel(both) > 3000, both = both(round(linspace(1, numel(both), 3000))); end
        X = reshape(V, [], T);  Y = reshape(C, [], T);
        X = X(both, :);  Y = Y(both, :);
        ok = all(isfinite(X), 2) & all(isfinite(Y), 2);
        X = X(ok, :);  Y = Y(ok, :);
        info.mode = 'per-pixel';  min_r = 0.3;
    else
        X = masked_mean(V, mv)';  Y = masked_mean(C, mc)';
        ok = isfinite(X) & isfinite(Y);  X = X(ok);  Y = Y(ok);
        info.mode = 'tissue mean';  min_r = 0.5;
    end
    if size(X, 1) < 1 || size(X, 2) < 50, return; end
    X = detrend(X')';  Y = detrend(Y')';
    df = cardiacSpectralMetrics(median(X, 1)', fs);
    hi = 0.015;
    if isfinite(df) && df > 0, hi = min(hi, 0.3 / df); end
    lags = round(-0.002 * fs):max(1, round(hi * fs));
    r = nan(size(lags));
    for k = 1:numel(lags)
        L = lags(k);
        if L >= 0, a = X(:, 1:end-L); b = Y(:, 1+L:end);
        else,      a = X(:, 1-L:end); b = Y(:, 1:end+L); end
        a = a - mean(a, 2);  b = b - mean(b, 2);
        rp = sum(a .* b, 2) ./ sqrt(sum(a.^2, 2) .* sum(b.^2, 2));
        r(k) = median(rp, 'omitnan');
    end
    [info.score_keep, ik] = max(r);
    [info.score_flip, iff] = max(-r);
    info.lag_keep_ms = 1000 * lags(ik) / fs;
    info.lag_flip_ms = 1000 * lags(iff) / fs;
    margin = 0.2;
    flip    = info.score_flip >= min_r && info.score_flip - info.score_keep >= margin;
    keep_up = info.score_keep >= min_r && info.score_keep - info.score_flip >= margin;
end

function tr = masked_mean(X, mask)
    keep = mask(:) > 0 & isfinite(mask(:));
    F = reshape(X, size(X,1) * size(X,2), []);
    tr = mean(F(keep, :), 1, 'omitnan')';
end

function v = get_num(d, f)
    if isfield(d, f) && ~isempty(d.(f)), v = double(d.(f)); else, v = NaN; end
end

function s = get_str(d, f)
    if isfield(d, f) && ~isempty(d.(f)), s = char(d.(f)); else, s = ''; end
end


function pr = prior_for(prior, cam)
    pr = [];
    if isempty(prior) || ~isstruct(prior) || ~isfield(prior, 'cam'), return; end
    k = find([prior.cam] == cam, 1);
    if ~isempty(k), pr = prior(k); end
end

function s = tern_str(c, a, b)
    if c, s = a; else, s = b; end
end
