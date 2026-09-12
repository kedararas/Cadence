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
%     Motion       false     G) MOTION ARTIFACT CORRECTION
%     Normalize    true      E) DATA NORMALIZATION
%     Ensemble     true      F) ENSEMBLE AVERAGING
%     MaskFloor    2         SNR floor of the final tissue mask (app: fixed 2)
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
                if info.applied
                    d.(f) = den;
                    say(['SVD Denoising of data completed for: ', f]);
                else
                    d.(f) = binning(d.(f), 5);
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
                [~, locs] = auto_detect_peaks(d.CAM1);
            end
            ref_frame = locs(1,1) - 1;
            for i = 1:d.num_files
                f = sprintf('CAM%d', i);
                d.(f) = trackCardiacMotion(d.(f), 'RefFrame', ref_frame);
            end
            say('Motion correction done to remove motion artifacts.');
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
        for i = 1:d.num_files
            cam_f = sprintf('CAM%d', i);
            avg_f = sprintf('CAM%d_average', i);
            snr_f = sprintf('CAM%d_SNR', i);
            snr_mask = create_snr_mask(d.(snr_f), 3, true);

            % Stimulus-anchored first: after a stimulus a correctly oriented
            % signal must deflect UP.  Immune to duty cycle.
            [inv, conf] = check_polarity_paced(d.(cam_f), pacing, d.acqFreq, snr_mask);
            if isnan(conf) || abs(conf) < 0.20
                [inv, conf] = check_signal_inversion(d.(cam_f), d.acqFreq, 'Mask', snr_mask);
                method = 'signal shape';
                why = sprintf('signal shape, confidence %.2f (no usable pacing)', conf);
            else
                method = 'post-stimulus deflection';
                why = sprintf('post-stimulus deflection, confidence %.2f', conf);
            end
            d.(sprintf('CAM%d_polarity_confidence', i)) = conf;   % provenance, saved with the file

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
               'FilterHz', 50, 'Motion', false, 'Normalize', true, 'Ensemble', true, ...
               'MaskFloor', 2, 'Name', '', 'Log', @(s) fprintf('%s\n', s));
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
end
