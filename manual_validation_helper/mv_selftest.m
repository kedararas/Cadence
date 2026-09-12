function mv_selftest(keep_dir)
%MV_SELFTEST  End-to-end check of the marking harness on synthetic data.
%
%   mv_selftest
%   mv_selftest(true)      % keep the scratch directory for inspection
%
%   Builds a synthetic conditioned recording in which the true APD80 of every
%   pixel is known by construction, simulates a set of reviewers marking it,
%   and runs the whole non-interactive path:
%
%       mv_sample_pixels -> (simulated marking) -> mv_derive -> mv_compare
%
%   It then asserts that mv_compare recovers a bias deliberately injected into
%   the "CADENCE" map.  Run this before recruiting reviewers: it verifies the
%   harness end to end without anyone clicking anything.
%
%   Only mv_mark's interactive loop is not exercised (it needs a human).  The
%   simulated marks are written in exactly the format mv_mark produces, so if
%   this passes, everything downstream of marking is known good.
%
%   This function is also a worked example of the synthetic ground-truth idea:
%   the same AP model, propagated at a known velocity, is what you would use to
%   validate conduction velocity, where manual pixel-by-pixel marking is not
%   possible.

    if nargin < 1; keep_dir = false; end

    root = fullfile(tempdir, sprintf('mv_selftest_%s', ...
                    char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
    mkdir(root);
    fprintf('Scratch: %s\n\n', root);

    % ================= 1. synthetic recording with known truth =================
    fs   = 1000;            % Hz
    nr   = 32; nc = 32;
    nt   = 2000;            % 2 s
    cl   = 200;             % pacing cycle length, ms
    amp  = 1;

    % Truth, varying in both axes so the comparison is not degenerate:
    %   activation sweeps across columns (a planar wave)
    %   APD80 varies down rows
    [RR, CC]   = ndgrid(1:nr, 1:nc);
    delay_ms   = 0.5 * (CC - 1);               % 0 .. 15.5 ms across the tissue
    true_apd   = 100 + 0.5 * (RR - 1);         % 100 .. 115.5 ms

    % Noise in three row-bands so the SNR strata have something to separate.
    sigma = 0.02 * ones(nr, nc);
    sigma(RR > nr/3)     = 0.06;
    sigma(RR > 2*nr/3)   = 0.14;

    t_ms  = (0:nt-1) / fs * 1000;
    stim  = 0:cl:(nt/fs*1000 - cl);            % stimulus times, ms
    X     = zeros(nr, nc, nt);
    rng(7, 'twister');

    for s = 1:numel(stim)
        for r = 1:nr
            for c = 1:nc
                ta = stim(s) + delay_ms(r, c);
                X(r, c, :) = squeeze(X(r, c, :))' + ap_waveform(t_ms, ta, true_apd(r, c), amp);
            end
        end
    end
    X = X + sigma .* randn(nr, nc, nt);

    % Pacing channel: 2 ms pulses at each stimulus.
    analog1 = zeros(1, nt);
    for s = 1:numel(stim)
        i0 = round(stim(s) / 1000 * fs) + 1;
        analog1(i0:min(nt, i0 + 1)) = 1;
    end

    % The ensemble-averaged beat, as Signal Conditioning would leave it: one
    % representative beat starting a pre-window before onset.  This is the array
    % feature extraction actually uses under the recommended workflow, and the
    % one reviewers should therefore mark.
    T_avg     = round(cl * 1.1);
    onset_avg = round(0.1 * cl);
    t_avg     = (0:T_avg-1) / fs * 1000;
    Xavg      = zeros(nr, nc, T_avg);
    for r = 1:nr
        for c = 1:nc
            Xavg(r, c, :) = ap_waveform(t_avg, onset_avg + delay_ms(r, c), true_apd(r, c), amp);
        end
    end

    cmos_all_data = struct();
    cmos_all_data.CAM1         = X;
    cmos_all_data.CAM1_average = Xavg;
    cmos_all_data.CAM1_SNR     = amp ./ sigma;  % the SNR convention: >0 is tissue
    cmos_all_data.acqFreq      = fs;
    cmos_all_data.analog1      = analog1;

    cond_file = fullfile(root, 'synthetic-conditioned.mat');
    save(cond_file, 'cmos_all_data');
    fprintf('1. Synthetic recording: %dx%dx%d, true APD80 %.1f-%.1f ms, SNR %.1f-%.1f\n\n', ...
            nr, nc, nt, min(true_apd(:)), max(true_apd(:)), ...
            min(cmos_all_data.CAM1_SNR(:)), max(cmos_all_data.CAM1_SNR(:)));

    % ============================ 2. sample pixels ============================
    fprintf('2. mv_sample_pixels\n');
    manifest = mv_sample_pixels(cond_file, 'CAM1', 90, 'BeatIndex', 5, 'Seed', 11, ...
                                'Source', 'raw', 'Candidates', 1, 'SNRFloor', 0);
    mf = manifest.manifest_file;
    fprintf('\n');

    % The map panels need pixels laid out so they can be interpolated, so the
    % grid path has to work as well as the stratified one.  Checked at step 6.
    fprintf('2b. mv_sample_pixels (grid)\n');
    gman = mv_sample_pixels(cond_file, 'CAM1', 90, 'BeatIndex', 5, 'Seed', 11, ...
                            'Sampling', 'grid', 'Source', 'raw', 'Candidates', 1, 'SNRFloor', 0, ...
                            'Output', fullfile(root, 'grid_manifest.mat'));
    fprintf('\n');

    % ====================== 3. simulate three reviewers ========================
    % Each reviewer marks the true fiducials plus independent jitter, which is
    % what inter-observer spread is.  Jitter grows as SNR falls, exactly as a
    % real reviewer's precision would.
    fprintf('3. Simulating 3 reviewers\n');
    beat_t   = stim(5);                                  % the beat in the window
    reviewers = {'R1', 'R2', 'R3'};
    marks_files = cell(1, numel(reviewers));

    for q = 1:numel(reviewers)
        rng(100 + q, 'twister');
        marks_files{q} = simulate_reviewer(mf, reviewers{q}, beat_t, ...
                                           delay_ms, true_apd, cmos_all_data.CAM1_SNR, root);
    end
    fprintf('\n');

    % ============================== 4. derive =================================
    fprintf('4. mv_derive\n');
    tbl = mv_derive(marks_files);
    check('simulated reviewers leave no slips', height(tbl.Properties.UserData.dropped) == 0, ...
          sprintf('%d dropped', height(tbl.Properties.UserData.dropped)));

    % A downstroke activation click (after the peak) gives a short positive
    % APD that the APD <= 0 rule misses.  It must be dropped and listed.
    S = load(marks_files{1}); marks = S.marks;
    slip = find(marks.done, 1);
    marks.t_act_ms(slip) = marks.t_peak_ms(slip) + 5;
    slip_file = fullfile(root, 'sim_slip_R1_marks.mat');
    save(slip_file, 'marks');
    stbl = mv_derive(slip_file);
    dr   = stbl.Properties.UserData.dropped;
    check('click-order slip is dropped and listed', ...
          ~any(stbl.pixel == slip) && height(dr) == 1 && dr.pixel(1) == slip && dr.reason(1) == "click order", ...
          sprintf('pixel %d, %d dropped', slip, height(dr)));
    fprintf('\n');

    % ================= 5. a "CADENCE" map with a known bias ===================
    % Inject a deliberate systematic offset; mv_compare must recover it.
    injected_bias = 1.5;                                   % ms
    apd_map = true_apd + injected_bias;

    m2 = struct();
    m2.CAM1       = X;
    m2.acqFreq    = fs;
    m2.num_files  = 1;
    m2.analog1    = analog1;
    m2.ep_metrics = struct('apd_data', {{apd_map}});
    cmos_all_data = m2;
    metrics_file  = fullfile(root, 'synthetic-metrics.mat');
    save(metrics_file, 'cmos_all_data');

    fprintf('5. mv_compare (injected bias %+.2f ms)\n', injected_bias);
    stats = mv_compare(tbl, metrics_file, 'Field', 'apd_data', ...
                       'ManualVar', 'apd_ms', 'Tolerance', [5 10], 'Plot', false);

    % ============================== 6. assert =================================
    fprintf('\n6. Checks\n');
    check('bias recovers the injected offset', abs(stats.bias - injected_bias) < 1.0, ...
          sprintf('bias %+.2f vs injected %+.2f', stats.bias, injected_bias));
    check('correlation is high',               stats.r > 0.8, sprintf('r = %.3f', stats.r));
    check('inter-observer was computed',       stats.interobserver.n_pairs > 0, ...
          sprintf('%d human-human pairs', stats.interobserver.n_pairs));
    % Stratum 1 is the LOWEST SNR band, so its spread should be the largest.
    % This is the operating envelope the stratified sample exists to expose.
    check('agreement degrades as SNR falls', ...
          stats.by_stratum(1).sd > stats.by_stratum(end).sd, ...
          sprintf('SD %.2f at SNR %.1f vs %.2f at SNR %.1f', ...
                  stats.by_stratum(1).sd,   stats.by_stratum(1).snr_median, ...
                  stats.by_stratum(end).sd, stats.by_stratum(end).snr_median));

    % Grid sampling.  The count is checked loosely on purpose: a square lattice
    % can only yield certain counts, so n_pixels is a target and hitting it
    % exactly is neither possible nor wanted.
    ng = size(gman.pixels, 1);
    check('grid sampling lands near its target', ng >= 0.6 * 90 && ng <= 1.5 * 90, ...
          sprintf('%d pixels at stride %d (target 90)', ng, gman.stride));
    check('grid pixels lie on one lattice', ...
          isscalar(unique(mod(gman.pixels(:,1), gman.stride))) && ...
          isscalar(unique(mod(gman.pixels(:,2), gman.stride))), ...
          sprintf('stride %d in both axes', gman.stride));
    check('grid strata are all populated', ...
          isequal(sort(unique(gman.stratum))', 1:gman.num_strata), ...
          sprintf('%d of %d strata used', numel(unique(gman.stratum)), gman.num_strata));

    % ============ 7. window-relative TIME: origin correction ==================
    % apd_data is a duration, so the window origin cancels and everything above
    % works without knowing where the window starts.  act_times does NOT: it is
    % measured in frames from the start of CADENCE's own analysis window, while
    % a reviewer clicks in absolute recording time.  Uncorrected, the bias is
    % the whole window offset and Bland-Altman still looks perfect.  This checks
    % that mv_compare shifts it back, and refuses when it cannot.
    fprintf('\n7. Window-relative time (act_times)\n');
    MMf = load(mf); mani = MMf.manifest;
    sf  = mani.window_frames(1);
    ef  = mani.window_frames(2);
    dt  = 1000 / fs;

    % CADENCE's stored value: 1-based fractional frame index inside its window.
    F_abs   = (beat_t + delay_ms) / dt + 1;          % absolute fractional frame
    act_map = (F_abs - sf + 1) * dt;                 % what compute_lat_50*dt gives

    m3 = struct();
    m3.CAM1       = X;
    m3.acqFreq    = fs;
    m3.num_files  = 1;
    m3.analog1    = analog1;
    m3.window     = {[sf ef]};                        % as compile_data saves it
    m3.ep_metrics = struct('act_times', {{act_map}});
    cmos_all_data = m3;
    time_file = fullfile(root, 'synthetic-time-metrics.mat');
    save(time_file, 'cmos_all_data');

    st = mv_compare(tbl, time_file, 'Field', 'act_times', ...
                    'ManualVar', 'act_ms', 'Tolerance', [5 10], 'Plot', false);

    % Same file with the window stripped, i.e. what an ensemble-average run
    % leaves behind: the correction is undefinable and must be refused, not
    % guessed.
    cmos_all_data = rmfield(m3, 'window');
    nowin_file = fullfile(root, 'synthetic-nowindow-metrics.mat');
    save(nowin_file, 'cmos_all_data');
    refused_nowin = false;
    try
        mv_compare(tbl, nowin_file, 'Field', 'act_times', 'ManualVar', 'act_ms', 'Plot', false);
    catch ME
        refused_nowin = strcmp(ME.identifier, 'mv_compare:noWindow');
    end

    % And a window on a different beat entirely.
    m4 = m3; m4.window = {[ef + 500, ef + 700]};
    cmos_all_data = m4;
    otherbeat_file = fullfile(root, 'synthetic-otherbeat-metrics.mat');
    save(otherbeat_file, 'cmos_all_data');
    refused_beat = false;
    try
        mv_compare(tbl, otherbeat_file, 'Field', 'act_times', 'ManualVar', 'act_ms', 'Plot', false);
    catch ME
        refused_beat = strcmp(ME.identifier, 'mv_compare:differentBeat');
    end

    % A duration compared against an absolute time must also be refused.
    refused_kind = false;
    try
        mv_compare(tbl, metrics_file, 'Field', 'apd_data', 'ManualVar', 'act_ms', 'Plot', false);
    catch ME
        refused_kind = strcmp(ME.identifier, 'mv_compare:kindMismatch');
    end

    fprintf('\n8. Checks (origin handling)\n');
    check('origin correction removes the window offset', abs(st.bias) < 1.0, ...
          sprintf('bias %+.2f ms after shifting by %+.1f ms', st.bias, (sf-2)*dt));
    check('uncorrected bias would have been large', abs((sf - 2) * dt) > 5, ...
          sprintf('window starts at frame %d -> %+.1f ms', sf, (sf-2)*dt));
    check('refuses when no analysis window is recorded', refused_nowin, 'mv_compare:noWindow');
    check('refuses when the beats differ',               refused_beat,  'mv_compare:differentBeat');
    check('refuses a duration vs an absolute time',      refused_kind,  'mv_compare:kindMismatch');

    % ============ 9. ENSEMBLE arm — the recommended workflow ==================
    % Feature extraction uses CAM<n>_average whenever the ensemble box is
    % ticked, and every metric except alternans derives from it.  So reviewers
    % mark that same array: both sides then start at frame 1, the only residual
    % is the 1-based index convention, and act_times becomes validatable with no
    % window, no beat matching and no recording timeline involved at all.
    fprintf('\n9. Ensemble-averaged arm\n');
    snr_map = amp ./ sigma;
    eman = mv_sample_pixels(cond_file, 'CAM1', 90, 'Seed', 11, 'Source', 'ensemble', ...
                            'Candidates', 1, 'SNRFloor', 0, ...
                            'Output', fullfile(root, 'ens_manifest.mat'));
    emf = eman.manifest_file;

    emarks = cell(1, numel(reviewers));
    for q = 1:numel(reviewers)
        rng(200 + q, 'twister');
        emarks{q} = simulate_reviewer(emf, reviewers{q}, onset_avg, ...
                                      delay_ms, true_apd, snr_map, root, 'ens_');
    end
    etbl = mv_derive(emarks);

    dtm     = 1000 / fs;
    ta_avg  = onset_avg + delay_ms;
    act_ens = (ta_avg / dtm + 1) * dtm;     % compute_lat_50's 1-based index * dt

    m5 = struct();
    m5.CAM1       = X;
    m5.acqFreq    = fs;
    m5.num_files  = 1;
    m5.analog1    = analog1;
    m5.ep_metrics = struct('act_times', {{act_ens}});
    % Deliberately NO 'window' field: an ensemble run records none, and under
    % this arm none is needed.
    cmos_all_data = m5;
    ens_file = fullfile(root, 'synthetic-ensemble-metrics.mat');
    save(ens_file, 'cmos_all_data');

    est = mv_compare(etbl, ens_file, 'Field', 'act_times', 'ManualVar', 'act_ms', ...
                     'Tolerance', [5 10], 'Plot', false);

    % Alternans cannot come from an averaged beat — averaging is precisely what
    % removes the alternation.
    refused_alt = false;
    try
        mv_sample_pixels(cond_file, 'CAM1', 20, 'Source', 'ensemble', 'NumBeats', 2, ...
                         'Output', fullfile(root, 'bad_alt.mat'));
    catch ME
        refused_alt = strcmp(ME.identifier, 'mv_sample_pixels:ensembleAlternans');
    end

    % And a conditioned file with no average must say so plainly.
    cmos_all_data = rmfield(m5, 'ep_metrics');
    cmos_all_data.CAM1_SNR = snr_map;
    noavg_file = fullfile(root, 'synthetic-noavg-conditioned.mat');
    save(noavg_file, 'cmos_all_data');
    refused_noavg = false;
    try
        mv_sample_pixels(noavg_file, 'CAM1', 20, 'Source', 'ensemble', ...
                         'Output', fullfile(root, 'bad_noavg.mat'));
    catch ME
        refused_noavg = strcmp(ME.identifier, 'mv_sample_pixels:noEnsemble');
    end

    fprintf('\n10. Checks (ensemble arm)\n');
    check('reviewers marked the averaged array', ...
          strcmp(eman.data_field, 'CAM1_average') && strcmp(eman.source, 'ensemble'), ...
          sprintf('data_field = %s', eman.data_field));
    check('act_times validates with no window at all', abs(est.bias) < 1.0, ...
          sprintf('bias %+.2f ms (shared frame-1 origin)', est.bias));
    check('refuses alternans on the ensemble average', refused_alt, ...
          'mv_sample_pixels:ensembleAlternans');
    check('refuses when no ensemble average exists',   refused_noavg, ...
          'mv_sample_pixels:noEnsemble');

    % ============ 11. candidate pool built by a lead reviewer ================
    % The draw oversamples; the lead accepts or skips in manifest order until
    % the target is reached; the accepted set becomes the pool everyone else
    % marks.  Also exercises the adaptive SNR floor: this synthetic heart has
    % three SNR bands (50, 16.7, 7.1), so the pipeline's rule — half the median
    % tissue SNR, never below 1.5 — lands at 8.3 and excludes the dimmest band,
    % exactly the pixels Feature Extraction would mask out.
    fprintf('\n11. Lead-marker pool\n');
    target = 60;
    pman = mv_sample_pixels(cond_file, 'CAM1', target, 'Seed', 11, 'Source', 'ensemble', ...
                            'Candidates', 2, 'NumStrata', 2, ...
                            'Output', fullfile(root, 'pool_manifest.mat'));
    pmf = pman.manifest_file;
    exp_floor = max(1.5, 0.5 * median(snr_map(isfinite(snr_map) & snr_map > 1.5)));
    check('adaptive SNR floor matches the pipeline rule', ...
          abs(pman.snr_floor - exp_floor) < 1e-9 && all(pman.snr >= pman.snr_floor), ...
          sprintf('floor %.2f, min candidate SNR %.1f', pman.snr_floor, min(pman.snr)));
    ncand = size(pman.pixels, 1);
    check('pool starts open with 2x candidates', ...
          strcmp(pman.pool_status, 'open') && ncand == 2 * target && ~any(pman.active), ...
          sprintf('%d candidates for target %d', ncand, pman.n_target));

    % A non-lead reviewer's marks cannot enter an open pool.
    rng(301, 'twister');
    early = simulate_reviewer(pmf, 'R2', onset_avg, delay_ms, true_apd, snr_map, root, 'pool_early_');
    refused_open = false;
    try
        mv_derive(early);
    catch ME
        refused_open = strcmp(ME.identifier, 'mv_derive:poolOpen');
    end
    check('refuses to derive from an open pool', refused_open, 'mv_derive:poolOpen');

    % Lead: every 4th candidate is unmarkable; the session stops at the target.
    skip = mod((1:ncand)', 4) == 0;
    rng(300, 'twister');
    lead_file = simulate_reviewer(pmf, 'LEAD', onset_avg, delay_ms, true_apd, snr_map, root, ...
                                  'pool_', 'Skip', skip, 'Stop', target);
    pman = mv_finalize_pool(pmf, 'Marks', lead_file);
    % Expected pool: per stratum, the first quota markable candidates in
    % manifest order.  Every 4th candidate falls in the same stratum under
    % 2-way round-robin, so a global stop rule would have left that stratum
    % short — the per-stratum quota is what keeps the pool balanced.
    acc   = find(pman.active);
    quota = target / pman.num_strata;
    exp_acc = [];
    for q = 1:pman.num_strata
        cq = find(~skip & pman.stratum == q);
        exp_acc = [exp_acc; cq(1:quota)]; %#ok<AGROW>
    end
    check('pool is the first markable candidates per stratum', ...
          isequal(acc, sort(exp_acc)) && numel(acc) == target, ...
          sprintf('%d active, %d candidates visited', numel(acc), pman.lead_visited));
    check('skipped candidates are excluded and recorded', ...
          ~any(pman.active & skip) && isequal(pman.lead_skipped, find(skip & (1:ncand)' <= pman.lead_visited)), ...
          sprintf('%d skipped = %.0f%% of visited', numel(pman.lead_skipped), ...
                  100 * numel(pman.lead_skipped) / pman.lead_visited));
    cnt = accumarray(pman.stratum(acc), 1);
    check('strata stay balanced at the stopping point', max(cnt) - min(cnt) <= 1, mat2str(cnt'));

    % Everyone else marks the pool only.
    others = cell(1, 2);
    for q = 1:2
        rng(310 + q, 'twister');
        others{q} = simulate_reviewer(pmf, sprintf('R%d', q + 1), onset_avg, delay_ms, true_apd, ...
                                      snr_map, root, 'pool_', 'ActiveOnly', true);
    end
    ptbl = mv_derive([{lead_file}, others]);
    check('derived rows cover the pool only', ...
          height(ptbl) == 3 * target && numel(unique(ptbl.pixel)) == target && all(pman.active(ptbl.pixel)), ...
          sprintf('%d rows, %d pixels, 3 reviewers', height(ptbl), numel(unique(ptbl.pixel))));
    refused_refinal = false;
    try
        mv_finalize_pool(pmf, 'Marks', lead_file);
    catch ME
        refused_refinal = strcmp(ME.identifier, 'mv_finalize_pool:alreadyFinal');
    end
    check('a final pool is not re-opened', refused_refinal, 'mv_finalize_pool:alreadyFinal');

    % Grid: tiered candidates, tier 1 filled first, strata assigned on the pool.
    gp = mv_sample_pixels(cond_file, 'CAM1', target, 'Seed', 11, 'Source', 'ensemble', ...
                          'Sampling', 'grid', 'Candidates', 2, ...
                          'Output', fullfile(root, 'pool_grid_manifest.mat'));
    n1 = sum(gp.tier == 1);
    check('grid candidates come in tiers, strata pending', ...
          n1 >= 0.6 * target && n1 <= 1.5 * target && any(gp.tier == 2) && all(isnan(gp.stratum)), ...
          sprintf('tier 1 %d, later tiers %d', n1, sum(gp.tier > 1)));
    ncg   = size(gp.pixels, 1);
    skipg = mod((1:ncg)', 5) == 0;
    rng(320, 'twister');
    glead = simulate_reviewer(gp.manifest_file, 'LEAD', onset_avg, delay_ms, true_apd, snr_map, ...
                              root, 'pool_grid_', 'Skip', skipg, 'Stop', target);
    gp = mv_finalize_pool(gp.manifest_file, 'Marks', glead);
    ga = find(gp.active);
    t1_ok = gp.tier == 1 & ~skipg;
    if sum(t1_ok) >= target
        fill_ok = ~any(gp.active & gp.tier > 1);
    else
        fill_ok = all(gp.active(t1_ok));
    end
    check('grid pool fills tier 1 before later tiers', fill_ok, ...
          sprintf('%d from tier 1, %d from later tiers', ...
                  sum(gp.active & gp.tier == 1), sum(gp.active & gp.tier > 1)));
    check('grid strata assigned over the pool only', ...
          isequal(sort(unique(gp.stratum(ga)))', 1:gp.num_strata) && all(isnan(gp.stratum(~gp.active))), ...
          sprintf('%d strata over %d pixels', numel(unique(gp.stratum(ga))), numel(ga)));

    fprintf('\nAll checks passed.\n');
    if keep_dir
        fprintf('Scratch kept: %s\n', root);
    else
        rmdir(root, 's');
    end
end


% ======================================================================
function y = ap_waveform(t_ms, ta, apd80, amp)
%AP_WAVEFORM  One action potential whose APD80 is exactly apd80 by construction.
%
%   Rises linearly over 2 ms from ta, then decays exponentially with tau chosen
%   so the trace sits at 20% of peak (i.e. 80% repolarized) exactly apd80 ms
%   after ta.  That makes the ground truth analytic rather than measured.

    t_up = ta + 2;
    tau  = (apd80 - 2) / log(5);            % exp(-(apd80-2)/tau) = 0.2
    y    = zeros(size(t_ms));

    rise = t_ms >= ta & t_ms < t_up;
    y(rise) = amp * (t_ms(rise) - ta) / 2;

    dec = t_ms >= t_up;
    y(dec) = amp * exp(-(t_ms(dec) - t_up) / tau);
end


function f = simulate_reviewer(manifest_file, reviewer, beat_t, delay_ms, true_apd, snr_map, root, tag, varargin)
%SIMULATE_REVIEWER  Write a marks file in mv_mark's format, without a human.
%
%   Name-value (after tag)
%     'Skip'        logical over candidates: the reviewer finds these unmarkable
%     'Stop'        stop once this many pixels are accepted (a lead session)
%     'ActiveOnly'  mark only manifest.active pixels (an ordinary reviewer on
%                   a finalised pool)
    if nargin < 8 || isempty(tag); tag = ''; end
    q = inputParser;
    q.addParameter('Skip',       [],    @(v) isempty(v) || islogical(v));
    q.addParameter('Stop',       Inf,   @(v) isnumeric(v) && isscalar(v));
    q.addParameter('ActiveOnly', false, @islogical);
    q.parse(varargin{:});
    so = q.Results;

    M = load(manifest_file); manifest = M.manifest;
    n = size(manifest.pixels, 1);
    if isempty(so.Skip); so.Skip = false(n, 1); end
    if so.ActiveOnly; present = logical(manifest.active(:)); else; present = true(n, 1); end
    % A lead session on a stratified draw fills a per-stratum quota (mv_mark).
    if isfinite(so.Stop) && strcmp(manifest.sampling, 'stratified')
        quota = floor(so.Stop / manifest.num_strata);
        str   = manifest.stratum(:);
    else
        quota = Inf; str = ones(n, 1);
    end

    marks = struct();
    marks.manifest_file    = manifest_file;
    marks.conditioned_file = manifest.conditioned_file;
    marks.cam              = manifest.cam;
    marks.reviewer         = reviewer;
    marks.mode             = 'apd';
    marks.percent          = 80;
    marks.act_percent      = 50;
    marks.baseline_win_ms  = 5;
    marks.order_seed       = 0;
    marks.lead             = isfinite(so.Stop);
    marks.order            = (1:n)';
    marks.done             = false(n, 1);
    marks.skipped          = false(n, 1);
    marks.clicks           = repmat({zeros(0, 2)}, n, 1);
    marks.seconds          = 8 + 3 * rand(n, 1);
    marks.t_base_ms        = nan(n, 1);  marks.v_base    = nan(n, 1);
    marks.t_peak_ms        = nan(n, 1);  marks.v_peak    = nan(n, 1);
    marks.t_act_ms         = nan(n, 1);  marks.t_rep_ms  = nan(n, 1);
    marks.t_peak2_ms       = nan(n, 1);  marks.v_peak2   = nan(n, 1);
    marks.created          = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

    for k = 1:n
        if ~present(k); continue; end
        if isfinite(quota) && sum(marks.done(str == str(k))) >= quota; continue; end
        if so.Skip(k)
            marks.skipped(k) = true;
            continue;
        end
        marks.done(k) = true;
        r = manifest.pixels(k, 1); c = manifest.pixels(k, 2);
        ta  = beat_t + delay_ms(r, c);

        % Marking precision falls as SNR falls — this is what produces the
        % SNR-dependent agreement the harness is designed to expose.
        jit = 0.6 + 6 / max(snr_map(r, c), 1);

        marks.t_act_ms(k)  = ta + jit * randn;
        marks.t_rep_ms(k)  = ta + true_apd(r, c) + jit * randn;
        % A reviewer clicks the peak after their own activation click, so the
        % simulated peak follows it (mv_mark refuses the reverse order).  Only
        % ordering depends on t_peak here; v_peak is fixed below.
        marks.t_peak_ms(k) = max(ta + 2, marks.t_act_ms(k) + 0.5);
        marks.t_base_ms(k) = ta - 8;
        marks.v_peak(k)    = 1;
        marks.v_base(k)    = 0;
        marks.clicks{k}    = [marks.t_base_ms(k) 0; marks.t_peak_ms(k) 1; ...
                              marks.t_act_ms(k) 0.5; marks.t_rep_ms(k) 0.2];
        if isfinite(quota)
            if all(arrayfun(@(q) sum(marks.done(str == q)), 1:manifest.num_strata) >= quota); break; end
        elseif sum(marks.done) >= so.Stop
            break;
        end
    end

    f = fullfile(root, sprintf('sim_%s%s_marks.mat', tag, reviewer));
    save(f, 'marks');
    fprintf('   %s: %d marked, %d skipped\n', reviewer, sum(marks.done), sum(marks.skipped));
end


function check(name, cond, detail)
    if cond
        fprintf('   PASS  %-38s (%s)\n', name, detail);
    else
        error('mv_selftest:failed', 'FAIL  %s (%s)', name, detail);
    end
end
