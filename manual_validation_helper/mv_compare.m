function stats = mv_compare(tbl, metrics_file, varargin)
%MV_COMPARE  Agreement between CADENCE and manual marking, against inter-observer spread.
%
%   stats = mv_compare(tbl, metrics_file, 'Field', 'apd_data')
%   stats = mv_compare(..., 'ManualVar', 'apd_ms', 'Tolerance', [5 10])
%
%   This is the only file in the harness that opens a *-metrics.mat.  It runs
%   after all marking is finished.
%
%   It reports two things, and the second is what makes the first meaningful:
%
%     Software vs manual  — bias, limits of agreement, and the fraction of
%       pixels within each tolerance, overall and per SNR stratum.
%
%     Inter-observer      — the same spread among the human reviewers.  If
%       software-vs-human disagreement is comparable to human-vs-human
%       disagreement, the software sits inside the noise floor of expert manual
%       measurement.  That claim does not depend on anyone accepting that 5 ms
%       is the right tolerance, which is why it is the one worth leading with.
%
%   Inputs
%     tbl           output of mv_derive (all reviewers pooled)
%     metrics_file  the *-metrics.mat produced by Feature Extraction
%
%   DURATIONS VS TIMES — read this before comparing anything new
%   CADENCE metrics come in two kinds, and only one of them can be differenced
%   against a manual mark directly.
%
%     DURATIONS (apd_data, ca_data, rise times, ca_tau, vc_delay) are
%       differences of two times, so the window origin cancels on both sides.
%       These pair directly.
%
%     WINDOW-RELATIVE TIMES (act_times, rep_data, ca_rep_data) are measured in
%       frames from the start of CADENCE's OWN analysis window, whereas a
%       reviewer's click is an absolute time in the recording.  Differencing
%       them without correction yields a constant bias of (2 - start_frame)
%       frames — hundreds of milliseconds — and Bland-Altman still looks
%       immaculate, because the offset is constant.  This function now applies
%       the correction and REFUSES the comparison when it cannot verify it.
%
%   Name-value
%     'Field'      ep_metrics field to compare against (default 'apd_data')
%     'ManualVar'  column of tbl to compare (default: 'apd_ms' for apd mode,
%                  'alt_ratio' for alternans mode)
%     'Camera'     which camera's map to use.  Default: taken from the
%                  manifest's cam ("CAM2" -> 2), which is the camera the
%                  reviewers actually marked.
%     'Slot'       which of the six data slots holds the map (default 1).
%     'Tolerance'  tolerances to report, in the metric's units (default [5 10])
%     'Plot'       true (default) — Bland-Altman + agreement vs SNR
%
%   CELL LAYOUT.  ep_metrics entries are {camera x slot}: the ROW is the camera
%   (the extract_* drivers loop `for i = 1:num_files` and pass `i` straight in
%   as the camera argument), and the COLUMN is a data slot 1-6 — slot 1 is the
%   masked map, 2 the background image, 5 the camera label, 6 the unmasked map.
%   The old FileIndex/CamIndex names described this backwards and are rejected:
%   passing CamIndex=2 returned the BACKGROUND IMAGE, which is numeric and
%   correctly sized, so it produced confident nonsense.

    p = inputParser;
    p.addParameter('Field',     'apd_data', @(v) ischar(v) || isstring(v));
    p.addParameter('ManualVar', '',         @(v) ischar(v) || isstring(v));
    p.addParameter('Camera',    [],         @(v) isempty(v) || isscalar(v));
    p.addParameter('Slot',      1,          @(v) isnumeric(v) && isscalar(v));
    p.addParameter('Tolerance', [5 10],     @isnumeric);
    p.addParameter('Plot',      true,       @islogical);
    p.addParameter('FileIndex', [],         @(v) isempty(v));
    p.addParameter('CamIndex',  [],         @(v) isempty(v));
    p.parse(varargin{:});
    o = p.Results;
    o.Field = char(o.Field);

    if ~isempty(o.FileIndex) || ~isempty(o.CamIndex)
        error('mv_compare:renamedIndex', ...
              ['FileIndex/CamIndex have been replaced by Camera/Slot, because the old ' ...
               'names had the layout backwards. ep_metrics cells are {camera, slot}: ' ...
               'pass ''Camera'', <camera number> and ''Slot'', 1 for the masked map.']);
    end

    ud = tbl.Properties.UserData;
    if isempty(o.ManualVar)
        if isfield(ud, 'mode') && strcmp(ud.mode, 'alternans')
            o.ManualVar = 'alt_ratio';
        else
            o.ManualVar = 'apd_ms';
        end
    end
    o.ManualVar = char(o.ManualVar);
    if ~ismember(o.ManualVar, tbl.Properties.VariableNames)
        error('mv_compare:noVar', '''%s'' is not a column of the derived table.', o.ManualVar);
    end

    % ---- CADENCE's map ----
    S = load(metrics_file);
    if isfield(S, 'cmos_all_data'); d = S.cmos_all_data; else
        fn = fieldnames(S); is_struct = structfun(@isstruct, S);
        d = S.(fn{find(is_struct, 1)});
    end
    if ~isfield(d, 'ep_metrics') || ~isfield(d.ep_metrics, o.Field)
        error('mv_compare:noField', 'ep_metrics.%s not present in %s.', o.Field, metrics_file);
    end

    % Recording geometry and camera.  The manifest is authoritative because it
    % was written from the conditioned file the reviewers actually marked.
    manifest = [];
    nr = NaN; nc = NaN;
    mf = char(ud_field(ud, 'manifest_file', ''));
    if ~isempty(mf) && isfile(mf)
        MM = load(mf);
        manifest = MM.manifest;
        nr = manifest.size(1); nc = manifest.size(2);
    end
    if ~isfinite(nr)
        cam = char(ud_field(ud, 'cam', 'CAM1'));
        if isfield(d, cam)
            [nr, nc, ~] = size(d.(cam));
        else
            nr = max(tbl.row); nc = max(tbl.col);
            warning('mv_compare:geometry', ...
                    ['Could not read geometry from the manifest or the metrics file; ' ...
                     'inferring %dx%d from the marked pixels. Verify the resolved map below.'], nr, nc);
        end
    end

    % Which camera?  Auto-resolving by size cannot tell CAM1's map from CAM2's
    % when both are populated and the same shape, and it silently takes the
    % first — so on a rig with two voltage cameras it would compare the wrong
    % one.  The manifest records the camera the reviewers marked; use it.
    camera = o.Camera;
    if isempty(camera)
        camstr = char(ud_field(ud, 'cam', ''));
        if isempty(camstr) && ~isempty(manifest); camstr = char(manifest.cam); end
        camera = str2double(regexprep(camstr, '\D', ''));
        if ~isfinite(camera) || camera < 1
            error('mv_compare:noCamera', ...
                  ['Could not determine the camera from the manifest (cam = ''%s''). ' ...
                   'Pass ''Camera'', <number> explicitly.'], camstr);
        end
    end

    [map, how] = resolve_map(d.ep_metrics.(o.Field), nr, nc, camera, o.Slot);
    fprintf('Using ep_metrics.%s -> %s\n', o.Field, how);

    % ---- pair up ----
    lin  = sub2ind([nr nc], tbl.row, tbl.col);
    soft = double(map(lin));
    man  = double(tbl.(o.ManualVar));

    % Kind check and origin correction.  See the header: a window-relative time
    % cannot be differenced against an absolute manual click without shifting
    % it into the recording's timeline first, and the failure is silent.
    [soft, kind_note] = align_origin(soft, o.Field, o.ManualVar, manifest, d);
    if ~isempty(kind_note); fprintf('%s\n', kind_note); end

    % act_times is masked with ZEROS rather than NaN (act .* mask, where the
    % mask is 1/NaN elsewhere but 0 here), so off-tissue pixels arrive as a
    % finite 0 and would pass the isfinite filter and drag the bias negative.
    % A genuine LAT is at least one frame, so 0 means "no value".
    if any(strcmp(o.Field, {'act_times'}))
        nz = (soft == 0);
        if any(nz)
            fprintf('  %d pixel(s) with act_times == 0 treated as missing (zero-masked background).\n', sum(nz));
            soft(nz) = NaN;
        end
    end

    ok = isfinite(soft) & isfinite(man);
    if ~any(ok)
        error('mv_compare:noPairs', ...
              ['No pixel has both a manual and a CADENCE value. The map may be all-NaN ' ...
               'at the sampled pixels, or the cell element resolved above is the wrong one.']);
    end
    if sum(~ok) > 0
        fprintf('  %d of %d rows dropped (no CADENCE value at that pixel).\n', sum(~ok), numel(ok));
    end

    T = tbl(ok, :); soft = soft(ok); man = man(ok);
    diffs = soft - man;

    % ---- software vs manual ----
    stats = struct();
    stats.field       = o.Field;
    stats.manual_var  = o.ManualVar;
    stats.n_pairs     = numel(diffs);
    stats.n_pixels    = numel(unique(T.pixel));
    stats.n_reviewers = numel(unique(T.reviewer));
    stats.bias        = mean(diffs);
    stats.sd          = std(diffs);
    stats.loa         = stats.bias + [-1.96 1.96] * stats.sd;
    stats.tolerance   = o.Tolerance;
    stats.within      = arrayfun(@(t) 100 * mean(abs(diffs) <= t), o.Tolerance);

    r = corrcoef(soft, man);
    stats.r = r(1, 2);

    % ---- per-stratum (the operating envelope) ----
    strata = unique(T.stratum);
    stats.by_stratum = struct('stratum', {}, 'snr_median', {}, 'n', {}, ...
                              'bias', {}, 'sd', {}, 'within', {});
    for k = 1:numel(strata)
        sel = T.stratum == strata(k);
        stats.by_stratum(k).stratum    = strata(k);
        stats.by_stratum(k).snr_median = median(T.snr(sel));
        stats.by_stratum(k).n          = sum(sel);
        stats.by_stratum(k).bias       = mean(diffs(sel));
        stats.by_stratum(k).sd         = std(diffs(sel));
        stats.by_stratum(k).within     = arrayfun(@(t) 100 * mean(abs(diffs(sel)) <= t), o.Tolerance);
    end

    % ---- inter-observer ----
    % Computed on ALL marked rows, not the software-resolved subset.  This is
    % the noise floor of manual measurement, which exists whether or not
    % CADENCE returned a value; conditioning it on the pixels CADENCE handled
    % removes exactly the hard ones from the human side and flatters the
    % comparison the whole report leads with.
    stats.interobserver = interobserver(tbl, o.ManualVar, o.Tolerance);

    % ---- report ----
    fprintf('\n--- %s vs manual %s ---\n', o.Field, o.ManualVar);
    fprintf('  pairs        : %d  (%d pixels x %d reviewers)\n', ...
            stats.n_pairs, stats.n_pixels, stats.n_reviewers);
    fprintf('  bias (SW-manual): %+.2f\n', stats.bias);
    fprintf('  95%% LoA      : [%+.2f, %+.2f]\n', stats.loa(1), stats.loa(2));
    fprintf('  r            : %.3f\n', stats.r);
    for k = 1:numel(o.Tolerance)
        fprintf('  within %-5g : %.1f%%\n', o.Tolerance(k), stats.within(k));
    end

    if ~isempty(stats.interobserver.n_pairs) && stats.interobserver.n_pairs > 0
        fprintf('\n--- inter-observer (human vs human) ---\n');
        fprintf('  pairs        : %d\n', stats.interobserver.n_pairs);
        fprintf('  SD of diffs  : %.2f   (software vs manual: %.2f)\n', ...
                stats.interobserver.sd, stats.sd);
        for k = 1:numel(o.Tolerance)
            fprintf('  within %-5g : %.1f%%   (software: %.1f%%)\n', ...
                    o.Tolerance(k), stats.interobserver.within(k), stats.within(k));
        end
        if stats.sd <= stats.interobserver.sd
            fprintf('  => software disagreement is within the inter-observer spread.\n');
        else
            fprintf('  => software disagreement EXCEEDS the inter-observer spread by %.2fx.\n', ...
                    stats.sd / stats.interobserver.sd);
        end
    else
        fprintf('\n  (only one reviewer — no inter-observer comparison; use >=2 for the\n');
        fprintf('   benchmark that makes the agreement numbers interpretable)\n');
    end

    fprintf('\n--- by SNR stratum ---\n');
    for k = 1:numel(stats.by_stratum)
        b = stats.by_stratum(k);
        fprintf('  stratum %d (median SNR %5.1f, n=%4d): bias %+6.2f, SD %5.2f, within %g: %5.1f%%\n', ...
                b.stratum, b.snr_median, b.n, b.bias, b.sd, o.Tolerance(1), b.within(1));
    end

    if o.Plot
        plot_agreement(soft, man, diffs, T, stats, o);
    end
end


% ======================================================================
function v = ud_field(ud, f, dflt)
    if isstruct(ud) && isfield(ud, f) && ~isempty(ud.(f)); v = ud.(f); else; v = dflt; end
end


function [map, how] = resolve_map(C, nr, nc, camera, slot)
%RESOLVE_MAP  Pull the [nr x nc] map out of an ep_metrics {camera, slot} entry.
%
%   Addressed, not guessed.  The previous version scanned for any element whose
%   size matched [nr nc] and took the first — which also matches the background
%   image in slot 2, and matches camera 1's map when the caller wanted camera 2.
%   Both failures are silent and produce plausible statistics.

    if isnumeric(C)
        if ~isequal(size(C), [nr nc])
            error('mv_compare:sizeMismatch', ...
                  'ep_metrics entry is numeric %s but the recording is %dx%d.', ...
                  mat2str(size(C)), nr, nc);
        end
        map = C; how = sprintf('numeric %s', mat2str(size(C)));
        return;
    end
    if ~iscell(C)
        error('mv_compare:badType', 'ep_metrics entry is %s; expected numeric or cell.', class(C));
    end

    if camera > size(C, 1) || slot > size(C, 2)
        error('mv_compare:outOfRange', ...
              'Requested {camera %d, slot %d} but the entry is %s.', ...
              camera, slot, mat2str(size(C)));
    end

    map = C{camera, slot};
    if isempty(map)
        error('mv_compare:emptyCell', ...
              ['ep_metrics{%d,%d} is empty — that metric was not extracted for camera %d. ' ...
               'Check the camera number, or that the box was ticked when the metrics were made.'], ...
              camera, slot, camera);
    end
    if ~isnumeric(map) || ~ismatrix(map) || ~isequal(size(map), [nr nc])
        error('mv_compare:notAMap', ...
              ['ep_metrics{%d,%d} is %s %s, not a %dx%d map. Slot 1 is the masked map ' ...
               'and slot 6 the unmasked one; slot 2 is the background image.'], ...
              camera, slot, class(map), mat2str(size(map)), nr, nc);
    end
    how = sprintf('cell{camera %d, slot %d} %s', camera, slot, mat2str(size(map)));
end


function kind = metric_kind(field)
%METRIC_KIND  'duration' | 'window_time' — see the header of mv_compare.
%
%   A duration is a difference of two times, so any window origin cancels.  A
%   window_time is measured in frames from the start of CADENCE's own analysis
%   window and must be shifted into the recording's timeline before it can be
%   compared with a reviewer's click.

    switch field
        case {'apd_data', 'ca_data', 'ca_apd_data', 'ap_rise_times', ...
              'ca_rise_times', 'ca_tau', 'vc_delay', 'alternans_data', 'local_cv'}
            kind = 'duration';
        case {'act_times', 'rep_data', 'ca_rep_data'}
            kind = 'window_time';
        otherwise
            kind = 'unknown';
    end
end


function [soft, note] = align_origin(soft, field, manual_var, manifest, d)
%ALIGN_ORIGIN  Shift a window-relative CADENCE time into absolute recording ms.
%
%   The manual value is an absolute time: mv_mark builds its axis as
%   ((wf(1):wf(2)) - 1)/acqFreq*1000, so a click at absolute frame F is
%   (F-1)*dt.  CADENCE stores compute_lat_50's 1-BASED fractional frame index
%   inside its own window, times dt: (F - sf + 1)*dt for window start sf.  The
%   uncorrected difference is therefore a constant (2 - sf) frames, which for a
%   window starting at frame 1001 at 1 kHz is -999 ms on every pixel, with
%   perfect-looking scatter.  Correct by adding (sf - 2)*dt.

    note = '';
    kind = metric_kind(field);

    if strcmp(kind, 'duration')
        if any(strcmp(manual_var, {'act_ms'}))
            error('mv_compare:kindMismatch', ...
                  ['''%s'' is a duration but ManualVar ''%s'' is an absolute time. ' ...
                   'Comparing them is meaningless.'], field, manual_var);
        end
        return;
    end

    if strcmp(kind, 'unknown')
        warning('mv_compare:unknownKind', ...
                ['''%s'' is not in the duration/time registry, so no origin check was ' ...
                 'made. If it stores a TIME rather than a DURATION, the bias below is ' ...
                 'the window offset, not a measurement error.'], field);
        return;
    end

    % --- window_time: correction required ---
    if isempty(manifest)
        error('mv_compare:noManifest', ...
              ['''%s'' is a window-relative time and needs the marking window to be ' ...
               'shifted into absolute time, but the manifest could not be loaded.'], field);
    end

    src = 'raw';
    if isfield(manifest, 'source') && ~isempty(manifest.source)
        src = char(manifest.source);
    end
    fs = double(manifest.acqFreq);
    dt = 1000 / fs;

    if strcmp(src, 'ensemble')
        % Both sides are on the SAME averaged beat.  Feature extraction takes
        % CAM<n>_average whenever the ensemble box is ticked, and the reviewers
        % marked that same array, so the origin is frame 1 on both sides and
        % nothing about the recording timeline enters.  The only residual is
        % the index convention: CADENCE stores compute_lat_50's 1-BASED
        % fractional frame index times dt, while mv_mark's axis is (frame-1)*dt.
        sf = 1;
        note = sprintf(['  origin: ensemble-averaged beat, shared frame-1 origin; ' ...
                        'correcting the 1-based convention by %+.1f ms'], (sf - 2) * dt);
    else
        if ~isfield(d, 'window') || isempty(d.window) || ~iscell(d.window) || isempty(d.window{1,1})
            error('mv_compare:noWindow', ...
                  ['''%s'' is a window-relative time and the marks were made on a RAW single ' ...
                   'beat, but the metrics file records no analysis window. Either re-extract ' ...
                   'with a drawn window, or — better — mark the ensemble average instead ' ...
                   '(mv_sample_pixels default), which shares its origin with feature ' ...
                   'extraction and needs no correction.'], field);
        end
        sf = double(d.window{1,1}(1));
        wf = manifest.window_frames;
        if sf > wf(2) || double(d.window{1,1}(2)) < wf(1)
            error('mv_compare:differentBeat', ...
                  ['CADENCE analysed frames %d-%d but the reviewers marked frames %d-%d — ' ...
                   'different beats, so a shared origin does not exist. Re-run ' ...
                   'mv_sample_pixels with ''MetricsFile'' pointing at this metrics file.'], ...
                  sf, double(d.window{1,1}(2)), wf(1), wf(2));
        end
        note = sprintf(['  origin: shifted CADENCE window-relative time by %+.1f ms ' ...
                        '(window starts at frame %d)'], (sf - 2) * dt, sf);
    end

    soft = soft + (sf - 2) * dt;

    % Sanity check the result rather than trusting the declared source.  The
    % metrics file does NOT record whether feature extraction ran on the
    % average or on a windowed raw beat (`ensemble || isempty(oap_window)`
    % chooses at run time and nothing is saved), so this is the only available
    % guard against the two sides having looked at different arrays.
    span = manifest.window_ms;
    frac_out = mean(soft < span(1) - 5*dt | soft > span(2) + 5*dt, 'omitnan');
    if frac_out > 0.5
        warning('mv_compare:originSuspect', ...
                ['%.0f%% of corrected values fall outside the marked window (%.1f-%.1f ms). ' ...
                 'Feature extraction was probably run on a different array than the reviewers ' ...
                 'marked — check the ensemble checkbox state that produced this metrics file.'], ...
                100 * frac_out, span(1), span(2));
    end
end


function io = interobserver(T, var, tol)
%INTEROBSERVER  Pairwise differences between reviewers on the same pixel.

    io = struct('n_pairs', 0, 'sd', NaN, 'bias', NaN, 'within', nan(size(tol)));
    T = T(isfinite(T.(var)), :);
    revs = unique(T.reviewer);
    if numel(revs) < 2; return; end

    dd = [];
    for a = 1:numel(revs)
        for b = a+1:numel(revs)
            Ta = T(T.reviewer == revs(a), :);
            Tb = T(T.reviewer == revs(b), :);
            [~, ia, ib] = intersect(Ta.pixel, Tb.pixel);
            if isempty(ia); continue; end
            dd = [dd; Ta.(var)(ia) - Tb.(var)(ib)]; %#ok<AGROW>
        end
    end
    if isempty(dd); return; end

    io.n_pairs = numel(dd);
    io.bias    = mean(dd);
    io.sd      = std(dd);
    io.within  = arrayfun(@(t) 100 * mean(abs(dd) <= t), tol);
end


function plot_agreement(soft, man, diffs, T, stats, o)
    figure('Name', 'mv_compare — agreement', 'Color', 'w', 'Position', [100 100 1100 420]);

    % Bland-Altman
    subplot(1, 3, 1);
    mn = (soft + man) / 2;
    scatter(mn, diffs, 14, T.snr, 'filled'); hold on;
    yline(stats.bias,  '-',  sprintf('bias %+.2f', stats.bias), 'LineWidth', 1.3);
    yline(stats.loa(1), '--', sprintf('%+.2f', stats.loa(1)));
    yline(stats.loa(2), '--', sprintf('%+.2f', stats.loa(2)));
    xlabel('Mean of software and manual'); ylabel('Software - manual');
    title('Bland-Altman'); grid on;
    cb = colorbar; cb.Label.String = 'SNR';

    % Identity
    subplot(1, 3, 2);
    scatter(man, soft, 14, T.snr, 'filled'); hold on;
    lims = [min([man; soft]) max([man; soft])];
    plot(lims, lims, 'k--', 'LineWidth', 1.1);
    xlabel(sprintf('Manual %s', strrep(o.ManualVar, '_', '\_')));
    ylabel(sprintf('CADENCE %s', strrep(o.Field, '_', '\_')));
    title(sprintf('r = %.3f', stats.r)); grid on; axis square;
    xlim(lims); ylim(lims);

    % Agreement vs SNR stratum — the operating envelope
    subplot(1, 3, 3);
    b = stats.by_stratum;
    x = arrayfun(@(s) s.snr_median, b);
    y = arrayfun(@(s) s.within(1),  b);
    e = arrayfun(@(s) s.sd,         b);
    yyaxis left;  plot(x, y, 'o-', 'LineWidth', 1.5); ylabel(sprintf('%% within %g', o.Tolerance(1)));
    ylim([0 100]);
    yyaxis right; plot(x, e, 's--', 'LineWidth', 1.2); ylabel('SD of difference');
    xlabel('Median SNR of stratum'); title('Agreement vs signal quality'); grid on;
end
