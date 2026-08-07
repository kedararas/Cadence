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
%   Name-value
%     'Field'      ep_metrics field to compare against (default 'apd_data')
%     'ManualVar'  column of tbl to compare (default: 'apd_ms' for apd mode,
%                  'alt_ratio' for alternans mode)
%     'FileIndex'  index into the ep_metrics cell (default: auto-detect)
%     'CamIndex'   index into the ep_metrics cell (default: auto-detect)
%     'Tolerance'  tolerances to report, in the metric's units (default [5 10])
%     'Plot'       true (default) — Bland-Altman + agreement vs SNR
%
%   Cell orientation: ep_metrics cells are indexed inconsistently across the
%   codebase ({file,cam} in some places, {cam,file} in others), so by default
%   this resolves the map POSITIONALLY — it looks for the cell element whose
%   size matches the recording's [rows cols].  It prints what it chose.  Check
%   that line before you trust the numbers; override with FileIndex/CamIndex.

    p = inputParser;
    p.addParameter('Field',     'apd_data', @(v) ischar(v) || isstring(v));
    p.addParameter('ManualVar', '',         @(v) ischar(v) || isstring(v));
    p.addParameter('FileIndex', [],         @(v) isempty(v) || isscalar(v));
    p.addParameter('CamIndex',  [],         @(v) isempty(v) || isscalar(v));
    p.addParameter('Tolerance', [5 10],     @isnumeric);
    p.addParameter('Plot',      true,       @islogical);
    p.parse(varargin{:});
    o = p.Results;
    o.Field = char(o.Field);

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

    % Recording geometry, used to resolve which cell element is the map.  The
    % manifest is authoritative because it was written from the conditioned file
    % the reviewers actually marked.
    nr = NaN; nc = NaN;
    mf = char(ud_field(ud, 'manifest_file', ''));
    if ~isempty(mf) && isfile(mf)
        MM = load(mf);
        nr = MM.manifest.size(1); nc = MM.manifest.size(2);
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

    [map, how] = resolve_map(d.ep_metrics.(o.Field), nr, nc, o.FileIndex, o.CamIndex);
    fprintf('Using ep_metrics.%s -> %s\n', o.Field, how);

    % ---- pair up ----
    lin  = sub2ind([nr nc], tbl.row, tbl.col);
    soft = double(map(lin));
    man  = double(tbl.(o.ManualVar));

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
    stats.interobserver = interobserver(T, o.ManualVar, o.Tolerance);

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


function [map, how] = resolve_map(C, nr, nc, fi, ci)
%RESOLVE_MAP  Pull the [nr x nc] map out of an ep_metrics entry.

    if isnumeric(C)
        map = C; how = sprintf('numeric %s', mat2str(size(C)));
        return;
    end
    if ~iscell(C)
        error('mv_compare:badType', 'ep_metrics entry is %s; expected numeric or cell.', class(C));
    end

    if ~isempty(fi) && ~isempty(ci)
        map = C{fi, ci};
        how = sprintf('cell{%d,%d} %s (user-specified)', fi, ci, mat2str(size(map)));
        return;
    end

    % Positional resolution: which elements are [nr nc] maps?
    cand = [];
    for a = 1:size(C, 1)
        for b = 1:size(C, 2)
            e = C{a, b};
            if isnumeric(e) && ismatrix(e) && isequal(size(e), [nr nc])
                cand(end+1, :) = [a b]; %#ok<AGROW>
            end
        end
    end

    if isempty(cand)
        error('mv_compare:noMap', ...
              ['No element of the %s cell is a %dx%d map. Pass FileIndex/CamIndex ' ...
               'explicitly, or check that the metrics file matches the conditioned file.'], ...
              mat2str(size(C)), nr, nc);
    end
    if size(cand, 1) > 1
        warning('mv_compare:ambiguous', ...
                ['%d elements of the cell are %dx%d maps; using {%d,%d}. Pass ' ...
                 'FileIndex/CamIndex to disambiguate.'], ...
                size(cand, 1), nr, nc, cand(1,1), cand(1,2));
    end
    map = C{cand(1,1), cand(1,2)};
    how = sprintf('cell{%d,%d} %s (auto-resolved by size)', cand(1,1), cand(1,2), mat2str(size(map)));
end


function io = interobserver(T, var, tol)
%INTEROBSERVER  Pairwise differences between reviewers on the same pixel.

    io = struct('n_pairs', 0, 'sd', NaN, 'bias', NaN, 'within', nan(size(tol)));
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
