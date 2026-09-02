function results = mv_cv_envelope(varargin)
%MV_CV_ENVELOPE  CV accuracy envelope across the real preparation geometries.
%
%   results = mv_cv_envelope
%   results = mv_cv_envelope('OutDir', 'figures')
%   results = mv_cv_envelope('Results', 'mv_cv_envelope.mat')   % re-plot only
%
%   MV_SYNTH_CV establishes that the CV engines are exact on analytic maps, so
%   in practice CV accuracy is set entirely by how precisely activation times
%   can be measured.  This script turns that into the number a user actually
%   needs: for MY preparation, what CV error should I expect?
%
%   It runs mv_synth_cv at the pixel size of each real CADENCE geometry and
%   writes the supplementary tables and the jitter figures.  Reporting one
%   envelope for all preparations would be wrong by up to 5x, because CV error
%   scales as sigma_LAT * v / ds -- the jitter measured against the activation
%   step between ADJACENT pixels (ds/v), not against the frame period.
%
%   Consequence worth stating explicitly in any write-up: FINER pixels make CV
%   MORE sensitive to activation-time error, because the per-pixel LAT step
%   shrinks.  The 10 mm slice is the hardest case, not the easiest.
%
%   OUTPUTS WRITTEN (into OutDir)
%     cv_jitter_<n>preps.csv   test C, long format, one row per engine x sigma
%                              x preparation; carries lat_step_ms and
%                              jitter_per_step for the collapsed plot
%     cv_aniso_<n>preps.csv    test B, L/T/anisotropy-ratio recovery
%     cv_planar_<n>preps.csv   test A, planar exactness
%     fig_cv_jitter_panels.png small multiples, one panel per preparation
%     fig_cv_jitter_collapse.png  all series vs sigma_LAT/(ds/v)
%     mv_cv_envelope.mat       raw mv_synth_cv results per preparation
%
%   Tests A and B are geometry-independent (the analytic maps scale with the
%   grid), so those two tables are identical across preparations and should be
%   reported once rather than four times.  Only the jitter sweep (C) moves.
%
%   PLOTTING CONVENTION
%     Points where fewer than 95% of pixels survive are drawn as OPEN markers.
%     Past its rejection cliff the circle engine discards a large fraction of
%     pixels, so its median describes only the survivors and can even fall as
%     jitter rises.  That is survivor bias, not accuracy, and it must not be
%     plotted as though it were.
%
%   Name-value
%     'FOV'        field of view per preparation, mm (default [10 20 30 50])
%     'Labels'     display names (default slice / rat WH / rabbit WH / wedge)
%     'GridSize'   pixels per side; scalar, or one value per preparation
%                  (default 256). The human LV wedge is acquired at 100x100
%                  while the other preparations are 256x256, so pixel pitch is
%                  NOT simply FOV/256 for every row -- pass a vector.
%     'JitterV'    planar speed for the jitter sweep, cm/s (default 60)
%     'OutDir'     where to write tables and figures (default pwd)
%     'Results'    path to a previous .mat; skips the sweep and re-plots only
%     'Save'       write files (default true); false returns results in memory
%
%   Output
%     results  struct with .jitter, .aniso, .planar tables (all preparations
%              concatenated) and .raw (per-preparation mv_synth_cv output)

    p = inputParser;
    p.addParameter('FOV',      [10 20 30 50], @isnumeric);
    p.addParameter('Labels',   {}, @iscell);
    p.addParameter('GridSize', 256, @isnumeric);
    p.addParameter('JitterV',  60,  @isscalar);
    p.addParameter('OutDir',   pwd, @(s) ischar(s) || isstring(s));
    p.addParameter('Results',  '',  @(s) ischar(s) || isstring(s));
    p.addParameter('Save',     true, @islogical);
    p.parse(varargin{:});
    o = p.Results;

    ensure_path();

    labels = o.Labels;
    if isempty(labels)
        defaults = {'Slice (10 mm)','Rat WH (20 mm)','Rabbit WH (30 mm)','Wedge (50 mm)'};
        if numel(o.FOV) == 4
            labels = defaults;
        else
            labels = arrayfun(@(f) sprintf('FOV %g mm', f), o.FOV, 'UniformOutput', false);
        end
    end
    assert(numel(labels) == numel(o.FOV), 'Labels must match FOV in length.');

    outdir = char(o.OutDir);
    if o.Save && ~isempty(outdir) && ~exist(outdir, 'dir'), mkdir(outdir); end

    % ---------- sweep (or reload) ----------
    if ~isempty(o.Results)
        S = load(char(o.Results));
        raw = S.raw;  v_jit = S.v_jit;
        assert(numel(raw) == numel(o.FOV), 'Saved results do not match FOV list.');
    else
        v_jit = o.JitterV;
        raw = cell(1, numel(o.FOV));
        for i = 1:numel(o.FOV)
            g  = grid_for(o.GridSize, i);
            ds = o.FOV(i) / g;
            fprintf('\n======== %s : FOV %g mm, %d x %d px, %.4f mm/px ========\n', ...
                    labels{i}, o.FOV(i), g, g, ds);
            raw{i} = mv_synth_cv('GridSize', [g g], ...
                                 'PixelSize', ds, 'JitterV', v_jit, 'EndToEnd', false);
        end
    end

    % ---------- assemble long-format tables ----------
    [J, A, P] = deal(table());
    for i = 1:numel(o.FOV)
        g    = grid_for(o.GridSize, i);
        ds   = o.FOV(i) / g;
        step = ds / (v_jit / 100);                 % adjacent-pixel LAT step, ms

        t = raw{i}.jitter;
        t.prep            = repmat(string(labels{i}), height(t), 1);
        t.fov_mm          = repmat(o.FOV(i), height(t), 1);
        t.px_mm           = repmat(ds, height(t), 1);
        t.grid_px         = repmat(g, height(t), 1);
        t.lat_step_ms     = repmat(step, height(t), 1);
        t.jitter_per_step = t.lat_jitter_ms / step;
        J = [J; t];                                                     %#ok<AGROW>

        t = raw{i}.aniso;
        t.prep  = repmat(string(labels{i}), height(t), 1);
        t.px_mm = repmat(ds, height(t), 1);
        A = [A; t];                                                     %#ok<AGROW>

        t = raw{i}.planar;
        t.prep  = repmat(string(labels{i}), height(t), 1);
        t.px_mm = repmat(ds, height(t), 1);
        P = [P; t];                                                     %#ok<AGROW>
    end

    results = struct('jitter', J, 'aniso', A, 'planar', P, 'raw', {raw}, ...
                     'labels', {labels}, 'fov', o.FOV, 'jitter_v_cms', v_jit);

    if ~o.Save, return; end

    n = numel(o.FOV);
    writetable(J, fullfile(outdir, sprintf('cv_jitter_%dpreps.csv', n)));
    writetable(A, fullfile(outdir, sprintf('cv_aniso_%dpreps.csv',  n)));
    writetable(P, fullfile(outdir, sprintf('cv_planar_%dpreps.csv', n)));
    save(fullfile(outdir, 'mv_cv_envelope.mat'), 'raw', 'v_jit', 'labels');

    plot_panels(J, labels, outdir);
    plot_collapse(J, labels, outdir);
    fprintf('\nWrote tables and figures to %s\n', outdir);
end


% ======================================================================
function plot_panels(J, labels, outdir)
%PLOT_PANELS  Small multiples: one panel per preparation, one line per engine.
%   Panels share a y-axis on purpose -- the whole point is that sensitivity
%   grows as pixels get finer, and that comparison is unreadable if each panel
%   is autoscaled.
    engines = unique(J.engine, 'stable');
    ecol = engine_colors(engines);
    np = numel(labels);

    f = figure('Position', [100 100 300*np 330], 'Color', 'w');
    for i = 1:np
        subplot(1, np, i); hold on; box on
        t = J(J.prep == string(labels{i}), :);
        for e = 1:numel(engines)
            s = sortrows(t(t.engine == engines(e), :), 'lat_jitter_ms');
            plot(s.lat_jitter_ms, s.med_abs_err_pct, '-', ...
                 'Color', ecol(e,:), 'LineWidth', 1.6);
            ok = s.pixels_kept_pct >= 95;
            plot(s.lat_jitter_ms(ok), s.med_abs_err_pct(ok), 'o', ...
                 'MarkerFaceColor', ecol(e,:), 'MarkerEdgeColor', ecol(e,:), 'MarkerSize', 5);
            plot(s.lat_jitter_ms(~ok), s.med_abs_err_pct(~ok), 'o', ...
                 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ecol(e,:), ...
                 'MarkerSize', 5, 'LineWidth', 1.1);
        end
        yline(5, 'k:', 'LineWidth', 1);
        xlim([0 max(J.lat_jitter_ms)]); ylim([0 100]);
        xlabel('\sigma_{LAT} (ms)');
        if i == 1, ylabel('median |CV error| (%)'); else, set(gca, 'YTickLabel', []); end
        title(sprintf('%s\n%.3f mm/px', labels{i}, t.px_mm(1)), ...
              'FontWeight', 'normal', 'FontSize', 9);
        if i == np
            h = gobjects(numel(engines), 1);
            for e = 1:numel(engines)
                h(e) = plot(NaN, NaN, '-o', 'Color', ecol(e,:), ...
                            'MarkerFaceColor', ecol(e,:), 'LineWidth', 1.6);
            end
            legend(h, cellstr(engines), 'Location', 'southeast', 'Box', 'off', 'FontSize', 8);
        end
        set(gca, 'FontSize', 9);
    end
    exportgraphics(f, fullfile(outdir, 'fig_cv_jitter_panels.png'), 'Resolution', 300);
end


function plot_collapse(J, labels, outdir)
%PLOT_COLLAPSE  All series against sigma_LAT/(ds/v), which removes geometry.
%   If the scaling law holds, the twelve series become one line per engine and
%   a reader with any other FOV can read their own envelope off it.  Deviations
%   are informative rather than embarrassing: the medians roll off below the
%   line once noise dominates, and the circle engine's rejected points sit
%   above it.
    engines = unique(J.engine, 'stable');
    ecol = engine_colors(engines);
    mk = {'o','s','^','d','v','p'};

    f = figure('Position', [100 100 470 390], 'Color', 'w'); hold on; box on
    for e = 1:numel(engines)
        for i = 1:numel(labels)
            s = sortrows(J(J.prep == string(labels{i}) & J.engine == engines(e), :), ...
                         'jitter_per_step');
            s = s(s.jitter_per_step > 0, :);
            ok = s.pixels_kept_pct >= 95;
            m = mk{min(i, numel(mk))};
            plot(s.jitter_per_step(ok), s.med_abs_err_pct(ok), m, ...
                 'Color', ecol(e,:), 'MarkerFaceColor', ecol(e,:), 'MarkerSize', 5);
            plot(s.jitter_per_step(~ok), s.med_abs_err_pct(~ok), m, ...
                 'Color', ecol(e,:), 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1);
        end
    end
    xr = [max(min(J.jitter_per_step(J.jitter_per_step > 0)) * 0.7, eps), ...
          max(J.jitter_per_step) * 1.4];
    xx = logspace(log10(xr(1)), log10(xr(2)), 100);
    plot(xx, 2.05 * xx, 'k--', 'LineWidth', 1);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim(xr); ylim([0.8 150]);
    xlabel('\sigma_{LAT} / (ds/v)   [jitter per adjacent-pixel LAT step]');
    ylabel('median |CV error| (%)');
    title('All preparations collapse onto one law', 'FontWeight', 'normal', 'FontSize', 10);

    h  = gobjects(numel(engines), 1);
    for e = 1:numel(engines)
        h(e) = plot(NaN, NaN, 'o', 'Color', ecol(e,:), 'MarkerFaceColor', ecol(e,:));
    end
    h2 = gobjects(numel(labels), 1);
    for i = 1:numel(labels)
        h2(i) = plot(NaN, NaN, mk{min(i, numel(mk))}, 'Color', 'k');
    end
    legend([h; h2], [cellstr(engines); labels(:)], 'Location', 'northwest', ...
           'Box', 'off', 'FontSize', 7, 'NumColumns', 2);
    set(gca, 'FontSize', 9);
    exportgraphics(f, fullfile(outdir, 'fig_cv_jitter_collapse.png'), 'Resolution', 300);
end


function c = engine_colors(engines)
%ENGINE_COLORS  Stable colour per engine so figures agree across the paper.
    map = struct('gradient', [0.85 0.33 0.10], ...
                 'polyfit',  [0.00 0.45 0.74], ...
                 'circle',   [0.47 0.67 0.19]);
    fallback = lines(numel(engines));
    c = zeros(numel(engines), 3);
    for i = 1:numel(engines)
        name = char(engines(i));
        if isfield(map, name), c(i,:) = map.(name); else, c(i,:) = fallback(i,:); end
    end
end


function g = grid_for(gs, i)
%GRID_FOR  Grid size for preparation i: scalar applies to all, else per-row.
    if isscalar(gs), g = gs; else, g = gs(min(i, numel(gs))); end
end


function ensure_path()
%ENSURE_PATH  Add sibling helper folders if the CV engines are not visible.
%   The other mv_* scripts assume the caller has set the path; this keeps the
%   figure scripts runnable on their own without changing that convention.
    if exist('conduction_velocity', 'file') == 2, return; end
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    for f = {'conduction_velocity_helper', 'feature_extraction_helper', ...
             'signal_conditioning_helper'}
        d = fullfile(root, f{1});
        if exist(d, 'dir'), addpath(d); end
    end
end
