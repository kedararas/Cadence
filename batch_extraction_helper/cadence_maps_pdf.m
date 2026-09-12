function [n_pages, pdf_file] = cadence_maps_pdf(d, pdf_file, varargin)
%CADENCE_MAPS_PDF  Write the cardiac maps of one recording to a multi-page PDF.
%
%   n_pages = cadence_maps_pdf(metrics, pdf_file)
%   n_pages = cadence_maps_pdf(metrics, pdf_file, Name, Value, ...)
%
%   Headless counterpart of the map windows the Feature Extraction app pops
%   up after EXTRACT FEATURES and saves with save_figures_to_pdf: one page per
%   feature and camera, in the app's order and with the app's panels.
%     feature maps   (activation, APD, repolarization, AP rise, V-Ca delay,
%                     CaTD, Ca decay, Ca rise, Ca tau, plus the batch-only
%                     APD50 / CaTD50 maps)
%                    map over the camera image with the representative pixel
%                    marked, robust colour scale (median +/- 3 MAD-sigma, jet),
%                    histogram, representative trace with the measured
%                    location / curve, and the mean / SD / median / sample
%     CV             map with propagation arrows (~1 mm spacing), speed histogram
%     complexity     DF / RI / OI maps and histograms, trace with pacing
%     alternans      the six substrate panels of the voltage or calcium
%                    branch, trace with pacing, per-panel statistics
%   `metrics` is the struct written to *-metrics.mat (cmos_all_data with
%   ep_metrics), as produced by cadence_extract_features or the app.  Pages
%   are rendered with hidden classic figures and exportgraphics, so this
%   works in a -nodisplay session; they are raster, like the app's PDF.
%
%   Name/Value
%     'Resolution'  dpi of the pages (default 150)
%     'Pixel'       [row col] representative pixel (default image centre, the
%                   batch extractor's choice)
%     'Name'        recording name printed on every page (default from pdf_file)
%     'Log'         function handle for messages (default @(s) fprintf('%s\n', s))
%     'PNGDir'      folder: also write every page as <name>_<nn>_<feature>_CAM<n>.png
%                   there (default '' = PDF only)
%
%   Returns the number of pages written (0 = nothing to plot; no file is left
%   behind).  A page that fails to render is logged and skipped.
%
%   See also cadence_batch_extract ('SaveMapsPDF'), cadence_recompute_medians,
%   save_figures_to_pdf.

    p = inputParser;
    p.addParameter('Resolution', 150, @isnumeric);
    p.addParameter('Pixel', [], @isnumeric);
    p.addParameter('Name', '', @(x) ischar(x) || isstring(x));
    p.addParameter('Log', @(s) fprintf('%s\n', s), @(x) isa(x, 'function_handle'));
    p.addParameter('PNGDir', '', @(x) ischar(x) || isstring(x));
    p.parse(varargin{:});
    o = p.Results;
    logf = o.Log;
    o.PNGDir = char(o.PNGDir);
    if ~isempty(o.PNGDir) && ~isfolder(o.PNGDir), mkdir(o.PNGDir); end

    pdf_file = char(pdf_file);
    [folder, base, ext] = fileparts(pdf_file);
    if ~strcmpi(ext, '.pdf'), pdf_file = fullfile(folder, [base '.pdf']); end
    if isempty(o.Name), o.Name = regexprep(base, '-maps$', ''); end
    name = char(o.Name);
    if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
    if isfile(pdf_file), delete(pdf_file); end     % exportgraphics appends; drop stale pages

    n_pages = 0;
    if ~isfield(d, 'ep_metrics') || isempty(d.ep_metrics)
        logf('maps PDF: no ep_metrics in this recording; nothing to plot.');
        return;
    end
    em = d.ep_metrics;
    if isempty(o.Pixel)
        [r, c] = size(d.CAM1(:,:,1));
        o.Pixel = [round(r/2), round(c/2)];
    end
    pacing = [];
    if isfield(d, 'analog1') && ~isempty(d.analog1), pacing = d.analog1; end
    fov = NaN;
    if isfield(d, 'FOV') && ~isempty(d.FOV), fov = double(d.FOV); end

    % the app's order (EXTRACT FEATURES callback), then the batch-only maps
    pages = { ...
        'act_times',       'Act Times',              'feature'; ...
        'apd_data',        'APD Data',               'feature'; ...
        'rep_data',        'Repolarization Data',    'feature'; ...
        'ap_rise_times',   'AP Rise Time',           'feature'; ...
        'local_cv',        'CV',                     'cv'; ...
        'vc_delay',        'Voltage-Ca Delay',       'feature'; ...
        'ca_data',         'Ca Transient Duration',  'feature'; ...
        'ca_rep_data',     'Ca Decay Time',          'feature'; ...
        'ca_rise_times',   'Ca Rise Time',           'feature'; ...
        'ca_tau',          'Ca Decay Constant',      'feature'; ...
        'apd50_data',      'APD50',                  'feature'; ...
        'ca50_data',       'CaTD50',                 'feature'; ...
        'complexity_data', 'Arrhythmia Complexity',  'complexity'; ...
        'alternans_data',  'Alternans',              'alternans'};

    for k = 1:size(pages, 1)
        fld = pages{k,1};
        if ~isfield(em, fld) || isempty(em.(fld)), continue; end
        md = em.(fld);
        for cam = 1:size(md, 1)
            if size(md, 2) < 1 || isempty(md{cam,1}), continue; end
            fig = [];
            try
                switch pages{k,3}
                    case 'feature',    fig = page_feature(md(cam,:), pages{k,2}, o.Pixel, name, cam);
                    case 'cv',         fig = page_cv(md(cam,:), pages{k,2}, fov, name);
                    case 'complexity', fig = page_complexity(md(cam,:), pages{k,2}, o.Pixel, pacing, name);
                    case 'alternans',  fig = page_alternans(md(cam,:), pages{k,2}, o.Pixel, pacing, name);
                end
                exportgraphics(fig, pdf_file, 'ContentType', 'image', 'Resolution', o.Resolution, ...
                               'Append', n_pages > 0);
                n_pages = n_pages + 1;
                if ~isempty(o.PNGDir)
                    png = fullfile(o.PNGDir, sprintf('%s_%02d_%s_CAM%d.png', base, n_pages, ...
                          regexprep(pages{k,2}, '[^A-Za-z0-9]+', '_'), cam));
                    exportgraphics(fig, png, 'Resolution', o.Resolution);
                end
            catch ME
                logf(sprintf('maps PDF: %s CAM%d page skipped (%s)', pages{k,2}, cam, ME.message));
            end
            if ~isempty(fig) && isvalid(fig), delete(fig); end
        end
    end
    if n_pages == 0
        logf('maps PDF: no pages written.');
        if isfile(pdf_file), delete(pdf_file); end
    end
end


% =========================================================================
%  Pages
% =========================================================================
function fig = page_feature(m, label, px, name, cam_index)
% {map, bg, trace, x_loc | curve, cam, unmasked}
    map = m{1,1};  bg = m{1,2};  tr = m{1,3};  xloc = m{1,4};  cam = m{1,5};
    row = px(1);  col = px(2);
    camlbl = char(string(cam));
    if isempty(regexp(camlbl, '^\d+$', 'once')), camlbl = sprintf('%d', cam_index); end   % 'VCDelay' slot
    ttl = sprintf('%s: CAM%s', label, camlbl);
    fig = new_page(1100, 800);
    t = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(t, sprintf('%s   |   %s', name, ttl), 'Interpreter', 'none', 'FontWeight', 'bold');

    ax1 = nexttile(t);
    draw_map(ax1, bg, map, robust_clim(map), [row col]);
    title(ax1, ttl, 'Interpreter', 'none');

    ax2 = nexttile(t);
    vals = map(isfinite(map));
    x = NaN;
    if row <= size(map,1) && col <= size(map,2), x = map(row, col); end
    if ~isempty(vals)
        histogram(ax2, vals);
        hold(ax2, 'on');
        if isfinite(x), xline(ax2, x, 'r-', 'LineWidth', 1.5); end
    end
    title(ax2, 'Distribution');  xlabel(ax2, label);  ylabel(ax2, 'Count');

    ax3 = nexttile(t);
    if strcmp(char(string(cam)), 'VCDelay') && iscell(tr)
        plot(ax3, tr{1,1}, 'LineWidth', 2.5);  hold(ax3, 'on');
        plot(ax3, tr{2,1}, 'LineWidth', 2.5, 'LineStyle', '-.', 'Color', 'r');
        legend(ax3, {'Vm', 'Ca'}, 'Location', 'best');
    elseif ~isempty(tr)
        plot(ax3, tr, 'LineWidth', 2.5);  hold(ax3, 'on');
        if isscalar(xloc) && isfinite(xloc) && round(xloc) >= 1 && round(xloc) <= numel(tr)
            plot(ax3, round(xloc), tr(round(xloc)), 'Marker', 'square', 'MarkerSize', 6, ...
                 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'LineStyle', 'none');
        elseif ~isscalar(xloc) && ~isempty(xloc)
            plot(ax3, xloc, 'LineWidth', 2.5, 'Color', 'r');
        end
    end
    title(ax3, 'Representative signal');  xlabel(ax3, 'Frame');

    ax4 = nexttile(t);
    stats_text(ax4, {'Mean', 'STD', 'Median', 'Sample'}, {ttl}, ...
        {sprintf('%.1f', mean(vals)); sprintf('%.1f', std(vals)); sprintf('%.1f', median(vals)); sprintf('%.1f', x)});
end

function fig = page_cv(m, label, fov, name)
% {cv_data {mag, Vx, Vy, ..., meta(6), T(7)}, bg, act map, cam}
    cv = m{1,1};  bg = m{1,2};  cam = m{1,4};
    cvmag = cv{1,1};
    Vx = [];  Vy = [];  meta = [];  T = [];
    if numel(cv) >= 3, Vx = cv{1,2};  Vy = cv{1,3}; end
    if numel(cv) >= 6, meta = cv{1,6}; end
    if numel(cv) >= 7, T = cv{1,7}; end
    if isempty(T), T = cvmag; end
    ttl = sprintf('%s: CAM%s', label, char(string(cam)));
    fig = new_page(1100, 700);
    t = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(t, sprintf('%s   |   %s', name, ttl), 'Interpreter', 'none', 'FontWeight', 'bold');

    % colour scale from interior pixels (8 px erosion), capped at 100 cm/s, as the app
    valid = isfinite(cvmag) & cvmag > 1 & cvmag < 249.9;
    try
        tissue = imfill(isfinite(cvmag), 'holes');
        interior = imerode(tissue, true(17));
        vi = valid & interior;
        if nnz(vi) >= 50, valid = vi; end
    catch
    end
    vals = cvmag(valid & cvmag > 0);
    if nnz(vals) < 50, vals = cvmag(isfinite(cvmag) & cvmag > 0); end
    if isempty(vals), cl = [0 1]; else, cl = [floor(prctile(vals, 2)), min(100, ceil(prctile(vals, 98)))]; end
    if cl(2) <= cl(1), cl(2) = cl(1) + 1; end

    ax1 = nexttile(t);
    draw_map(ax1, bg, cvmag, cl, []);
    cb = colorbar(ax1, 'southoutside');  cb.Label.String = 'CV (cm/s)';
    title(ax1, ttl, 'Interpreter', 'none');
    if ~isempty(Vx) && ~isempty(Vy) && isequal(size(Vx), size(T))
        if isfinite(fov) && fov > 0, px_mm = fov / size(T,1); else, px_mm = 20 / size(T,1); end
        step = max(2, round(1.0 / max(px_mm, eps)));
        [xx, yy] = meshgrid(1:size(T,2), 1:size(T,1));
        if isstruct(meta) && isfield(meta, 'good') && ~isempty(meta.good), good = meta.good;
        else, good = isfinite(Vx) & isfinite(Vy); end
        mask = good & mod(xx, step) == 0 & mod(yy, step) == 0;
        mag = hypot(Vx, Vy);
        vref = prctile(mag(mask), 90);
        if ~isfinite(vref) || vref <= 0, vref = 1; end
        hold(ax1, 'on');
        quiver(ax1, xx(mask), yy(mask), Vx(mask) ./ vref, Vy(mask) ./ vref, 0.9, 'm', 'LineWidth', 0.5);
    end

    ax2 = nexttile(t);
    allv = cvmag(isfinite(cvmag));
    if ~isempty(allv), histogram(ax2, allv, 40); end
    xlabel(ax2, 'Speed (cm/s)');  ylabel(ax2, 'Count');
    title(ax2, sprintf('Conduction velocity distribution (median %.1f cm/s, interior)', median(vals)));
end

function fig = page_complexity(m, label, px, pacing, name)
% {{DF, RI, OI}, bg, trace, samples [DF RI OI], cam}
    maps = m{1,1};  bg = m{1,2};  tr = m{1,3};  xs = m{1,4};  cam = m{1,5};
    row = px(1);  col = px(2);
    ttl = sprintf('%s: CAM%s', label, char(string(cam)));
    titles = {'Dominant Frequency (DF) Hz', 'Regularity Index (RI)', 'Organization Index (OI)'};
    fig = new_page(1200, 880);
    t = tiledlayout(fig, 3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(t, sprintf('%s   |   %s', name, ttl), 'Interpreter', 'none', 'FontWeight', 'bold');
    st = nan(3, 4);
    for k = 1:3
        ax = nexttile(t, k);
        dmap = maps{1,k};
        if k == 1, cl = robust_clim(dmap); else, cl = [0 1]; end
        draw_map(ax, bg, dmap, cl, [row col]);
        title(ax, titles{k});
        v = dmap(isfinite(dmap));
        st(k,:) = [mean(v), std(v), median(v), sample_or(xs, k)];
    end
    for k = 1:3
        ax = nexttile(t, 3 + k);
        dmap = maps{1,k};  v = dmap(isfinite(dmap));
        if ~isempty(v), histogram(ax, v);  hold(ax, 'on'); end
        if isfinite(st(k,4)), xline(ax, st(k,4), 'r-', 'LineWidth', 1.5); end
        xlabel(ax, titles{k});  ylabel(ax, 'Count');
    end
    ax7 = nexttile(t, 7, [1 2]);
    trace_with_pacing(ax7, tr, pacing);
    ax9 = nexttile(t, 9);
    stats_text(ax9, {'Mean', 'STD', 'Median', 'Sample'}, {'DF', 'RI', 'OI'}, fmt_matrix(st', '%.2f'));
end

function fig = page_alternans(m, label, px, pacing, name)
% {alternans struct, bg, trace, cam, substrate, 'AP' | 'Ca'}
    a = m{1,1};  bg = m{1,2};  tr = m{1,3};  cam = m{1,4};
    sub = [];  if size(m, 2) >= 5, sub = m{1,5}; end
    branch = 'AP';  if size(m, 2) >= 6 && ~isempty(m{1,6}), branch = char(string(m{1,6})); end
    row = px(1);  col = px(2);
    ttl = sprintf('%s: CAM%s', label, char(string(cam)));
    if strcmp(branch, 'AP')
        panel = {sfield(sub, 'gradient_mag'), sfield(sub, 'restitution_map'), sfield(a, 'APD50_ratio_map'), ...
                 sfield(a, 'dvdt_ratio_map'), sfield(sub, 'risk_map'), sfield(a, 'phase_map')};
        ptitle = {'APD50 Spatial Gradient', 'APD50 Restitution Slope', 'APD50 Alternans Ratio Map', ...
                  'Conduction (dv/dt ratio map)', 'AP Alternans Risk Map', 'Alternans Phase Map'};
        rows = {'APD50 Spatial Gradient'; 'APD50 Restitution Slope'; 'APD50 Alt Ratio'; ...
                'Conduction'; 'AP Amp Ratio'; 'APD50 Phase'};
    else
        panel = {sfield(sub, 'coupling_map'), sfield(a, 'release_alt_map'), sfield(a, 'load_alt_map'), ...
                 sfield(a, 'D50_alt_map'), sfield(a, 'D50_ratio_map'), sfield(a, 'phase_map')};
        ptitle = {'Ca-AP Pearson r map', 'SR Ca Release Alternans Map', 'SR Ca Load Alternans Map', ...
                  'CaTD50 Alternans Map', 'CaTD50 Alternans Ratio Map', 'CaT Amplitude Alternans Phase Map'};
        rows = {'AP-Ca Coupling'; 'SR Ca Release'; 'SR Ca Load'; 'CaD50 Alt'; 'CaD50 Alt Ratio'; 'CaT Amp Phase'};
    end
    fig = new_page(1200, 880);
    t = tiledlayout(fig, 3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(t, sprintf('%s   |   %s', name, ttl), 'Interpreter', 'none', 'FontWeight', 'bold');
    st = nan(6, 4);
    for k = 1:6
        ax = nexttile(t, k);
        dmap = panel{k};
        if isempty(dmap)
            draw_map(ax, bg, [], [], [row col]);
            title(ax, [ptitle{k} ' (n/a)']);
            continue;
        end
        draw_map(ax, bg, dmap, robust_clim(dmap), [row col]);
        title(ax, ptitle{k});
        v = dmap(isfinite(dmap));
        smp = NaN;
        if row <= size(dmap,1) && col <= size(dmap,2), smp = dmap(row, col); end
        if ~isempty(v), st(k,:) = [mean(v), std(v), median(v), smp]; end
    end
    ax7 = nexttile(t, 7, [1 2]);
    trace_with_pacing(ax7, tr, pacing);
    ax9 = nexttile(t, 9);
    stats_text(ax9, rows, {'Mean', 'STD', 'Median', 'Sample'}, fmt_matrix(st, '%.3g'));
end


% =========================================================================
%  Drawing helpers
% =========================================================================
function fig = new_page(w, h)
    fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', 'Position', [50 50 w h], ...
                 'InvertHardcopy', 'off');
end

function draw_map(ax, bg, map, cl, marker)
    if ~isempty(bg)
        image(ax, real2rgb(double(bg), 'gray'));
        hold(ax, 'on');
    end
    if ~isempty(map)
        imagesc(ax, map, 'AlphaData', isfinite(map));
        hold(ax, 'on');
        colormap(ax, 'jet');
        if ~isempty(cl) && all(isfinite(cl)) && cl(2) > cl(1), clim(ax, cl); end
        colorbar(ax, 'southoutside');
    end
    if ~isempty(marker)
        plot(ax, marker(2), marker(1), 'Marker', 'square', 'MarkerSize', 8, ...
             'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'r', 'LineStyle', 'none');
    end
    set(ax, 'XTick', [], 'YTick', [], 'YDir', 'reverse');
    axis(ax, 'image');
end

function cl = robust_clim(map)
% median +/- 3 robust-sigma (MAD-based), as the app's display functions
    vals = map(isfinite(map));
    if isempty(vals), cl = []; return; end
    med  = median(vals);
    rsig = 1.4826 * median(abs(vals - med));
    if rsig > 0
        cl = [floor(med - 3*rsig), ceil(med + 3*rsig)];
    else
        cl = [min(vals), max(vals)];
    end
    if cl(2) <= cl(1), cl = [cl(1) - 0.5, cl(1) + 0.5]; end
end

function trace_with_pacing(ax, tr, pacing)
    if ~isempty(tr)
        plot(ax, tr, 'LineWidth', 1.5);
        hold(ax, 'on');
    end
    if ~isempty(pacing)
        plot(ax, 1:numel(pacing), pacing, 'LineWidth', 1, 'Color', [0.47 0.47 0.47]);
    end
    title(ax, 'Representative signal (grey: pacing)');  xlabel(ax, 'Frame');
end

function stats_text(ax, row_names, col_names, cells)
% A small text table in place of the app's uitable (exportgraphics ignores
% uicontrols).  cells is {rows x cols} of strings.
    axis(ax, 'off');
    nr = numel(row_names);  nc = numel(col_names);
    x0 = 0.02;                                   % row-name column
    xv = 0.42 + (0:nc-1) * (0.56 / nc);          % value columns (names get ~40 % of the width)
    y  = linspace(0.9, 0.1, nr + 1);
    fs = 9;  if nr > 4, fs = 8; end
    for c = 1:nc
        text(ax, xv(c), y(1), col_names{c}, 'Units', 'normalized', 'FontWeight', 'bold', ...
             'Interpreter', 'none', 'FontSize', fs, 'HorizontalAlignment', 'left');
    end
    for r = 1:nr
        text(ax, x0, y(r+1), row_names{r}, 'Units', 'normalized', 'FontWeight', 'bold', ...
             'Interpreter', 'none', 'FontSize', fs);
        for c = 1:nc
            text(ax, xv(c), y(r+1), cells{r,c}, 'Units', 'normalized', 'Interpreter', 'none', 'FontSize', fs);
        end
    end
    title(ax, 'Statistics');
end

function C = fmt_matrix(M, fmt)
    C = cell(size(M));
    for i = 1:numel(M)
        if isfinite(M(i)), C{i} = sprintf(fmt, M(i)); else, C{i} = 'n/a'; end
    end
end

function v = sfield(s, f)
    v = [];
    if isstruct(s) && isfield(s, f), v = s.(f); end
end

function v = sample_or(xs, k)
    v = NaN;
    if isnumeric(xs) && numel(xs) >= k, v = xs(k); end
end
