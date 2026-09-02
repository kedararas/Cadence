function results = mv_lat_quantization(varargin)
%MV_LAT_QUANTIZATION  What integer-frame activation times cost, per geometry.
%
%   results = mv_lat_quantization
%   results = mv_lat_quantization('OutDir', 'figures')
%
%   CADENCE times activation by linear interpolation of the 50% upstroke
%   crossing (compute_lat_50) rather than by taking the max-dV/dt frame.  This
%   script measures what that buys, by running the CV engines on the same
%   analytic activation map twice: once rounded to whole frames, once with the
%   sub-frame precision the pipeline actually achieves.
%
%   WHY A JITTER SWEEP DOES NOT ANSWER THIS
%     Rounding to 1 kHz frames gives a standard deviation of 0.289 ms (error
%     uniform on +/-0.5 ms), so it is tempting to read the cost off the
%     mv_synth_cv jitter curve at sigma = 0.289.  That understates it badly.
%     Quantization error is STRUCTURED, not random: rounding a smooth
%     activation ramp produces flat plateaus separated by one-frame risers.
%     Inside a plateau grad(T) = 0, so CV diverges and the pixel is rejected by
%     the velocity cap; all the gradient piles up at the risers.  On a 10 mm
%     slice ~15 adjacent pixels share one frame, and the result is ~6x worse
%     than random noise of the same magnitude while losing a third of the map.
%
%     The damage therefore scales with pixel fineness, and vanishes once the
%     adjacent-pixel LAT step (ds/v) approaches the frame period: sub-frame
%     interpolation is ENABLING for slice and rat whole heart, and close to
%     cosmetic for rabbit whole heart and wedge.  Reporting it as a single
%     number across preparations would be wrong.
%
%   OUTPUTS WRITTEN (into OutDir)
%     cv_lat_quantization.csv     one row per preparation x engine
%     fig_lat_staircase.png       LAT maps, CV maps and profiles for one
%                                 preparation (default the finest), showing the
%                                 staircase and the pixels it costs
%
%   Name-value
%     'FOV'        field of view per preparation, mm (default [10 20 30 50])
%     'Labels'     display names (default slice / rat WH / rabbit WH / wedge)
%     'GridSize'   pixels per side (default 256)
%     'Velocity'   planar wave speed, cm/s (default 60)
%     'SubFrameSD' sigma_LAT achieved by interpolated LAT, ms (default 0.09,
%                  measured end-to-end at SNR 20 in mv_synth_cv test E)
%     'Engines'    any of {'gradient','polyfit','circle'} (default all three)
%     'FigurePrep' index into FOV to illustrate (default 1, the finest)
%     'OutDir'     where to write table and figure (default pwd)
%     'Save'       write files (default true)
%
%   Output
%     results  table: preparation, engine, integer-frame vs sub-frame median
%              |error| and pixels retained, plus the measured rounding sigma

    p = inputParser;
    p.addParameter('FOV',        [10 20 30 50], @isnumeric);
    p.addParameter('Labels',     {}, @iscell);
    p.addParameter('GridSize',   256, @isscalar);
    p.addParameter('Velocity',   60,  @isscalar);
    p.addParameter('SubFrameSD', 0.09, @isscalar);
    p.addParameter('Engines',    {'gradient','polyfit','circle'}, @iscell);
    p.addParameter('FigurePrep', 1, @isscalar);
    p.addParameter('EdgeMargin', 12, @isscalar);
    p.addParameter('OutDir',     pwd, @(s) ischar(s) || isstring(s));
    p.addParameter('Save',       true, @islogical);
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

    outdir = char(o.OutDir);
    if o.Save && ~isempty(outdir) && ~exist(outdir, 'dir'), mkdir(outdir); end

    n  = o.GridSize;
    v  = o.Velocity / 100;                       % mm/ms
    em = o.EdgeMargin;
    interior = false(n, n);
    interior(em+1:n-em, em+1:n-em) = true;

    [prep, eng, e_int, k_int, e_sub, k_sub, sd_q, step_ms] = deal(...
        strings(0,1), strings(0,1), [], [], [], [], [], []);

    fprintf('\n%-18s %-9s | %-22s | %-22s\n', 'preparation', 'engine', ...
            'INTEGER-FRAME', sprintf('SUB-FRAME (sd %.2f ms)', o.SubFrameSD));
    fprintf('%s\n', repmat('-', 1, 78));

    rng(7, 'twister');
    for i = 1:numel(o.FOV)
        ds = o.FOV(i) / n;
        [cc, ~] = meshgrid(0:n-1, 0:n-1);
        Ttrue = (cc * ds) / v;  Ttrue = Ttrue - min(Ttrue(:));

        Tint = round(Ttrue);                                 % integer-frame LAT
        Tsub = Ttrue + o.SubFrameSD * randn(n, n);           % interpolated LAT
        sdq  = std(Tint(interior) - Ttrue(interior));

        for e = 1:numel(o.Engines)
            [ei, ki] = score(Tint, ds, o.Engines{e}, interior, o.Velocity);
            [es, ks] = score(Tsub, ds, o.Engines{e}, interior, o.Velocity);

            prep(end+1,1)   = string(labels{i});             %#ok<AGROW>
            eng(end+1,1)    = string(o.Engines{e});          %#ok<AGROW>
            e_int(end+1,1)  = ei;   k_int(end+1,1) = ki;     %#ok<AGROW>
            e_sub(end+1,1)  = es;   k_sub(end+1,1) = ks;     %#ok<AGROW>
            sd_q(end+1,1)   = sdq;                           %#ok<AGROW>
            step_ms(end+1,1)= ds / v;                        %#ok<AGROW>

            if e == 1, nm = labels{i}; else, nm = ''; end
            fprintf('%-18s %-9s | %6.2f%% kept %5.1f%%   | %6.2f%% kept %5.1f%%\n', ...
                    nm, o.Engines{e}, ei, ki, es, ks);
        end
        fprintf('%18s (rounding sigma %.3f ms; adjacent-pixel LAT step %.3f ms)\n', ...
                '', sdq, ds/v);
    end

    results = table(prep, eng, step_ms, sd_q, e_int, k_int, e_sub, k_sub, ...
        'VariableNames', {'preparation','engine','lat_step_ms','rounding_sd_ms', ...
                          'int_frame_err_pct','int_frame_kept_pct', ...
                          'sub_frame_err_pct','sub_frame_kept_pct'});

    if ~o.Save, return; end
    writetable(results, fullfile(outdir, 'cv_lat_quantization.csv'));
    plot_staircase(o, labels, interior, outdir);
    fprintf('\nWrote table and figure to %s\n', outdir);
end


% ======================================================================
function [err, kept] = score(T, ds, eng, interior, v_cms)
%SCORE  Median |CV error| and surviving-pixel fraction for one activation map.
%   Pixels the engine rejects (CV above the physiological cap) are excluded
%   from the error, which is exactly why `kept` must be reported alongside it:
%   a median over survivors is not an accuracy figure on its own.
    if strcmpi(eng, 'circle')
        [Vmag, ~, ~, ~, ~, ~] = circle_method_cv(T, ds, ds, 10);
    else
        [Vmag, ~, ~, ~, ~, ~] = conduction_velocity(T, ds, ds, struct('method', eng));
    end
    cv = Vmag * 100;
    m  = interior & isfinite(cv);
    kept = 100 * sum(m(:)) / sum(interior(:));
    if ~any(m(:)), err = NaN; else, err = median(abs(100*(cv(m) - v_cms)/v_cms)); end
end


function plot_staircase(o, labels, interior, outdir)
%PLOT_STAIRCASE  Show the mechanism, not just the penalty.
%   Six panels: activation maps and a profile on top, the CV fields they
%   produce and a CV profile below.  The profiles are what make the argument --
%   the LAT staircase is obvious, and the CV trace swings between far too fast
%   and far too slow across every tread.
    i  = min(max(round(o.FigurePrep), 1), numel(o.FOV));
    n  = o.GridSize;
    ds = o.FOV(i) / n;
    v  = o.Velocity / 100;

    [cc, ~] = meshgrid(0:n-1, 0:n-1);
    Ttrue = (cc * ds) / v;  Ttrue = Ttrue - min(Ttrue(:));
    rng(7, 'twister');
    Tsub = Ttrue + o.SubFrameSD * randn(n, n);
    Tint = round(Ttrue);

    [Vs, ~, ~, ~, ~, ~] = conduction_velocity(Tsub, ds, ds, struct('method','polyfit'));
    [Vi, ~, ~, ~, ~, ~] = conduction_velocity(Tint, ds, ds, struct('method','polyfit'));
    CVs = Vs * 100;  CVi = Vi * 100;
    [es, ks] = score(Tsub, ds, 'polyfit', interior, o.Velocity);
    [ei, ki] = score(Tint, ds, 'polyfit', interior, o.Velocity);

    xmm = (0:n-1) * ds;  row = round(n/2);
    xzoom = [0 min(4, max(xmm))];
    f = figure('Position', [80 80 1150 640], 'Color', 'w');

    ax = subplot(2,3,1); imagesc(xmm, xmm, Tsub); axis image; colormap(ax, parula(256));
    caxis([0 max(Ttrue(:))]); c = colorbar; c.Label.String = 'LAT (ms)';
    title({'Sub-frame LAT (50% crossing)', sprintf('\\sigma_{LAT} = %.2f ms', o.SubFrameSD)}, ...
          'FontWeight','normal','FontSize',10);
    xlabel('mm'); ylabel('mm'); set(gca,'FontSize',9);

    ax = subplot(2,3,2); imagesc(xmm, xmm, Tint); axis image; colormap(ax, parula(256));
    caxis([0 max(Ttrue(:))]); c = colorbar; c.Label.String = 'LAT (ms)';
    title({'Integer-frame LAT', 'quantized to whole frames'}, ...
          'FontWeight','normal','FontSize',10);
    xlabel('mm'); ylabel('mm'); set(gca,'FontSize',9);

    subplot(2,3,3); hold on; box on
    plot(xmm, Tsub(row,:), '-', 'Color',[0.00 0.45 0.74], 'LineWidth',1.2);
    plot(xmm, Tint(row,:), '-', 'Color',[0.85 0.33 0.10], 'LineWidth',1.4);
    plot(xmm, Ttrue(row,:), 'k--', 'LineWidth',0.8);
    xlim(xzoom); xlabel('mm'); ylabel('LAT (ms)');
    title('Activation profile','FontWeight','normal','FontSize',10);
    legend({'sub-frame','integer-frame','truth'},'Location','northwest','Box','off','FontSize',8);
    set(gca,'FontSize',9);

    % White marks rejected pixels: give NaN its own slot at the bottom.
    cmapV = [1 1 1; parula(255)];
    Ds = CVs; Ds(~isfinite(Ds)) = -1;
    Di = CVi; Di(~isfinite(Di)) = -1;
    hi = 2 * o.Velocity;

    ax = subplot(2,3,4); imagesc(xmm, xmm, Ds, [-1 hi]); axis image; colormap(ax, cmapV);
    c = colorbar; c.Label.String = 'CV (cm/s)';
    title(sprintf('CV from sub-frame LAT\nerror %.1f%%, %.0f%% pixels kept', es, ks), ...
          'FontWeight','normal','FontSize',10);
    xlabel('mm'); ylabel('mm'); set(gca,'FontSize',9);

    ax = subplot(2,3,5); imagesc(xmm, xmm, Di, [-1 hi]); axis image; colormap(ax, cmapV);
    c = colorbar; c.Label.String = 'CV (cm/s)';
    title(sprintf('CV from integer-frame LAT\nerror %.1f%%, %.0f%% pixels kept', ei, ki), ...
          'FontWeight','normal','FontSize',10);
    xlabel('mm'); ylabel('mm'); set(gca,'FontSize',9);

    subplot(2,3,6); hold on; box on
    plot(xmm, CVs(row,:), '-', 'Color',[0.00 0.45 0.74], 'LineWidth',1.2);
    plot(xmm, CVi(row,:), '-', 'Color',[0.85 0.33 0.10], 'LineWidth',1.4);
    yline(o.Velocity, 'k--', 'LineWidth',0.9);
    xlim(xzoom); ylim([0 hi + o.Velocity]); xlabel('mm'); ylabel('CV (cm/s)');
    title('CV profile (gaps = rejected)','FontWeight','normal','FontSize',10);
    legend({'sub-frame','integer-frame',sprintf('true %g cm/s', o.Velocity)}, ...
           'Location','northwest','Box','off','FontSize',8);
    set(gca,'FontSize',9);

    exportgraphics(f, fullfile(outdir, 'fig_lat_staircase.png'), 'Resolution', 300);
end


function ensure_path()
%ENSURE_PATH  Add sibling helper folders if the CV engines are not visible.
    if exist('conduction_velocity', 'file') == 2, return; end
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    for f = {'conduction_velocity_helper', 'feature_extraction_helper', ...
             'signal_conditioning_helper'}
        d = fullfile(root, f{1});
        if exist(d, 'dir'), addpath(d); end
    end
end
