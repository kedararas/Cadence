function save_figures_to_pdf(figs, filename, closeAfterSave)
% SAVE_FIGURES_TO_PDF  Save an array of uifigure handles into one PDF file.
%
%   save_figures_to_pdf(figs, filename)
%   save_figures_to_pdf(figs, filename, closeAfterSave)
%
%   figs           - array of uifigure handles, e.g. [fig1, fig2, fig3]
%   filename       - output path, e.g. 'results.pdf'
%   closeAfterSave - (optional, default true) close/delete the figures once
%                    saved. Batch processing passes true so windows do not pile
%                    up across files; single-file extraction passes false so the
%                    user can read/copy on-screen values (e.g. median APD80)
%                    before closing the windows manually.
%
%   Notes
%   -----
%   - Each window is captured with EXPORTAPP, which renders the whole uifigure
%     exactly as displayed — including uiaxes, colorbars, and uitable. This
%     replaces the old exportgraphics path: exportgraphics cannot capture a
%     multi-panel uifigure (it errors with "Figure has more than one
%     container..."), and it could not capture uitable at all.
%   - exportapp writes a single image per call and has no append option, so
%     each capture is assembled into the multi-page PDF via a hidden image
%     figure + exportgraphics(...,'Append',true).
%   - Output pages are raster (exportapp is a window screenshot). This is the
%     expected behavior for capturing app-style figures.
%
%   IMPORTANT: figs must be UIFIGURES (created with uifigure). exportapp does
%   not support classic figures created with figure(); if you still produce
%   any classic figures for export, they will be skipped with a warning.
%
%   Example:
%     f1 = uifigure; ax = uiaxes(f1); plot(ax, rand(10,1));
%     save_figures_to_pdf(f1, 'my_report.pdf');

if nargin < 2 || isempty(filename)
    filename = 'figures.pdf';
end
if nargin < 3 || isempty(closeAfterSave)
    closeAfterSave = true;   % preserve prior behavior for any other caller
end

% Ensure .pdf extension
[folder, name, ext] = fileparts(filename);
if ~strcmpi(ext, '.pdf')
    filename = fullfile(folder, [name '.pdf']);
end

% Delete existing file — exportgraphics appends, so stale pages would remain
if isfile(filename)
    delete(filename);
end

first_page = true;   % first page passes Append=false, the rest Append=true

for i = 1:numel(figs)
    fig = figs(i);
    if ~isvalid(fig)
        warning('save_figures_to_pdf: figure %d is not a valid handle — skipping.', i);
        continue;
    end

    drawnow;   % flush so the window is fully rendered before capture

    % ---- capture the whole window to a temp image ---- %
    tmp = [tempname '.png'];
    try
        exportapp(fig, tmp);
    catch ME
        warning(['save_figures_to_pdf: exportapp failed on figure %d (%s).\n' ...
                 'Is it a uifigure? Classic figures are not supported here. Skipping.'], ...
                 i, ME.message);
        if isfile(tmp), delete(tmp); end
        continue;
    end

    % ---- place the captured image on a PDF page ---- %
    try
        img = imread(tmp);
        pf  = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
                     'Position', [100 100 size(img,2) size(img,1)]);
        ax  = axes(pf, 'Units', 'normalized', 'Position', [0 0 1 1]);
        image(ax, img);
        axis(ax, 'image', 'off');

        exportgraphics(pf, filename, ...
            'ContentType', 'image', ...
            'Resolution',  200, ...
            'Append',      ~first_page);
        first_page = false;
        delete(pf);
    catch ME
        warning('save_figures_to_pdf: failed to add figure %d to PDF (%s).', i, ME.message);
        if exist('pf', 'var') && isvalid(pf), delete(pf); end
    end

    if isfile(tmp), delete(tmp); end
end

if first_page
    warning('save_figures_to_pdf: no pages were written (no valid uifigures exported).');
else
    fprintf('Saved figures → %s\n', filename);
end

% Close all figures now that they have been saved — but only when the caller
% asked for it (batch runs). Single-file extraction passes closeAfterSave=false
% so the windows stay open for the user to read/copy on-screen values.
% Use delete, not close: close() routes through CloseRequestFcn / handle
% visibility and frequently leaves uifigures open; delete() destroys them
% unconditionally.
if closeAfterSave
    delete(figs(isvalid(figs)));
end
end
