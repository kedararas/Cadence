function marks_file = mv_mark(manifest_file, reviewer_id, varargin)
%MV_MARK  Blinded manual fiducial marking for software validation.
%
%   marks_file = mv_mark(manifest_file, reviewer_id)
%   marks_file = mv_mark(..., 'Mode', 'alternans', 'Percent', 50)
%
%   Presents one pixel's optical trace at a time and records where a human
%   places the fiducial points.  It computes no metric and detects no feature:
%   the reviewer supplies every judgement, and mv_derive does the arithmetic
%   afterwards from the saved clicks.
%
%   WHAT IS BLINDED, AND WHY
%   This function never loads a *-metrics.mat and never displays a CADENCE
%   value.  That is the point.  If the automated answer is visible — even as a
%   marker "for reference" — reviewers anchor to it and the agreement you
%   measure is manufactured.  Mark first; compare in mv_compare, afterwards.
%   Pixel identity (row/col) is also withheld, so a reviewer who has seen the
%   APD map for this recording cannot place the pixel on it from memory.
%
%   MODES
%     'apd' (default)  Four clicks per pixel: diastolic baseline, peak,
%                      activation, and the crossing at 'Percent' repolarization.
%                      After baseline and peak are set, a guide line is drawn at
%                      the repolarization level so the reviewer marks WHERE the
%                      trace crosses it rather than eyeballing a percentage.
%                      Yields activation time, APD/CaTD at 'Percent', amplitude.
%     'alternans'      Three clicks: baseline, peak of beat 1, peak of beat 2.
%                      Yields the two amplitudes and their ratio.  Build the
%                      manifest with 'NumBeats', 2.
%
%   HOW A CLICK IS INTERPRETED
%   The reviewer's click fixes a TIME.  The harness then reads the trace value
%   at that time — it does not search for a nearby maximum, which would be
%   exactly the automated detection under test.  The one exception is the
%   baseline, taken as the median of the trace over +/-'BaselineWin' ms around
%   the clicked time, because a single diastolic sample is noise-dominated.
%   Both behaviours are recorded in the output so they can be reported.
%
%   KEYS (while marking)
%     u  undo last click      r  restart this pixel     s  skip (unmarkable)
%     b  back one pixel       q  save and quit
%
%   Progress is written to disk after every pixel, so the session is resumable:
%   re-run with the same reviewer_id and marking continues where it stopped.
%
%   Inputs
%     manifest_file  output of mv_sample_pixels
%     reviewer_id    short string, e.g. "AA" — identifies the reviewer and seeds
%                    their presentation order
%
%   Name-value
%     'Mode'         'apd' (default) | 'alternans'
%     'Percent'      repolarization percentage for 'apd' (default 80)
%     'BaselineWin'  half-width in ms for the baseline median (default 5)
%     'Output'       marks path (default alongside the manifest)

    p = inputParser;
    p.addParameter('Mode',        'apd', @(v) any(strcmpi(v, {'apd','alternans'})));
    p.addParameter('Percent',     80,    @(v) isnumeric(v) && isscalar(v) && v > 0 && v < 100);
    p.addParameter('BaselineWin', 5,     @(v) isnumeric(v) && isscalar(v) && v > 0);
    p.addParameter('Output',      '',    @(v) ischar(v) || isstring(v));
    p.parse(varargin{:});
    o = p.Results;
    o.Mode = lower(o.Mode);

    reviewer = char(string(reviewer_id));
    if isempty(reviewer)
        error('mv_mark:reviewer', 'reviewer_id must be a non-empty string.');
    end

    M = load(manifest_file);
    if ~isfield(M, 'manifest')
        error('mv_mark:badManifest', '%s does not contain a ''manifest'' struct.', manifest_file);
    end
    manifest = M.manifest;

    if strcmp(o.Mode, 'alternans') && manifest.num_beats < 2
        warning('mv_mark:oneBeat', ...
                ['Manifest window covers %d beat(s) but mode is ''alternans''. ' ...
                 'Re-run mv_sample_pixels with ''NumBeats'', 2.'], manifest.num_beats);
    end

    % ---- trace data for the sampled pixels, restricted to the marking window ----
    S = load(manifest.conditioned_file);
    if isfield(S, 'cmos_all_data'); d = S.cmos_all_data; else
        fn = fieldnames(S); is_struct = structfun(@isstruct, S);
        d = S.(fn{find(is_struct, 1)});
    end
    X  = d.(manifest.cam);
    nt = size(X, 3);
    wf = manifest.window_frames;
    T  = double(reshape(X, [], nt));
    T  = T(manifest.linear_index, wf(1):wf(2));            % nPix x nFrames
    t_ms = ((wf(1):wf(2)) - 1) / manifest.acqFreq * 1000;  % absolute ms
    npix = size(T, 1);

    % ---- resume or start ----
    if isempty(o.Output)
        [pdir, base] = fileparts(manifest_file);
        o.Output = fullfile(pdir, sprintf('%s_%s_marks.mat', base, reviewer));
    end
    o.Output = char(o.Output);

    if isfile(o.Output)
        R = load(o.Output);
        marks = R.marks;
        if ~strcmp(marks.mode, o.Mode) || marks.percent ~= o.Percent
            error('mv_mark:resumeMismatch', ...
                  ['%s was started with Mode=''%s'', Percent=%g but you passed ' ...
                   'Mode=''%s'', Percent=%g. Use the original settings or a new Output file.'], ...
                  o.Output, marks.mode, marks.percent, o.Mode, o.Percent);
        end
        fprintf('Resuming %s: %d of %d pixels already marked.\n', ...
                o.Output, sum(marks.done | marks.skipped), npix);
    else
        marks = new_marks(manifest, manifest_file, reviewer, o, npix);
        fprintf('New marking session for reviewer ''%s'': %d pixels.\n', reviewer, npix);
    end

    nclick = click_count(o.Mode);
    fig = figure('Name', sprintf('mv_mark — reviewer %s [BLINDED]', reviewer), ...
                 'NumberTitle', 'off', 'Color', 'w');

    % No onCleanup guard here on purpose.  onCleanup would capture `marks` by
    % VALUE at creation, so on an error it would write that stale pre-loop
    % snapshot over the file — destroying the very progress it looks like it is
    % protecting.  The loop below saves after every pixel instead, which is both
    % safer and what makes the session resumable.

    % ---- main loop over the reviewer's presentation order ----
    k = find(~(marks.done | marks.skipped), 1);
    if isempty(k); k = npix + 1; end

    while k >= 1 && k <= npix
        if ~ishandle(fig)
            fprintf('Figure closed — saving and exiting.\n'); break;
        end
        i = marks.order(k);                     % index into manifest.pixels
        y = T(i, :);

        [res, action] = mark_one(fig, t_ms, y, k, npix, o, nclick);

        switch action
            case 'done'
                marks = store(marks, i, res, o, t_ms, y);
                marks.done(i) = true; marks.skipped(i) = false;
                k = k + 1;
            case 'skip'
                marks.skipped(i) = true; marks.done(i) = false;
                marks.clicks{i}  = zeros(0, 2);
                k = k + 1;
            case 'back'
                k = max(1, k - 1);
            case 'quit'
                save_marks(o.Output, marks);
                fprintf('Saved %d marked / %d skipped to %s\n', ...
                        sum(marks.done), sum(marks.skipped), o.Output);
                if ishandle(fig); close(fig); end
                marks_file = o.Output;
                return;
        end
        save_marks(o.Output, marks);            % resumable after every pixel
    end

    if ishandle(fig); close(fig); end
    save_marks(o.Output, marks);
    fprintf('Complete: %d marked, %d skipped -> %s\n', ...
            sum(marks.done), sum(marks.skipped), o.Output);
    fprintf('Next: mv_derive(''%s'')\n', o.Output);
    marks_file = o.Output;
end


% ======================================================================
function marks = new_marks(manifest, manifest_file, reviewer, o, npix)
    % Presentation order is randomized PER REVIEWER: it removes order/fatigue
    % effects from the inter-observer comparison and makes it impractical for
    % two reviewers to sit together and mark "the same one" in step.
    seed = mod(sum(double(reviewer)) * 7919 + manifest.seed, 2^31 - 1);
    rng(seed, 'twister');

    marks = struct();
    marks.manifest_file   = char(manifest_file);
    marks.conditioned_file = manifest.conditioned_file;
    marks.cam             = manifest.cam;
    marks.reviewer        = reviewer;
    marks.mode            = o.Mode;
    marks.percent         = o.Percent;
    marks.baseline_win_ms = o.BaselineWin;
    marks.order_seed      = seed;
    marks.order           = randperm(npix)';
    marks.done            = false(npix, 1);
    marks.skipped         = false(npix, 1);
    marks.clicks          = repmat({zeros(0, 2)}, npix, 1);
    marks.seconds         = nan(npix, 1);
    marks.t_base_ms       = nan(npix, 1);
    marks.v_base          = nan(npix, 1);
    marks.t_peak_ms       = nan(npix, 1);
    marks.v_peak          = nan(npix, 1);
    marks.t_act_ms        = nan(npix, 1);
    marks.t_rep_ms        = nan(npix, 1);
    marks.t_peak2_ms      = nan(npix, 1);
    marks.v_peak2         = nan(npix, 1);
    marks.created         = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
end


function n = click_count(mode)
    switch mode
        case 'apd',       n = 4;   % baseline, peak, activation, repolarization
        case 'alternans', n = 3;   % baseline, peak 1, peak 2
    end
end


function lbl = click_label(mode, j, pct)
    if strcmp(mode, 'apd')
        switch j
            case 1, lbl = 'DIASTOLIC BASELINE (click in diastole, before the upstroke)';
            case 2, lbl = 'PEAK (click at the time of the peak)';
            case 3, lbl = 'ACTIVATION (click on the upstroke)';
            case 4, lbl = sprintf('REPOLARIZATION — click where the trace crosses the %g%% guide line', pct);
        end
    else
        switch j
            case 1, lbl = 'DIASTOLIC BASELINE (click in diastole, before beat 1)';
            case 2, lbl = 'PEAK of BEAT 1';
            case 3, lbl = 'PEAK of BEAT 2';
        end
    end
end


function [res, action] = mark_one(fig, t_ms, y, k, npix, o, nclick)
%MARK_ONE  Collect one pixel's clicks.  Returns raw click coordinates.

    res = zeros(0, 2);
    t0  = tic;

    while true
        if ~ishandle(fig); action = 'quit'; return; end
        figure(fig); clf(fig);
        ax = axes('Parent', fig);
        plot(ax, t_ms, y, '-', 'Color', [0.20 0.35 0.75], 'LineWidth', 1.2);
        hold(ax, 'on'); grid(ax, 'on');
        xlabel(ax, 'Time (ms)'); ylabel(ax, 'Fluorescence (a.u.)');
        xlim(ax, [t_ms(1) t_ms(end)]);

        j = size(res, 1) + 1;

        % Draw what has been placed so far.
        for m = 1:size(res, 1)
            plot(ax, res(m,1), res(m,2), 'o', 'MarkerSize', 9, 'LineWidth', 1.6, ...
                 'Color', [0.85 0.20 0.15]);
            text(ax, res(m,1), res(m,2), sprintf('  %d', m), 'Color', [0.85 0.20 0.15], ...
                 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
        end

        % The repolarization guide: only meaningful once baseline and peak exist.
        % It converts an unreliable judgement ("where is 80% repolarized?") into
        % a reliable one ("where does the trace cross this line?").
        if strcmp(o.Mode, 'apd') && size(res, 1) >= 2
            vb = trace_at(t_ms, y, res(1,1), o.BaselineWin);
            vp = trace_at(t_ms, y, res(2,1), 0);
            lv = vp - o.Percent / 100 * (vp - vb);
            plot(ax, [t_ms(1) t_ms(end)], [lv lv], '--', 'Color', [0.10 0.60 0.30], 'LineWidth', 1.2);
            text(ax, t_ms(1), lv, sprintf(' %g%% level ', o.Percent), 'Color', [0.10 0.60 0.30], ...
                 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
        end

        if j <= nclick
            title(ax, {sprintf('Pixel %d of %d   —   click %d of %d', k, npix, j, nclick), ...
                       click_label(o.Mode, j, o.Percent), ...
                       'u undo   r restart   s skip   b back   q save+quit'}, ...
                      'FontSize', 11);
        else
            title(ax, {sprintf('Pixel %d of %d   —   all %d clicks placed', k, npix, nclick), ...
                       'ENTER or click to accept', ...
                       'u undo   r restart   s skip   b back   q save+quit'}, ...
                      'FontSize', 11);
        end
        hold(ax, 'off');

        % ---- one interaction ----
        try
            was_key = waitforbuttonpress;
        catch
            action = 'quit'; return;               % figure closed mid-wait
        end
        if ~ishandle(fig); action = 'quit'; return; end

        if was_key
            ch = lower(get(fig, 'CurrentCharacter'));
            if isempty(ch); continue; end
            switch ch
                case 'u'
                    if ~isempty(res); res(end, :) = []; end
                case 'r'
                    res = zeros(0, 2);
                case 's'
                    action = 'skip'; return;
                case 'b'
                    action = 'back'; return;
                case 'q'
                    action = 'quit'; return;
                case {char(13), char(3)}                    % Enter
                    if size(res, 1) == nclick
                        action = 'done';
                        res(:, 3) = toc(t0);                % stash elapsed
                        return;
                    end
            end
        else
            cp = get(ax, 'CurrentPoint');
            x  = cp(1, 1);
            if x < t_ms(1) || x > t_ms(end)
                continue;                                    % click outside the axes
            end
            if size(res, 1) < nclick
                res(end+1, :) = [x, cp(1, 2)];               %#ok<AGROW>
                if size(res, 1) == nclick
                    % Redraw once so the reviewer sees the final placement,
                    % then require Enter — prevents an accidental 5th click
                    % from silently committing a bad mark.
                    continue;
                end
            else
                action = 'done';
                res(:, 3) = toc(t0);
                return;
            end
        end
    end
end


function marks = store(marks, i, res, o, t_ms, y)
%STORE  Turn raw clicks into recorded fiducials.  No searching, no detection.

    marks.clicks{i}  = res(:, 1:2);
    marks.seconds(i) = res(1, 3);

    % Baseline: median over a short window about the clicked time.  The only
    % smoothing anywhere in this harness, and it is applied to diastole only.
    marks.t_base_ms(i) = res(1, 1);
    marks.v_base(i)    = trace_at(t_ms, y, res(1, 1), o.BaselineWin);

    % Peak(s): value read AT the clicked time — no local-maximum search.
    marks.t_peak_ms(i) = res(2, 1);
    marks.v_peak(i)    = trace_at(t_ms, y, res(2, 1), 0);

    if strcmp(o.Mode, 'apd')
        marks.t_act_ms(i) = res(3, 1);      % the measurement is the clicked time
        marks.t_rep_ms(i) = res(4, 1);
    else
        marks.t_peak2_ms(i) = res(3, 1);
        marks.v_peak2(i)    = trace_at(t_ms, y, res(3, 1), 0);
    end
end


function v = trace_at(t_ms, y, t_click, half_win_ms)
%TRACE_AT  Trace value at a clicked time.
%   half_win_ms == 0 -> the sample nearest the click.
%   half_win_ms  > 0 -> median over +/- half_win_ms about the click.

    if half_win_ms <= 0
        [~, ix] = min(abs(t_ms - t_click));
        v = y(ix);
    else
        sel = abs(t_ms - t_click) <= half_win_ms;
        if ~any(sel)
            [~, ix] = min(abs(t_ms - t_click));
            sel = false(size(t_ms)); sel(ix) = true;
        end
        v = median(y(sel), 'omitnan');
    end
end


function save_marks(f, marks)
    marks.updated = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    save(f, 'marks');
end
