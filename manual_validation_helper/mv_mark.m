function marks_file = mv_mark(manifest_file, reviewer_id, varargin)
%MV_MARK  Blinded manual fiducial marking for software validation.
%
%   marks_file = mv_mark(manifest_file, reviewer_id)
%   marks_file = mv_mark(manifest_file, reviewer_id, 'Lead', true)
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
%                      Yields activation time, APD/CaTD at 'Percent', amplitude.
%     'alternans'      Three clicks: baseline, peak of beat 1, peak of beat 2.
%                      Yields the two amplitudes and their ratio.  Build the
%                      manifest with 'NumBeats', 2.
%
%   GUIDE LINES ('apd' mode)
%   Once baseline and peak are placed, the level for the click being asked for
%   is drawn across the axes, so the reviewer answers "where does the trace
%   cross this line?" instead of "where is 80% repolarized?".  The first is a
%   judgement humans make reliably; the second is not.
%
%   This applies to ACTIVATION as much as to repolarization, and that matters:
%   CADENCE times activation at the 50% baseline-to-peak crossing of the
%   upstroke (compute_lat_50).  Asking a reviewer for "the upstroke" invites
%   them to mark maximum dV/dt, which is a different definition — the resulting
%   offset would be booked as software-versus-human disagreement when it is
%   really the two sides measuring different things.  'ActPercent' therefore
%   defaults to 50 to match the code under test; change it only if the pipeline
%   definition changes.
%
%   HOW A CLICK IS INTERPRETED
%   The reviewer's click fixes a TIME.  The harness then reads the trace value
%   at that time — it does not search for a nearby maximum, which would be
%   exactly the automated detection under test.  The one exception is the
%   baseline, taken as the median of the trace over +/-'BaselineWin' ms around
%   the clicked time, because a single diastolic sample is noise-dominated.
%   Both behaviours are recorded in the output so they can be reported.
%
%   CLICK ORDER ('apd' mode)
%   An activation click after the reviewer's own peak click, or a
%   repolarization click before it, is refused on the spot and the reviewer is
%   asked to click again.  The 50% guide line is crossed twice, and the
%   downstroke crossing is an easy slip that yields a plausible short APD.  The
%   check compares the reviewer's clicks with each other only; it never judges
%   where on the upstroke they clicked.  mv_derive drops marks that break the
%   same rule, for files written before this check existed.
%
%   LEAD SESSION AND THE POOL
%   mv_sample_pixels oversamples candidates and leaves the pool OPEN.  The
%   lead reviewer runs with 'Lead', true: candidates are presented in manifest
%   order (not shuffled), and the session stops by itself once the target
%   number of pixels has been ACCEPTED, at which point the pool is finalised
%   (mv_finalize_pool) — accepted pixels in, skipped ones recorded, unreached
%   candidates discarded.  Everyone else marks only the finalised pool, in a
%   per-reviewer random order, and is refused while the pool is still open.
%   The lead is blinded exactly like everyone else.
%
%   KEYS (while marking)
%     u  undo last click      r  restart this pixel     s  skip (unmarkable)
%     b  back one pixel       q  save and quit
%
%   A progress line and bar under the trace show marked / skipped / remaining
%   (for the lead: accepted against the target, plus the candidate count).
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
%     'ActPercent'   upstroke percentage defining activation (default 50, to
%                    match compute_lat_50)
%     'BaselineWin'  half-width in ms for the baseline median (default 5)
%     'Lead'         true for the lead reviewer building the pool (default false)
%     'Preview'      draw the first pixel's screen and return without waiting
%                    for input or writing anything (default false; for checking
%                    the display, e.g. from a script)
%     'Output'       marks path (default alongside the manifest)

    p = inputParser;
    p.addParameter('Lead',        false, @(v) islogical(v) || isnumeric(v));
    p.addParameter('Preview',     false, @(v) islogical(v) || isnumeric(v));
    p.addParameter('Mode',        'apd', @(v) any(strcmpi(v, {'apd','alternans'})));
    p.addParameter('Percent',     80,    @(v) isnumeric(v) && isscalar(v) && v > 0 && v < 100);
    p.addParameter('ActPercent',  50,    @(v) isnumeric(v) && isscalar(v) && v > 0 && v < 100);
    p.addParameter('BaselineWin', 5,     @(v) isnumeric(v) && isscalar(v) && v > 0);
    p.addParameter('Output',      '',    @(v) ischar(v) || isstring(v));
    p.parse(varargin{:});
    o = p.Results;
    o.Mode = lower(o.Mode);
    o.Lead = logical(o.Lead);

    reviewer = char(string(reviewer_id));
    if isempty(reviewer)
        error('mv_mark:reviewer', 'reviewer_id must be a non-empty string.');
    end

    M = load(manifest_file);
    if ~isfield(M, 'manifest')
        error('mv_mark:badManifest', '%s does not contain a ''manifest'' struct.', manifest_file);
    end
    manifest = M.manifest;

    % ---- pool status: who may mark what ----
    % Manifests written before candidate pools existed are treated as final
    % with every pixel active.
    n_all = size(manifest.pixels, 1);
    if isfield(manifest, 'pool_status')
        pool_open = strcmp(manifest.pool_status, 'open');
        active    = logical(manifest.active(:));
        n_target  = manifest.n_target;
    else
        pool_open = false;
        active    = true(n_all, 1);
        n_target  = n_all;
    end
    if pool_open && ~o.Lead
        error('mv_mark:poolOpen', ...
              ['The pool in %s is still OPEN: the lead reviewer has not finished building ' ...
               'it. Wait for the lead session to complete (or, if you are the lead, pass ' ...
               '''Lead'', true).'], manifest_file);
    end
    if ~pool_open && o.Lead
        warning('mv_mark:poolFinal', ...
                'Pool already final (%d pixels); continuing as an ordinary reviewer session.', ...
                sum(active));
        o.Lead = false;
    end
    if o.Lead
        shown = (1:n_all)';                 % manifest order: tiers / strata interleaved
        % Stratified draws are balanced BY STRATUM, not just in total: the
        % lead fills a per-stratum quota, and once a stratum is full its
        % remaining candidates are passed over (a dim stratum is skipped more
        % often, and a global count would let the bright ones crowd it out).
        % Grid pools fill globally, tier order doing the spreading.
        if isfield(manifest, 'sampling') && strcmp(manifest.sampling, 'stratified')
            quota = floor(n_target / manifest.num_strata);
            stratum_of = manifest.stratum(:);
        else
            quota = Inf;
            stratum_of = ones(n_all, 1);
        end
    else
        shown = find(active);
        quota = Inf; stratum_of = ones(n_all, 1);
    end

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
    % The manifest names the array the reviewers see — CAM<n>_average under the
    % recommended ensemble workflow, CAM<n> for the raw single-beat arm.  Older
    % manifests predate the field and were always raw.
    if isfield(manifest, 'data_field') && ~isempty(manifest.data_field)
        data_field = char(manifest.data_field);
    else
        data_field = char(manifest.cam);
    end
    if ~isfield(d, data_field)
        error('mv_mark:noField', '%s has no field %s.', manifest.conditioned_file, data_field);
    end
    X  = d.(data_field);
    nt = size(X, 3);
    wf = manifest.window_frames;
    T  = double(reshape(X, [], nt));
    T  = T(manifest.linear_index, wf(1):wf(2));            % nPix x nFrames
    t_ms = ((wf(1):wf(2)) - 1) / manifest.acqFreq * 1000;  % absolute ms
    npix = size(T, 1);                                     % ALL candidates

    % ---- resume or start ----
    if isempty(o.Output)
        [pdir, base] = fileparts(manifest_file);
        o.Output = fullfile(pdir, sprintf('%s_%s_marks.mat', base, reviewer));
    end
    o.Output = char(o.Output);

    if isfile(o.Output)
        R = load(o.Output);
        marks = R.marks;
        % Marks files written before ActPercent existed carry the 50% behaviour
        % this default reproduces, so treat a missing field as 50 rather than
        % refusing to resume them.
        if isfield(marks, 'act_percent'); prev_act = marks.act_percent; else; prev_act = 50; end
        if ~strcmp(marks.mode, o.Mode) || marks.percent ~= o.Percent || ...
           (strcmp(o.Mode, 'apd') && prev_act ~= o.ActPercent)
            error('mv_mark:resumeMismatch', ...
                  ['%s was started with Mode=''%s'', Percent=%g, ActPercent=%g but you ' ...
                   'passed Mode=''%s'', Percent=%g, ActPercent=%g. Use the original ' ...
                   'settings or a new Output file.'], ...
                  o.Output, marks.mode, marks.percent, prev_act, o.Mode, o.Percent, o.ActPercent);
        end
        if isfield(marks, 'lead') && marks.lead ~= o.Lead
            error('mv_mark:resumeLead', ...
                  '%s was started as a %s session; pass the same ''Lead'' setting.', ...
                  o.Output, ternary(marks.lead, 'lead', 'non-lead'));
        end
        fprintf('Resuming %s: %s\n', o.Output, progress_line(marks, o.Lead, n_target, numel(shown)));
    else
        marks = new_marks(manifest, manifest_file, reviewer, o, npix, shown);
        if o.Lead
            fprintf('New LEAD session for reviewer ''%s'': %d candidates, target %d.\n', ...
                    reviewer, numel(shown), n_target);
        else
            fprintf('New marking session for reviewer ''%s'': %d pixels.\n', reviewer, numel(shown));
        end
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
    nshow = numel(marks.order);
    visited = marks.done(marks.order) | marks.skipped(marks.order);
    k = find(~visited, 1);
    if isempty(k); k = nshow + 1; end
    finalise = false;

    while k >= 1 && k <= nshow
        if ~ishandle(fig)
            fprintf('Figure closed — saving and exiting.\n'); break;
        end
        i = marks.order(k);                     % index into manifest.pixels
        if o.Lead && isfinite(quota) && ~marks.done(i) && ~marks.skipped(i) && ...
           sum(marks.done(stratum_of == stratum_of(i))) >= quota
            k = k + 1;                          % this stratum is full: pass over
            continue;
        end
        y = T(i, :);

        prog = progress_info(marks, o.Lead, n_target, nshow, k);
        [res, action] = mark_one(fig, t_ms, y, prog, o, nclick);
        if o.Preview
            marks_file = '';                    % nothing written, figure left open
            return;
        end

        switch action
            case 'done'
                marks = store(marks, i, res, o, t_ms, y);
                marks.done(i) = true; marks.skipped(i) = false;
                k = k + 1;
                if o.Lead
                    if isfinite(quota)
                        got = arrayfun(@(s) sum(marks.done(stratum_of == s)), 1:manifest.num_strata);
                        finalise = all(got >= quota);
                    else
                        finalise = sum(marks.done) >= n_target;
                    end
                end
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
        if finalise; break; end
    end

    if ishandle(fig); close(fig); end
    save_marks(o.Output, marks);
    marks_file = o.Output;

    if o.Lead
        if ~finalise
            % Candidates exhausted before the target: finalise with what there is,
            % loudly (mv_finalize_pool warns and records the shortfall).
            fprintf('Candidate list exhausted at %d accepted (target %d).\n', ...
                    sum(marks.done), n_target);
        end
        mv_finalize_pool(manifest_file, 'Marks', o.Output);
    else
        fprintf('Complete: %d marked, %d skipped -> %s\n', ...
                sum(marks.done), sum(marks.skipped), o.Output);
    end
    fprintf('Next: mv_derive(''%s'')\n', o.Output);
end


% ======================================================================
function prog = progress_info(marks, lead, n_target, nshow, k)
%PROGRESS_INFO  What the reviewer has done and how much is left.
    prog = struct();
    prog.done    = sum(marks.done);
    prog.skipped = sum(marks.skipped);
    prog.lead    = lead;
    prog.k       = k;
    prog.nshow   = nshow;
    if lead
        prog.target    = n_target;
        prog.remaining = max(0, n_target - prog.done);
    else
        prog.target    = nshow;
        prog.remaining = max(0, nshow - prog.done - prog.skipped);
    end
    prog.frac = min(1, prog.done / max(prog.target, 1));
    prog.line = progress_line(marks, lead, n_target, nshow, k);
end


function s = progress_line(marks, lead, n_target, nshow, k)
    nd = sum(marks.done); ns = sum(marks.skipped);
    if lead
        s = sprintf('accepted %d of %d target  |  skipped %d  |  %d to go', ...
                    nd, n_target, ns, max(0, n_target - nd));
        if nargin >= 5
            s = sprintf('%s  (candidate %d of %d)', s, k, nshow);
        end
    else
        s = sprintf('marked %d  |  skipped %d  |  %d of %d remaining', ...
                    nd, ns, max(0, nshow - nd - ns), nshow);
        if nargin >= 5
            s = sprintf('%s  (pixel %d of %d)', s, k, nshow);
        end
    end
end


function out = ternary(c, a, b)
    if c; out = a; else; out = b; end
end


% ======================================================================
function marks = new_marks(manifest, manifest_file, reviewer, o, npix, shown)
    % Presentation order is randomized PER REVIEWER: it removes order/fatigue
    % effects from the inter-observer comparison and makes it impractical for
    % two reviewers to sit together and mark "the same one" in step.
    % The LEAD is the exception: candidates go in manifest order, because that
    % order is what keeps strata balanced (stratified) or the lattice evenly
    % filled (grid) wherever the session stops.
    seed = mod(sum(double(reviewer)) * 7919 + manifest.seed, 2^31 - 1);
    rng(seed, 'twister');
    shown = shown(:);
    if o.Lead
        order = shown;
    else
        order = shown(randperm(numel(shown)));
    end

    marks = struct();
    marks.manifest_file   = char(manifest_file);
    marks.conditioned_file = manifest.conditioned_file;
    marks.cam             = manifest.cam;
    marks.reviewer        = reviewer;
    marks.mode            = o.Mode;
    marks.percent         = o.Percent;
    marks.act_percent     = o.ActPercent;
    marks.baseline_win_ms = o.BaselineWin;
    marks.lead            = o.Lead;
    marks.order_seed      = seed;
    marks.order           = order;             % indices into manifest.pixels
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


function lbl = click_label(mode, j, pct, act_pct)
    if strcmp(mode, 'apd')
        switch j
            case 1, lbl = 'DIASTOLIC BASELINE (click in diastole, before the upstroke)';
            case 2, lbl = 'PEAK (click at the time of the peak)';
            case 3, lbl = sprintf('ACTIVATION — click where the UPSTROKE crosses the %g%% guide line', act_pct);
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


function why = order_problem(mode, res, x)
%ORDER_PROBLEM  Message if a click at time x would be out of order, else ''.
%   Same rule mv_derive applies afterwards: activation before the peak click,
%   repolarization after it.  Compares clicks with clicks, nothing else.
    why = '';
    if ~strcmp(mode, 'apd') || size(res, 1) < 2, return; end
    t_peak = res(2, 1);
    switch size(res, 1) + 1
        case 3
            if x > t_peak
                why = 'Refused: activation must be on the UPSTROKE, before your peak click. Click again.';
            end
        case 4
            if x < t_peak
                why = 'Refused: repolarization must come after your peak click. Click again.';
            end
    end
end


function [res, action] = mark_one(fig, t_ms, y, prog, o, nclick)
%MARK_ONE  Collect one pixel's clicks.  Returns raw click coordinates.

    res  = zeros(0, 2);
    t0   = tic;
    note = '';                                   % refusal message for the last click

    while true
        if ~ishandle(fig); action = 'quit'; return; end
        % The window can be closed while a redraw is in flight (the handle
        % test above passes, then the axes call finds a deleted figure).
        % Treat that like any other close: save and quit, never error out.
        try
            figure(fig); clf(fig);
            set(fig, 'Units', 'normalized');
            ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.10 0.20 0.85 0.64]);
            progress_bar(fig, prog);
        catch
            action = 'quit'; return;
        end
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

        % Guide lines, drawn as each becomes relevant and then left up so the
        % reviewer can check the final placement.  Both are only meaningful
        % once baseline and peak exist, since both levels are defined from
        % them.  See the header for why activation gets one too.
        if strcmp(o.Mode, 'apd') && size(res, 1) >= 2
            vb = trace_at(t_ms, y, res(1,1), o.BaselineWin);
            vp = trace_at(t_ms, y, res(2,1), 0);
            if j >= 3
                la = vb + o.ActPercent / 100 * (vp - vb);
                guide_line(ax, t_ms, la, sprintf(' %g%% upstroke ', o.ActPercent), [0.85 0.45 0.05]);
            end
            if j >= 4
                lv = vp - o.Percent / 100 * (vp - vb);
                guide_line(ax, t_ms, lv, sprintf(' %g%% repolarization ', o.Percent), [0.10 0.60 0.30]);
            end
        end

        if j <= nclick
            title(ax, {sprintf('Click %d of %d', j, nclick), ...
                       click_label(o.Mode, j, o.Percent, o.ActPercent), ...
                       'u undo   r restart   s skip   b back   q save+quit'}, ...
                      'FontSize', 11);
        else
            title(ax, {sprintf('All %d clicks placed', nclick), ...
                       'ENTER or click to accept', ...
                       'u undo   r restart   s skip   b back   q save+quit'}, ...
                      'FontSize', 11);
        end
        if ~isempty(note)
            text(ax, 0.5, 0.96, note, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'top', 'Color', [0.80 0.10 0.10], 'FontWeight', 'bold', ...
                 'BackgroundColor', 'w');
        end
        hold(ax, 'off');
        if o.Preview
            drawnow; action = 'preview'; return;
        end

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
                    note = '';
                case 'r'
                    res = zeros(0, 2);
                    note = '';
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
            % Only clicks INSIDE the trace axes count: a click on the progress
            % bar or the margins must not be read as a fiducial.
            fp = get(fig, 'CurrentPoint');                   % normalized figure units
            ap = get(ax, 'Position');
            if fp(1) < ap(1) || fp(1) > ap(1) + ap(3) || fp(2) < ap(2) || fp(2) > ap(2) + ap(4)
                continue;
            end
            cp = get(ax, 'CurrentPoint');
            x  = cp(1, 1);
            if x < t_ms(1) || x > t_ms(end)
                continue;                                    % click outside the axes
            end
            if size(res, 1) < nclick
                note = order_problem(o.Mode, res, x);
                if ~isempty(note)
                    continue;                                % refused: redraw with the message
                end
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


function progress_bar(fig, prog)
%PROGRESS_BAR  Marked / skipped / remaining, as a line and a bar under the trace.
    pax = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.10 0.05 0.85 0.025], ...
               'XLim', [0 1], 'YLim', [0 1], 'XTick', [], 'YTick', [], 'Box', 'on');
    hold(pax, 'on');
    patch(pax, [0 1 1 0], [0 0 1 1], [0.94 0.94 0.94], 'EdgeColor', 'none');
    if prog.frac > 0
        col = [0.20 0.55 0.30];
        if prog.lead; col = [0.20 0.35 0.75]; end
        patch(pax, [0 prog.frac prog.frac 0], [0 0 1 1], col, 'EdgeColor', 'none');
    end
    hold(pax, 'off');
    text(pax, 0, 1.6, prog.line, 'Units', 'data', 'FontSize', 10, ...
         'VerticalAlignment', 'bottom', 'Interpreter', 'none');
    text(pax, 1, 0.5, sprintf(' %.0f%%', 100 * prog.frac), 'Units', 'data', ...
         'HorizontalAlignment', 'left', 'FontSize', 9, 'Clipping', 'off');
    set(pax, 'HitTest', 'off');
end


function guide_line(ax, t_ms, level, label, col)
%GUIDE_LINE  Horizontal reference line the reviewer marks a crossing against.
    plot(ax, [t_ms(1) t_ms(end)], [level level], '--', 'Color', col, 'LineWidth', 1.2);
    text(ax, t_ms(1), level, label, 'Color', col, ...
         'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
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
