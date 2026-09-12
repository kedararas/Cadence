function tbl = mv_derive(marks_files, varargin)
%MV_DERIVE  Turn saved fiducial clicks into manual metric values.
%
%   tbl = mv_derive(marks_file)
%   tbl = mv_derive({marks_A, marks_B, marks_C})
%   tbl = mv_derive(..., 'Output', 'manual_values.mat')
%
%   Pure arithmetic on what the reviewers marked.  No detection, no search, no
%   thresholds — every judgement in the result was made by a human in mv_mark.
%   Keeping this separate from marking is deliberate: the raw clicks are the
%   record, so a metric definition can be revised (APD80 -> APD50, a different
%   alternans normalisation) without asking anyone to mark anything again.
%
%   Returns one row per (reviewer, pixel), with the pixel's SNR and stratum
%   carried through from the manifest so mv_compare can stratify.
%
%   Columns
%     reviewer, pixel, row, col, snr, stratum, seconds
%     'apd' mode        : act_ms, apd_ms, amplitude
%     'alternans' mode  : amp1, amp2, alt_ratio
%
%   Note on row/col: they are attached HERE, after marking is finished.  The
%   marking tool withholds them on purpose (see mv_mark).
%
%   Slips are dropped, never corrected: in 'apd' mode a mark whose clicks are
%   out of order (activation after the peak, repolarization before it) or
%   whose APD is <= 0.  Both rules read the clicks only.  The dropped marks
%   are listed in tbl.Properties.UserData.dropped (reviewer, pixel, reason).

    p = inputParser;
    p.addParameter('Output', '', @(v) ischar(v) || isstring(v));
    p.parse(varargin{:});
    o = p.Results;

    if ischar(marks_files) || isstring(marks_files)
        marks_files = {char(marks_files)};
    end

    tbl = table();
    dropped = drop_rows('', [], "");       % marks removed as slips, for reporting
    mode_seen = '';

    for f = 1:numel(marks_files)
        R = load(char(marks_files{f}));
        if ~isfield(R, 'marks')
            error('mv_derive:badFile', '%s contains no ''marks'' struct.', marks_files{f});
        end
        m = R.marks;

        if isempty(mode_seen)
            mode_seen = m.mode;
            manifest_seen = m.manifest_file;
        elseif ~strcmp(mode_seen, m.mode)
            error('mv_derive:mixedModes', ...
                  'Cannot combine marks from different modes (''%s'' and ''%s'').', ...
                  mode_seen, m.mode);
        elseif ~strcmp(manifest_seen, m.manifest_file)
            % The pixel column is an index INTO THE MANIFEST.  Pooling marks
            % from two manifests makes pixel k mean different pixels for
            % different reviewers, so mv_compare's intersect() would pair
            % unrelated locations — inflating the inter-observer spread and
            % making the software look better by comparison.
            error('mv_derive:mixedManifests', ...
                  ['Marks files reference different manifests:\n  %s\n  %s\n' ...
                   'Reviewers must mark the SAME pixel list for the inter-observer ' ...
                   'comparison to mean anything.'], manifest_seen, m.manifest_file);
        end

        M = load(m.manifest_file);
        manifest = M.manifest;

        % Only the finalised POOL enters the study.  A lead reviewer's marks
        % file may also hold accepted candidates beyond the target (kept out by
        % mv_finalize_pool) — those are dropped here too.  Manifests from
        % before candidate pools existed have every pixel active.
        if isfield(manifest, 'pool_status')
            if strcmp(manifest.pool_status, 'open')
                error('mv_derive:poolOpen', ...
                      ['The pool in %s is still OPEN. Finish the lead session (mv_mark with ' ...
                       '''Lead'', true) or run mv_finalize_pool before deriving values.'], ...
                      m.manifest_file);
            end
            active = logical(manifest.active(:));
        else
            active = true(numel(m.done), 1);
        end

        keep = find(m.done(:) & active);
        if isempty(keep)
            warning('mv_derive:noMarks', '%s has no completed pixels.', marks_files{f});
            continue;
        end

        n = numel(keep);
        t = table();
        t.reviewer = repmat(string(m.reviewer), n, 1);
        t.pixel    = keep;
        t.row      = manifest.pixels(keep, 1);
        t.col      = manifest.pixels(keep, 2);
        t.snr      = manifest.snr(keep);
        t.stratum  = manifest.stratum(keep);
        t.seconds  = m.seconds(keep);

        switch m.mode
            case 'apd'
                t.act_ms    = m.t_act_ms(keep);
                t.apd_ms    = m.t_rep_ms(keep) - m.t_act_ms(keep);
                t.amplitude = m.v_peak(keep) - m.v_base(keep);

                % Clicks out of order are slips too.  Activation lies on the
                % upstroke, so it cannot follow the peak, and repolarization
                % cannot precede it.  The usual form is an activation click on
                % the DOWNSTROKE crossing of the 50% guide line: it gives a
                % short but positive APD, so the check below misses it.  The
                % rule reads only the clicks, never a CADENCE value, so it
                % cannot favour the software.  Baseline is not ordered: late
                % diastole is a legitimate place to click it.
                t_peak = m.t_peak_ms(keep);
                bad = t.act_ms > t_peak | m.t_rep_ms(keep) < t_peak;
                if any(bad)
                    warning('mv_derive:clickOrder', ...
                            ['%s: dropping %d pixel(s) with clicks out of order (activation after ' ...
                             'peak, or repolarization before peak): pixel %s.'], ...
                            m.reviewer, sum(bad), mat2str(t.pixel(bad)'));
                    dropped = [dropped; drop_rows(m.reviewer, t.pixel(bad), "click order")]; %#ok<AGROW>
                    t(bad, :) = [];
                end

                % A negative APD means repolarization was marked before
                % activation — a slip, not a measurement.  Drop it loudly
                % rather than letting it widen the limits of agreement.
                bad = t.apd_ms <= 0;
                if any(bad)
                    warning('mv_derive:nonPositiveAPD', ...
                            '%s: dropping %d pixel(s) with APD <= 0 (repolarization marked before activation).', ...
                            m.reviewer, sum(bad));
                    dropped = [dropped; drop_rows(m.reviewer, t.pixel(bad), "APD <= 0")]; %#ok<AGROW>
                    t(bad, :) = [];
                end

            case 'alternans'
                a1 = m.v_peak(keep)  - m.v_base(keep);
                a2 = m.v_peak2(keep) - m.v_base(keep);
                t.amp1      = a1;
                t.amp2      = a2;
                % |A2 - A1| / mean(A1, A2): dimensionless, sign-free, and the
                % same normalisation the alternans maps use, so the two are
                % directly comparable.
                t.alt_ratio = abs(a2 - a1) ./ ((a1 + a2) / 2);

                bad = ~isfinite(t.alt_ratio) | (a1 + a2) <= 0;
                if any(bad)
                    warning('mv_derive:badAlternans', ...
                            '%s: dropping %d pixel(s) with non-positive or non-finite amplitudes.', ...
                            m.reviewer, sum(bad));
                    t(bad, :) = [];
                end
        end

        % act_percent records which upstroke crossing the reviewers were shown
        % as the activation guide, so a table can be traced back to the
        % definition it was marked against.  Absent in pre-guide marks files.
        if isfield(m, 'act_percent'); act_pct = m.act_percent; else; act_pct = NaN; end
        ud = struct('mode', m.mode, 'percent', m.percent, 'act_percent', act_pct, ...
                    'cam', m.cam, 'manifest_file', m.manifest_file, 'metrics_hint', '');
        tbl = [tbl; t]; %#ok<AGROW>
    end

    if isempty(tbl)
        error('mv_derive:empty', 'No completed marks found in any input file.');
    end

    % Carried through for mv_compare: it needs the camera to size the maps and
    % the manifest to recover the recording geometry.
    ud.nreviewers = numel(unique(tbl.reviewer));
    ud.dropped    = dropped;
    tbl.Properties.UserData = ud;

    fprintf('Derived %d rows from %d reviewer file(s), mode ''%s''.\n', ...
            height(tbl), numel(marks_files), mode_seen);
    if height(dropped) > 0
        fprintf('  dropped   : %d mark(s) as slips (listed in UserData.dropped)\n', height(dropped));
    end
    fprintf('  reviewers : %s\n', strjoin(cellstr(unique(tbl.reviewer)), ', '));
    fprintf('  pixels    : %d unique\n', numel(unique(tbl.pixel)));
    fprintf('  median time per pixel: %.1f s\n', median(tbl.seconds, 'omitnan'));

    if ~isempty(o.Output)
        save(char(o.Output), 'tbl');
        fprintf('  saved     : %s\n', o.Output);
    end
end


function T = drop_rows(reviewer, pixels, reason)
    n = numel(pixels);
    T = table(repmat(string(reviewer), n, 1), pixels(:), repmat(string(reason), n, 1), ...
              'VariableNames', {'reviewer', 'pixel', 'reason'});
end
