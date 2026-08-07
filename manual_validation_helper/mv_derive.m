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

    p = inputParser;
    p.addParameter('Output', '', @(v) ischar(v) || isstring(v));
    p.parse(varargin{:});
    o = p.Results;

    if ischar(marks_files) || isstring(marks_files)
        marks_files = {char(marks_files)};
    end

    tbl = table();
    mode_seen = '';

    for f = 1:numel(marks_files)
        R = load(char(marks_files{f}));
        if ~isfield(R, 'marks')
            error('mv_derive:badFile', '%s contains no ''marks'' struct.', marks_files{f});
        end
        m = R.marks;

        if isempty(mode_seen)
            mode_seen = m.mode;
        elseif ~strcmp(mode_seen, m.mode)
            error('mv_derive:mixedModes', ...
                  'Cannot combine marks from different modes (''%s'' and ''%s'').', ...
                  mode_seen, m.mode);
        end

        M = load(m.manifest_file);
        manifest = M.manifest;

        keep = find(m.done);
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

                % A negative APD means repolarization was marked before
                % activation — a slip, not a measurement.  Drop it loudly
                % rather than letting it widen the limits of agreement.
                bad = t.apd_ms <= 0;
                if any(bad)
                    warning('mv_derive:nonPositiveAPD', ...
                            '%s: dropping %d pixel(s) with APD <= 0 (repolarization marked before activation).', ...
                            m.reviewer, sum(bad));
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

        ud = struct('mode', m.mode, 'percent', m.percent, 'cam', m.cam, ...
                    'manifest_file', m.manifest_file, 'metrics_hint', '');
        tbl = [tbl; t]; %#ok<AGROW>
    end

    if isempty(tbl)
        error('mv_derive:empty', 'No completed marks found in any input file.');
    end

    % Carried through for mv_compare: it needs the camera to size the maps and
    % the manifest to recover the recording geometry.
    ud.nreviewers = numel(unique(tbl.reviewer));
    tbl.Properties.UserData = ud;

    fprintf('Derived %d rows from %d reviewer file(s), mode ''%s''.\n', ...
            height(tbl), numel(marks_files), mode_seen);
    fprintf('  reviewers : %s\n', strjoin(cellstr(unique(tbl.reviewer)), ', '));
    fprintf('  pixels    : %d unique\n', numel(unique(tbl.pixel)));
    fprintf('  median time per pixel: %.1f s\n', median(tbl.seconds, 'omitnan'));

    if ~isempty(o.Output)
        save(char(o.Output), 'tbl');
        fprintf('  saved     : %s\n', o.Output);
    end
end
