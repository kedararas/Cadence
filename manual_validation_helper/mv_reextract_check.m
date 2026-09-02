function report = mv_reextract_check(old_file, new_file)
%MV_REEXTRACT_CHECK  Compare a metrics file before and after re-extraction.
%
%   report = mv_reextract_check('old-metrics.mat', 'new-metrics.mat')
%
%   After the 2026-09-01 fixes, re-extraction should change SOME metrics and
%   leave others bit-identical.  That split is the check: extract_apd,
%   extract_act_time and extract_cv were not touched, so if their maps move,
%   something changed that was not intended and it needs explaining before the
%   figures are trusted.
%
%   EXPECTED UNCHANGED   apd_data, ca_data, act_times, local_cv
%   EXPECTED CHANGED     ap_rise_times, ca_rise_times   (interpolation + the
%                                                        amplitude reference)
%                        ca_tau                         (rewritten)
%                        vc_delay                       (signed; may now be
%                                                        absent if it was only
%                                                        ever produced in error)
%                        alternans_data                 (interpolated crossings)
%   CONDITIONAL          rep_data, ca_rep_data          (change only if
%                                                        Repolarization was run
%                                                        without APD, or both
%                                                        with different levels)
%
%   Compares slot 1 (the masked map) per camera.  Reports the max absolute
%   difference over pixels finite in both, and separately how many pixels
%   changed finite/NaN status — a metric can move a lot simply by censoring
%   more pixels, and that is worth seeing on its own.
%
%   Returns a table; prints a verdict line per field.

    E = {'apd_data','ca_data','act_times','local_cv'};                       % expect same
    C = {'ap_rise_times','ca_rise_times','ca_tau','vc_delay','alternans_data'}; % expect diff
    K = {'rep_data','ca_rep_data'};                                          % conditional

    a = load_metrics(old_file);
    b = load_metrics(new_file);

    fields = [E C K];
    rows = {};
    fprintf('\n%-16s %-10s %12s %10s %10s   %s\n', ...
            'field', 'expect', 'max|diff|', 'n both', 'NaN delta', 'verdict');
    fprintf('%s\n', repmat('-', 1, 82));

    for k = 1:numel(fields)
        f = fields{k};
        if any(strcmp(f, E))
            expect = 'same';
        elseif any(strcmp(f, C))
            expect = 'differ';
        else
            expect = 'maybe';
        end

        [ma, mb] = deal(get_maps(a, f), get_maps(b, f));
        if isempty(ma) && isempty(mb)
            fprintf('%-16s %-10s %12s %10s %10s   absent in both\n', f, expect, '-', '-', '-');
            continue;
        end
        if isempty(ma) || isempty(mb)
            fprintf('%-16s %-10s %12s %10s %10s   PRESENT IN ONLY ONE FILE\n', f, expect, '-', '-', '-');
            rows(end+1,:) = {f, expect, NaN, 0, NaN, 'one-sided'}; %#ok<AGROW>
            continue;
        end

        maxd = 0; nboth = 0; nan_delta = 0; shape_ok = true;
        for c = 1:min(numel(ma), numel(mb))
            x = ma{c}; y = mb{c};
            if isempty(x) || isempty(y); continue; end
            if ~isequal(size(x), size(y)); shape_ok = false; continue; end
            fx = isfinite(x); fy = isfinite(y);
            both = fx & fy;
            nboth = nboth + nnz(both);
            if any(both(:))
                maxd = max(maxd, max(abs(double(x(both)) - double(y(both)))));
            end
            nan_delta = nan_delta + nnz(fx ~= fy);
        end

        if ~shape_ok
            verdict = 'SIZE MISMATCH';
        elseif nboth == 0
            verdict = 'no overlapping pixels';
        elseif maxd == 0 && nan_delta == 0
            verdict = 'identical';
        else
            verdict = 'changed';
        end

        flag = '';
        if strcmp(expect, 'same') && ~strcmp(verdict, 'identical')
            flag = '   <== UNEXPECTED, investigate';
        elseif strcmp(expect, 'differ') && strcmp(verdict, 'identical')
            flag = '   <== expected a change, got none';
        end

        fprintf('%-16s %-10s %12.6g %10d %10d   %s%s\n', ...
                f, expect, maxd, nboth, nan_delta, verdict, flag);
        rows(end+1,:) = {f, expect, maxd, nboth, nan_delta, verdict}; %#ok<AGROW>
    end

    report = cell2table(rows, 'VariableNames', ...
        {'field','expected','max_abs_diff','n_pixels_both','nan_status_changed','verdict'});

    unexpected = strcmp(report.expected, 'same') & ~strcmp(report.verdict, 'identical');
    fprintf('\n');
    if any(unexpected)
        fprintf('%d field(s) expected to be identical are NOT. Do not trust the figures\n', nnz(unexpected));
        fprintf('until that is explained — extract_apd, extract_act_time and extract_cv\n');
        fprintf('were not modified, so a change there means something else moved.\n');
    else
        fprintf('All untouched metrics are bit-identical. The re-extraction changed only\n');
        fprintf('what it was supposed to.\n');
    end
end


% ======================================================================
function d = load_metrics(f)
    S = load(char(f));
    if isfield(S, 'cmos_all_data')
        d = S.cmos_all_data;
    else
        fn = fieldnames(S); is_s = structfun(@isstruct, S);
        if ~any(is_s)
            error('mv_reextract_check:noStruct', '%s contains no cmos_all_data struct.', f);
        end
        d = S.(fn{find(is_s, 1)});
    end
    if ~isfield(d, 'ep_metrics')
        error('mv_reextract_check:noMetrics', '%s has no ep_metrics.', f);
    end
end


function maps = get_maps(d, field)
%GET_MAPS  Slot-1 map per camera, as a cell.  {} if the field is absent.
    maps = {};
    if ~isfield(d.ep_metrics, field) || isempty(d.ep_metrics.(field)); return; end
    v = d.ep_metrics.(field);
    if isnumeric(v); maps = {v}; return; end
    if ~iscell(v); return; end
    for c = 1:size(v, 1)
        e = v{c, 1};
        if isnumeric(e) && ismatrix(e) && ~isempty(e)
            maps{end+1} = e; %#ok<AGROW>
        else
            maps{end+1} = []; %#ok<AGROW>
        end
    end
end
