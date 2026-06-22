function [summary, detail] = qc_summarize(recs, out_csv)
%QC_SUMMARIZE  Roll a batch of qc_check records into a report.
%
%   [summary, detail] = qc_summarize(recs)
%   [summary, detail] = qc_summarize(recs, out_csv)
%
%   recs      array of qc_check records accumulated across a whole batch
%             (all recordings, all stages).
%   out_csv   optional path; if given, the full per-check detail table is
%             written there.
%
%   Outputs
%     summary  table: one row per recording, worst status per stage, and an
%              overall verdict — the at-a-glance "which recordings need a look".
%     detail   table: every individual check (the drill-down).
%
%   Also prints a console summary: counts per status, and the list of
%   recordings that FAILed or were flagged WARN, with their reasons.

    if isempty(recs)
        warning('qc_summarize:empty', 'No QC records to summarize.');
        summary = table(); detail = table();
        return;
    end

    detail = struct2table(recs);

    % Severity rank for "worst status" rollups.
    rank = @(s) double(s == "WARN") + 2*double(s == "FAIL");   % PASS=0,WARN=1,FAIL=2
    unrank = ["PASS","WARN","FAIL"];

    recordings = unique(detail.recording, 'stable');
    stages     = unique(detail.stage, 'stable');

    nR = numel(recordings);
    summary = table('Size', [nR, 2 + numel(stages)], ...
        'VariableTypes', ["string", repmat("string",1,numel(stages)), "string"], ...
        'VariableNames', ["recording", stages(:)', "overall"]);

    for r = 1:nR
        rec  = recordings(r);
        mask = detail.recording == rec;
        summary.recording(r) = rec;
        worst_overall = 0;
        for s = 1:numel(stages)
            sm = mask & (detail.stage == stages(s));
            if ~any(sm)
                summary.(stages(s))(r) = "-";
                continue;
            end
            w = max(arrayfun(rank, detail.status(sm)));
            summary.(stages(s))(r) = unrank(w+1);
            worst_overall = max(worst_overall, w);
        end
        summary.overall(r) = unrank(worst_overall+1);
    end

    % ---- console report ----
    nPASS = sum(detail.status == "PASS");
    nWARN = sum(detail.status == "WARN");
    nFAIL = sum(detail.status == "FAIL");
    fprintf('\n==================== CADENCE batch QC ====================\n');
    fprintf('Recordings: %d   Checks: %d   (PASS %d | WARN %d | FAIL %d)\n', ...
        nR, height(detail), nPASS, nWARN, nFAIL);

    failed  = summary.recording(summary.overall == "FAIL");
    flagged = summary.recording(summary.overall == "WARN");
    fprintf('Overall:  %d FAIL | %d WARN | %d PASS\n', ...
        numel(failed), numel(flagged), nR - numel(failed) - numel(flagged));

    if ~isempty(failed) || ~isempty(flagged)
        fprintf('\n-- recordings needing review --\n');
        bad = detail(detail.status ~= "PASS", :);
        for r = [failed(:)', flagged(:)']
            rows = bad(bad.recording == r, :);
            for k = 1:height(rows)
                fprintf('  [%-4s] %-28s %-16s %s\n', rows.status(k), r, ...
                    rows.stage(k) + "/" + rows.check(k), rows.message(k));
            end
        end
    end
    fprintf('==========================================================\n\n');

    % ---- optional CSV ----
    if nargin >= 2 && ~isempty(out_csv)
        writetable(detail, out_csv);
        fprintf('QC detail written to %s\n', out_csv);
    end
end
