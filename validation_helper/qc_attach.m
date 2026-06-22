function qc = qc_attach(qc, new_recs)
%QC_ATTACH  Merge QC records into cmos_all_data.qc, idempotently per stage.
%
%   cmos_all_data.qc = qc_attach(cmos_all_data.qc, new_recs)
%
%   Appends new_recs to the per-recording QC log carried inside cmos_all_data.
%   Records for any stage present in new_recs are REPLACED, so re-running a
%   stage (e.g. re-extracting features with a new min_cycles) refreshes that
%   stage's QC instead of duplicating it.  Records from other stages are kept.
%
%   Inputs
%     qc         existing QC log (struct array of qc_check records), or [].
%     new_recs   qc_check records for the stage(s) just produced.
%
%   Output
%     qc         merged log.

    if isempty(qc)
        qc = repmat(qc_check("","","" ,"PASS",NaN,""), 0, 1);
    end
    if isempty(new_recs)
        return;
    end

    new_recs = new_recs(:);
    stages   = unique([new_recs.stage]);

    % Drop existing records belonging to the stage(s) being refreshed.
    if ~isempty(qc)
        keep = arrayfun(@(r) ~ismember(r.stage, stages), qc);
        qc   = qc(keep);
    end

    qc = [qc; new_recs];
end
