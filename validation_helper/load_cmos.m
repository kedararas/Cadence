function [cmos_all_data, recs] = load_cmos(file, stage, on_fail)
%LOAD_CMOS  Load a cmos_all_data .mat and validate its schema at the boundary.
%
%   cmos_all_data            = load_cmos(file, stage)
%   [cmos_all_data, recs]    = load_cmos(file, stage, on_fail)
%
%   Use this EVERYWHERE you currently call load() on a pipeline .mat, so the
%   structural contract is checked before any stage touches the data.
%
%   stage     "raw" | "conditioned" | "metrics" — which contract to enforce.
%   on_fail   "error" (default) | "warn".
%               "error" — throw on a schema FAIL (fail fast; the data is not
%                         usable by the next stage anyway).
%               "warn"  — warn and return anyway (useful in a batch loop that
%                         wants to record the failure and move on).
%
%   recs      schema qc_check records (attach to cmos_all_data.qc / the batch
%             aggregate just like the other validators).
%
%   The .mat is expected to contain a variable 'cmos_all_data'; if not, the
%   first struct variable in the file is used.

    if nargin < 3 || isempty(on_fail), on_fail = "error"; end
    on_fail = lower(string(on_fail));

    S = load(file);
    if isfield(S, 'cmos_all_data')
        cmos_all_data = S.cmos_all_data;
    else
        fn = fieldnames(S);
        is_struct = structfun(@isstruct, S);
        if ~any(is_struct)
            error('load_cmos:noStruct', '%s contains no cmos_all_data struct.', file);
        end
        cmos_all_data = S.(fn{find(is_struct, 1)});
    end

    [~, base, ext] = fileparts(file);
    recs = validate_schema(cmos_all_data, stage, [base ext]);

    failed = recs(arrayfun(@(r) r.status == "FAIL", recs));
    if ~isempty(failed)
        msg = sprintf('Schema FAIL on %s (stage=%s):', [base ext], stage);
        for k = 1:numel(failed)
            msg = sprintf('%s\n  - %s', msg, failed(k).message);
        end
        if on_fail == "error"
            error('load_cmos:schema', '%s', msg);
        else
            warning('load_cmos:schema', '%s', msg);
        end
    end
end
