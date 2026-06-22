function [ok, msgs, recs] = schema_gate(d, stage, recording)
%SCHEMA_GATE  Run a schema check at a module boundary and format the result.
%
%   [ok, msgs, recs] = schema_gate(d, stage, recording)
%
%   Thin wrapper around VALIDATE_SCHEMA for use *inside* the CADENCE App
%   Designer modules, where errors/warnings must be surfaced through the app's
%   own console rather than thrown.  Paste two lines at each load boundary:
%
%       [ok, msgs] = schema_gate(cmos_all_data, "conditioned", file_name);
%       for k = 1:numel(msgs), update_console(app, msgs{k}); end
%       if ~ok, return; end     % optional: refuse to work on a broken struct
%
%   Inputs
%     d          the struct just loaded (cmos_all_data, or metrics).
%     stage      "raw" | "conditioned" | "metrics".
%     recording  file name, for the messages.
%
%   Outputs
%     ok    true if no FAIL records (WARN is allowed through).
%     msgs  cellstr of human-readable lines for the app console (FAIL + WARN
%           only; empty if everything PASSed).
%     recs  the raw qc_check records, if you want to embed them in d.qc.

    if nargin < 3, recording = ""; end

    recs = validate_schema(d, stage, recording);

    isFail = arrayfun(@(r) r.status == "FAIL", recs);
    ok = ~any(isFail);

    % Surface only FAILs (real structural problems) to the app console. The only
    % WARNs validate_schema emits are "optional field absent" — by definition not
    % a problem (e.g. a recording with no analog/pacing channel), so they would
    % just be noise. They remain in `recs` for anyone who wants the full log.
    msgs = {};
    flagged = recs(isFail);
    for k = 1:numel(flagged)
        msgs{end+1, 1} = sprintf('SCHEMA %s [%s/%s]: %s', ...
            flagged(k).status, stage, flagged(k).check, flagged(k).message); %#ok<AGROW>
    end
end
