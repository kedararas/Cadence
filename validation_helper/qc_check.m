function rec = qc_check(recording, stage, check, status, value, message)
%QC_CHECK  Build one standardised validation record.
%
%   rec = qc_check(recording, stage, check, status, value, message)
%
%   Every per-stage validator emits an array of these.  Keeping one fixed
%   shape lets all records from a whole batch be concatenated and turned into
%   a single table by qc_summarize.
%
%   Inputs
%     recording  identifier (file name) being validated.
%     stage      'conditioning' | 'features' | 'analysis' | ...
%     check      short name of the specific check, e.g. "CAM1_average".
%     status     "PASS" | "WARN" | "FAIL".
%                  PASS — output trustworthy.
%                  WARN — proceed, but flag for human review.
%                  FAIL — output is not trustworthy; do not use downstream.
%     value      numeric value behind the check (NaN if not applicable).
%     message    short human-readable explanation.
%
%   Output
%     rec        scalar struct with the fields above (strings + double).

    if nargin < 5, value   = NaN; end
    if nargin < 6, message = "";  end

    status = upper(string(status));
    assert(ismember(status, ["PASS","WARN","FAIL"]), ...
        'qc_check:badStatus', 'status must be PASS, WARN, or FAIL.');

    rec = struct( ...
        'recording', string(recording), ...
        'stage',     string(stage), ...
        'check',     string(check), ...
        'status',    status, ...
        'value',     double(value), ...
        'message',   string(message));
end
