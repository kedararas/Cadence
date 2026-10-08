function ok = cadence_write_table(T, file, logf)
%CADENCE_WRITE_TABLE  writetable that cannot destroy the file it replaces.
%
%   ok = cadence_write_table(T, file)
%   ok = cadence_write_table(T, file, logf)     % logf(msg) for the warning
%
%   The table is written to <file>.tmp and renamed over <file> only when the
%   write succeeded and produced a non-empty file.  A plain writetable on a
%   full disk truncates its target to 0 bytes: that wiped the medians CSV of a
%   two-day batch (2026-10-08), and the next run, finding no rows, started
%   extracting everything again.  On failure the previous file is kept, the
%   .tmp is removed, a warning is logged and ok is false.

    if nargin < 3 || isempty(logf), logf = @(s) warning('cadence:writeTable', '%s', s); end
    file = char(file);
    [folder, name, ext] = fileparts(file);
    tmp = fullfile(folder, [name '.tmp' ext]);      % keeps the extension writetable needs
    ok = false;
    try
        writetable(T, tmp);
        d = dir(tmp);
        if isempty(d) || (d.bytes == 0 && (height(T) > 0 || width(T) > 0))
            error('cadence:writeTable', 'the file came out empty');
        end
        movefile(tmp, file, 'f');
        ok = true;
    catch ME
        if isfile(tmp), try, delete(tmp); catch, end, end
        free = cadence_require_space(folder);
        logf(sprintf('WARNING: could not write %s (%s); the previous file is kept. Free space there: %.1f GB.', ...
            file, ME.message, free / 1e9));
    end
end
