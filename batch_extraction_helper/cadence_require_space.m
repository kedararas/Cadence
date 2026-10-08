function free = cadence_require_space(folder, need_bytes, what)
%CADENCE_REQUIRE_SPACE  Free space on the volume holding folder; error if short.
%
%   free = cadence_require_space(folder)                      % bytes free (NaN if unknown)
%   cadence_require_space(folder, need_bytes, 'metrics file') % throws when short
%
%   Throws 'cadence:diskFull' when fewer than need_bytes are free, so the batch
%   can stop the run with one clear message instead of failing every remaining
%   recording.  A save to a full disk fails inside MATLAB with messages such as
%   "Unable to write to file ... because it appears to be corrupt" or "No space
%   left on device", which do not say what is wrong.

    free = NaN;
    folder = char(folder);
    while ~isempty(folder) && ~isfolder(folder)      % folder not created yet: ask its parent
        parent = fileparts(folder);
        if strcmp(parent, folder), break; end
        folder = parent;
    end
    try
        free = double(java.io.File(folder).getUsableSpace());
    catch
    end
    if nargin < 2 || isempty(need_bytes) || isnan(free), return; end
    if nargin < 3, what = 'file'; end
    if free < need_bytes
        error('cadence:diskFull', ['Not enough disk space for the %s: %.1f GB needed, %.1f GB free on the ' ...
              'volume holding %s. Free space and run again (Resume continues where it stopped).'], ...
              what, need_bytes / 1e9, free / 1e9, folder);
    end
end
