function T = run_cadence_batch(input_root, output_root, varargin)
%RUN_CADENCE_BATCH  Re-extract every conditioned recording under a parent folder.
%
%   T = run_cadence_batch(input_root, output_root)
%   T = run_cadence_batch(input_root, output_root, 'Folders', {'ZT2_8am/R2'})
%
%   Thin wrapper around cadence_batch_extract: puts the CADENCE helpers on
%   the path, then processes <input_root>/**/*.mat into
%   <output_root>/<same relative folder>/<name>-metrics.mat and writes
%   <output_root>/CADENCE_median_summary.xlsx.  Any cadence_batch_extract
%   option can be passed through (Folders, Resume, DryRun, ExtractOpts, ...).
%
%   Shell one-liner (no MATLAB desktop needed):
%     matlab -nodisplay -batch "addpath('batch_extraction_helper'); run_cadence_batch('/in', '/out')"

    addpath(fileparts(mfilename('fullpath')));
    cadence_batch_paths();
    T = cadence_batch_extract(input_root, output_root, varargin{:});
end
