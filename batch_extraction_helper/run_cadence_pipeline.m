function T = run_cadence_pipeline(raw_root, processed_root, metrics_root, varargin)
%RUN_CADENCE_PIPELINE  Raw camera files -> converted -> conditioned -> metrics + summary.
%
%   T = run_cadence_pipeline(raw_root, processed_root, metrics_root)
%   T = run_cadence_pipeline(raw_root, processed_root, metrics_root, 'Folders', {'ZT2_8am/R2'})
%
%   Thin wrapper around cadence_batch_pipeline: puts the CADENCE helpers on
%   the path, then for every recording folder under <raw_root> (any depth;
%   .tif volumes or SciMedia .gsh/.gsd pairs) converts it to
%   <processed_root>/converted/<rel>/<name>.mat, conditions it with the
%   Signal Conditioning app's default settings (SVD denoising, [0 50] Hz
%   temporal filter, drift correction, normalization, ensemble averaging) to
%   <processed_root>/conditioned/<rel>/<name>-conditioned.mat, extracts
%   features with the adaptive SNR mask to <metrics_root>/<rel>/<name>-metrics.mat,
%   and writes <metrics_root>/CADENCE_median_summary.xlsx.  Any
%   cadence_batch_pipeline option can be passed through (Folders, DryRun,
%   Recondition, ConditionOpts, ExtractOpts, ...).
%
%   Shell one-liner (no MATLAB desktop needed):
%     matlab -nodisplay -batch "addpath('batch_extraction_helper'); run_cadence_pipeline('/raw', '/processed', '/metrics')"

    addpath(fileparts(mfilename('fullpath')));
    cadence_batch_paths();
    T = cadence_batch_pipeline(raw_root, processed_root, metrics_root, varargin{:});
end
