function root = cadence_batch_paths()
%CADENCE_BATCH_PATHS  Put the CADENCE helper folders on the MATLAB path.
%   Adds only the source helper folders (not distribution/release copies).
    root = fileparts(fileparts(mfilename('fullpath')));
    sub = {'batch_extraction_helper', 'feature_extraction_helper', ...
           'conduction_velocity_helper', 'signal_conditioning_helper', ...
           'arrhythmia_dynamics_helper', fullfile('arrhythmia_dynamics_helper','objects'), ...
           fullfile('arrhythmia_dynamics_helper','frechet'), 'utils', ...
           'validation_helper', 'manual_validation_helper'};
    for k = 1:numel(sub)
        p = fullfile(root, sub{k});
        if isfolder(p), addpath(p); end
    end
end
