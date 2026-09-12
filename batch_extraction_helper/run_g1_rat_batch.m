%RUN_G1_RAT_BATCH  Re-extract the G1 rat circadian study with the current algorithms.
%
%   Edit the three settings below, then run this script (or from a shell):
%     /Applications/MATLAB_R2025b.app/bin/matlab -nodisplay -batch "run('batch_extraction_helper/run_g1_rat_batch.m')"
%
%   TEST_FOLDERS = {'ZT2_8am/R2'}  processes one rat only.
%   TEST_FOLDERS = {}              processes the whole study.  Re-running
%   resumes: recordings already in cadence_recording_medians.csv whose
%   metrics file exists are skipped, so the run can be stopped and restarted.

INPUT_ROOT   = '/Volumes/Dayo_Aras /G1_Rat_Circadian_Mat';
OUTPUT_ROOT  = '/Volumes/Dayo_Aras /G1_Rat_Circadian_Metrics';
TEST_FOLDERS = {'ZT2_8am/R2'};          % {} = everything

opts = struct();
opts.FOV_mm = 20;                        % field of view used for CV (mm); the value stored in the old metrics files

cadence_batch_paths();
T = cadence_batch_extract(INPUT_ROOT, OUTPUT_ROOT, 'Folders', TEST_FOLDERS, 'ExtractOpts', opts);
disp(T(:, {'ZT','Experiment','Condition','CL_ms','apd_v','catd_ca','cv_v','df_v','Errors'}));
