% RELAX_residual_screen() - Residual-noise screen for the cleaned files of a RELAX_Wrapper_slurm.m run.
%
% Usage:
%   >> RELAX_residual_screen('score')   % score every cleaned file, write residual_screen.mat/.csv and README.md
%   >> RELAX_residual_screen('sheets')  % draw review sheets for the files to check by eye (needs 'score' first)
%   >> RELAX_residual_screen('all')     % score, then sheets
%   >> RELAX_residual_screen('sheets', {'m144_dlpfc_day1_stitched'})   % only these files
%   >> RELAX_residual_screen('all', {}, 'RELAXProcessed_muscleAfter')   % another run's output folder
%
% Output goes to <myPath>/<outputFolderName>/Residual_screen/.
% Edit myPath, outputFolderName, cutoff, holeDepthLimit, sampleSize, seed and the sheet settings below for a different dataset.

function RELAX_residual_screen(stage, onlyFiles, outputFolderName)

    myPath = '/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Oct26/rest_treat_joined/';
    if nargin < 3, outputFolderName = 'RELAXProcessed_broadbandGEDAI'; end % the RELAX_cfg.OutputFolderName of the run being screened
    cutoff = 3;
    holeDepthLimit = 1.4; % in electrode spacings
    sampleSize = 10;
    seed = 20261022;
    windowSeconds = 3;
    nWindows = 3;
    blinkHalfSeconds = 0.5;
    topoHalfSeconds = 0.05;

    if nargin < 1, stage = 'all'; end
    stage = validatestring(stage, {'score', 'sheets', 'all'});
    if nargin < 2, onlyFiles = {}; end

    outDir = fullfile(myPath, outputFolderName, 'Residual_screen');
    opts = struct('myPath', myPath, 'outputFolderName', outputFolderName, 'outDir', outDir, ...
                  'cutoff', cutoff, 'holeDepthLimit', holeDepthLimit, 'sampleSize', sampleSize, 'seed', seed, ...
                  'windowSeconds', windowSeconds, 'nWindows', nWindows, ...
                  'blinkHalfSeconds', blinkHalfSeconds, 'topoHalfSeconds', topoHalfSeconds);

    if any(strcmp(stage, {'score', 'all'}))
        RELAX_screen_run_score(opts);
    end

    if any(strcmp(stage, {'sheets', 'all'}))
        RELAX_screen_run_sheets(opts, onlyFiles);
    end
end
