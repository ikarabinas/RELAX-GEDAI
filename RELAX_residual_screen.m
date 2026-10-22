% RELAX_residual_screen() - Residual-noise screen for the cleaned files of a RELAX_Wrapper_slurm.m run.
%
% Usage:
%   >> RELAX_residual_screen('score')   % score every cleaned file, write residual_screen.mat/.csv and README.md
%   >> RELAX_residual_screen('all')     % every stage
%
% Output goes to <myPath>/<outputFolderName>/Residual_screen/.
% Edit myPath, outputFolderName, cutoff, sampleSize and seed below for a different dataset.

function RELAX_residual_screen(stage, onlyFiles)

    myPath = '/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Oct26/rest_treat_joined/';
    outputFolderName = 'RELAXProcessed_broadbandGEDAI'; % the RELAX_cfg.OutputFolderName of the run being screened
    cutoff = 3;
    sampleSize = 10;
    seed = 20261022;

    if nargin < 1, stage = 'all'; end
    if nargin < 2, onlyFiles = {}; end %#ok<NASGU>

    outDir = fullfile(myPath, outputFolderName, 'Residual_screen');
    opts = struct('myPath', myPath, 'outputFolderName', outputFolderName, 'outDir', outDir, ...
                  'cutoff', cutoff, 'sampleSize', sampleSize, 'seed', seed);

    if any(strcmp(stage, {'score', 'all'}))
        RELAX_screen_run_score(opts);
    end
end
