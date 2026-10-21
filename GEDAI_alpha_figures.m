function varargout = GEDAI_alpha_figures(stage, varargin)
% GEDAI_ALPHA_FIGURES  Broadband-GEDAI alpha figure pipeline: options and stage dispatch.
%
%   opts = GEDAI_alpha_figures('opts')   options struct; creates no directories, reads no data
%   GEDAI_alpha_figures(stage, ...)      run one stage (filled in by later tasks)
%
%   Every GEDAI call uses opts.gedaiArgs with refCOV (slot 4) filled by the caller:
%   GEDAI(EEG, opts.gedaiArgs{:}).

    if nargin < 1
        stage = 'opts';
    end

    switch lower(stage)
        case 'opts'
            varargout{1} = buildOpts();
        otherwise
            error('GEDAI_alpha:notImplemented', 'Stage "%s" is not implemented yet.', stage);
    end
end

function opts = buildOpts()
    home = getenv('HOME');
    opts.root = fullfile(home, 'scratch', 'gedai_review', 'alpha_figures');
    opts.cacheDir = fullfile(opts.root, 'cache');
    opts.figDir = fullfile(opts.root, 'figures');
    opts.eyeocDir = fullfile(opts.root, 'eyeOC');

    processedDir = fullfile(home, 'TMS_response_variability', 'store', 'mdd_dlpfc', ...
        'preprocessed_Oct26', 'rest_treat_joined', 'RELAXProcessed');
    opts.cleanedDir = fullfile(processedDir, 'Cleaned_Data');
    opts.extremesDir = fullfile(processedDir, 'Extremes_Rejected');
    opts.leadfieldDir = fullfile(home, 'TMS_response_variability', 'scratch', ...
        'eeg_sources_data', 'full_leadfield');

    opts.gedaiArgs = {'auto', 300, 4.0, [], false, false, [], [], [], [], 'auto-', 1.0, [], true};

    opts.edgeSec = 0.3;
    opts.minSegSec = 2;
    opts.winSec = 2;
    opts.maxHz = 200;
    opts.restTrainPadSec = 8;

    opts.bands = [1 4; 4 8; 8 13; 13 30; 30 45; 55 200];
    opts.bandNames = {'delta'; 'theta'; 'alpha'; 'beta'; 'gamma'; 'high'};
    opts.alphaBand = [8 13];
    opts.iafRange = [7 14];
    opts.nTopAlpha = 20;

    opts.injectScales = [0.25 0.5 1 2 4];

    opts.blinkHalfSec = 0.5;
    opts.muscleTopFraction = 0.10;
    opts.blinkElectrodes = {'E18'; 'E37'; 'E19'; 'E33'; 'E32'; 'E25'; 'E31'; 'E26'};

    opts.eyeocLeadfields = {'m121', 'm144'};
end
