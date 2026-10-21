function varargout = GEDAI_alpha_figures(stage, varargin)
% GEDAI_ALPHA_FIGURES  Broadband-GEDAI alpha figure pipeline: options and stage dispatch.
%
%   opts = GEDAI_alpha_figures('opts')   options struct; creates no directories, reads no data
%   GEDAI_alpha_figures('list')          print the numbered eligible sessions and their count N
%   GEDAI_alpha_figures('status')        count cached sessions by status; name the ones not ok
%   GEDAI_alpha_figures('extract')       summarise session SLURM_ARRAY_TASK_ID (submit --array=1-N)
%   GEDAI_alpha_figures(stage, ...)      other stages are filled in by later tasks
%
%   Every GEDAI call uses opts.gedaiArgs with refCOV (slot 4) filled by the caller:
%   GEDAI(EEG, opts.gedaiArgs{:}).

    if nargin < 1
        stage = 'opts';
    end

    switch lower(stage)
        case 'opts'
            varargout{1} = buildOpts();
        case 'list'
            listSessions(buildOpts(), true);
        case 'status'
            printStatus(buildOpts());
        case 'extract'
            extractSession(buildOpts());
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

function names = listSessions(opts, doPrint)
    cleaned = dir(fullfile(opts.cleanedDir, '*_RELAX.set'));
    candidates = sort(erase({cleaned.name}, '_RELAX.set'));
    hasExtremes = cellfun(@(n) isfile(fullfile(opts.extremesDir, [n '_Extremes_Rejected.set'])), candidates);
    hasLeadfield = cellfun(@(n) ~isempty(dir(fullfile(opts.leadfieldDir, [strtok(n, '_') '_*.mat']))), candidates);
    names = candidates(hasExtremes & hasLeadfield);
    if doPrint
        for k = 1:numel(names)
            fprintf('%3d  %s\n', k, names{k});
        end
        fprintf('N = %d\n', numel(names));
    end
end

function extractSession(opts)
    names = listSessions(opts, false);
    index = str2double(getenv('SLURM_ARRAY_TASK_ID'));
    if isnan(index) || index < 1 || index > numel(names)
        error('GEDAI_alpha:badTaskId', 'SLURM_ARRAY_TASK_ID %s is outside 1-%d', getenv('SLURM_ARRAY_TASK_ID'), numel(names));
    end
    name = names{index};
    cacheFile = fullfile(opts.cacheDir, ['session_' name '.mat']);
    if isfile(cacheFile)
        fprintf('%s already cached, skipping.\n', name);
        return
    end
    if ~isfolder(opts.cacheDir)
        mkdir(opts.cacheDir);
    end

    try
        S = GEDAI_alpha_session_summary(name, opts);
    catch err
        S = struct('name', name, 'participant', strtok(name, '_'), 'status', 'error', 'reason', err.message);
        fprintf(2, '%s\n', getReport(err));
    end
    save(cacheFile, 'S', '-v7.3');
    fprintf('%s: status %s\n', name, S.status);
    if strcmp(S.status, 'ok')
        fprintf(['alphaRatioTop %.4f  enova %.4f  threshold %.4f  gofInput %.4f  gofRemoved %.4f\n' ...
                 'corrBlink %.4f  corrMuscle %.4f  nExcluded %d  nSamplesUsed %d\n'], ...
            S.B.alphaRatioTop, S.enova, S.threshold, S.gofInput, S.gofRemoved, ...
            S.tmpl.corrBlink, S.tmpl.corrMuscle, S.M.nExcluded, S.M.nSamplesUsed);
    else
        fprintf('reason: %s\n', S.reason);
    end
end

function printStatus(opts)
    files = dir(fullfile(opts.cacheDir, 'session_*.mat'));
    statuses = cell(1, numel(files));
    reasons = cell(1, numel(files));
    for k = 1:numel(files)
        loaded = load(fullfile(files(k).folder, files(k).name), 'S');
        statuses{k} = loaded.S.status;
        reasons{k} = loaded.S.reason;
    end
    for status = {'ok', 'skipped', 'error'}
        fprintf('%-8s %d\n', status{1}, sum(strcmp(statuses, status{1})));
    end
    fprintf('cached   %d of %d eligible sessions\n', numel(files), numel(listSessions(opts, false)));
    for k = find(~strcmp(statuses, 'ok'))
        fprintf('  %s [%s]: %s\n', erase(files(k).name, {'session_', '.mat'}), statuses{k}, reasons{k});
    end
end
