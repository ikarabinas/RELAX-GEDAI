% RELAX_screen_run_sheets() - Sheets stage of the residual-noise screen: draw a review sheet for
%                             each file a person should check by eye.
%
% Usage:
%   >> screen = RELAX_screen_run_sheets(opts, onlyFiles);
%
% opts: built by RELAX_residual_screen (the score-stage fields plus .windowSeconds, .nWindows,
% .muscleSmoothSeconds, .blinkHalfSeconds and .topoHalfSeconds).
% onlyFiles: cell or string array of base names; empty means the files that fail, are
% unscreened or are in the review sample. An unknown name raises RELAX_screen:unknownFile.
%
% Loads <outDir>/residual_screen.mat, writes <outDir>/review/<File>.png/.pdf for each selected
% file, and rewrites residual_screen.mat and .csv with SheetStatus "ok" or "failed". A file that
% fails (SheetMessage = the error) does not stop the others. Rows that were not selected are
% left as they were.

function screen = RELAX_screen_run_sheets(opts, onlyFiles)

    if nargin < 2, onlyFiles = {}; end
    loaded = load(fullfile(opts.outDir, 'residual_screen.mat'), 'screen', 'meta');
    screen = loaded.screen;
    meta = loaded.meta;

    selected = find(selectedRows(screen, string(onlyFiles)))';
    if isempty(selected)
        fprintf('RELAX_screen_run_sheets: no files to draw\n');
        return
    end

    processed = fullfile(opts.myPath, opts.outputFolderName);
    [T, cfgs] = RELAX_QC_load_tables(fullfile(processed, 'PerParticipantStats'), fullfile(processed, 'Cleaned_Data'));
    reviewDir = fullfile(opts.outDir, 'review');

    n = numel(selected);
    files = screen.File(selected);
    fileCfgs = arrayfun(@(row) participantConfig(cfgs, screen.Participant(row)), selected, 'UniformOutput', false);
    fileOpts = arrayfun(@(file) withBlinkFallback(opts, T, file), files, 'UniformOutput', false);
    rows = arrayfun(@(row) screen(row, :), selected, 'UniformOutput', false);

    maxWorkers = startWorkerPool(n);
    outcomes = cell(n, 1);
    parfor (k = 1:n, maxWorkers)
        outcomes{k} = drawOneSheet(processed, rows{k}, screen, meta, fileCfgs{k}, fileOpts{k}, reviewDir);
    end

    for k = 1:n
        row = selected(k);
        screen.SheetStatus(row) = outcomes{k}.status;
        screen.SheetMessage(row) = messageAfter(screen.SheetStatus(row), screen.SheetMessage(row), outcomes{k});
        fprintf('RELAX_screen_run_sheets %d/%d %s: %s %s\n', k, n, files(k), outcomes{k}.status, outcomes{k}.message);
    end

    save(fullfile(opts.outDir, 'residual_screen.mat'), 'screen', 'meta');
    RELAX_screen_write_csv(screen, fullfile(opts.outDir, 'residual_screen.csv'));
    fprintf('RELAX_screen_run_sheets: %d ok, %d failed, sheets in %s\n', ...
            sum(screen.SheetStatus == "ok"), sum(screen.SheetStatus == "failed"), reviewDir);
end

function chosen = selectedRows(screen, onlyFiles)
    if isempty(onlyFiles)
        chosen = ismember(screen.ScreenOutcome, ["fail", "unscreened"]) | screen.InReviewSample;
        return
    end
    unknown = setdiff(onlyFiles, screen.File);
    if ~isempty(unknown)
        error('RELAX_screen:unknownFile', 'Not in the screen table: %s', strjoin(unknown, ', '));
    end
    chosen = ismember(screen.File, onlyFiles);
end

function cfg = participantConfig(cfgs, participant)
    cfg = [];
    if isKey(cfgs, char(participant)), cfg = cfgs(char(participant)); end
end

function fileOpts = withBlinkFallback(opts, T, file)
    fileOpts = opts;
    fileOpts.blinkElectrodeFallback = NaN;
    hit = find(T.File == file, 1);
    if ~isempty(hit), fileOpts.blinkElectrodeFallback = T.BlinkElectrodeFallback(hit); end
end

function outcome = drawOneSheet(processed, row, screen, meta, cfg, opts, reviewDir)
    outcome = struct('status', "ok", 'message', "");
    figuresBefore = findall(groot, 'Type', 'figure');
    try
        if isempty(cfg)
            error('RELAX_screen:noConfig', 'No saved RELAX_cfg for participant %s', row.Participant);
        end
        paths = sheetPaths(processed, row.File);
        sheet = RELAX_screen_sheet_data(paths, cfg, opts);
        RELAX_screen_draw_sheet(sheet, row, screen, meta, reviewDir);
    catch err
        close(setdiff(findall(groot, 'Type', 'figure'), figuresBefore));
        outcome = struct('status', "failed", 'message', string(err.message));
        warning('RELAX_screen:sheetFailed', '%s: %s', row.File, err.message);
    end
end

function paths = sheetPaths(processed, file)
    cleanedDir = fullfile(processed, 'Cleaned_Data');
    paths.cleaned  = char(fullfile(cleanedDir, file + RELAX_QC_cleaned_suffix(cleanedDir, file)));
    paths.extremes = char(fullfile(processed, 'Extremes_Rejected', file + "_Extremes_Rejected.set"));
end

% A failed message from an earlier run is stale once the file succeeds; a score-stage note stays.
function message = messageAfter(previousStatus, previousMessage, outcome)
    if outcome.status == "failed"
        message = outcome.message;
    elseif previousStatus == "failed"
        message = "";
    else
        message = previousMessage;
    end
end

% Each SLURM job needs its own job storage, or concurrent pools clobber each other. Returns the
% most workers parfor may use: 0 (run serially, no pool) outside SLURM, Inf with a pool already open.
function maxWorkers = startWorkerPool(nFiles)
    maxWorkers = Inf;
    if ~isempty(gcp('nocreate')), return, end
    maxWorkers = 0;
    nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));
    if nFiles == 0 || isnan(nWorkers), return, end
    nWorkers = min(nFiles, max(1, nWorkers));
    cluster = parcluster('local');
    storage = fullfile(getenv('TMPDIR'), sprintf('residual_screen_%s', getenv('SLURM_JOB_ID')));
    if ~isempty(getenv('TMPDIR')) && ~isempty(getenv('SLURM_JOB_ID'))
        mkdir(storage);
        cluster.JobStorageLocation = storage;
    end
    parpool(cluster, nWorkers);
    maxWorkers = nWorkers;
end
