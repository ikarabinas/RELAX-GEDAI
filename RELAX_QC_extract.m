% RELAX_QC_extract() - Stage 1 of the QC figure set: summarise every cleaned file, pick
%                      Figure 1 examples, extract their windows, and save the cache.
%
% Usage:
%   >> cache = RELAX_QC_extract(opts);
%
% opts: built by RELAX_QC_figures (.myPath, .cachePath, .overrides, .onlyFiles and the
% analysis settings). A file that fails is recorded as "skipped" with its error message;
% the run continues.

function cache = RELAX_QC_extract(opts)

    processed = fullfile(opts.myPath, opts.outputFolderName);
    [T, cfgs] = RELAX_QC_load_tables(fullfile(processed, 'PerParticipantStats'), fullfile(processed, 'Cleaned_Data'));
    cellfun(@RELAX_QC_check_config, values(cfgs));
    if ~isempty(opts.onlyFiles)
        T = T(ismember(T.File, string(opts.onlyFiles)), :);
        if isempty(T)
            error('RELAX_QC:noFilesMatched', 'None of the requested files were found: %s', strjoin(string(opts.onlyFiles), ', '));
        end
    end
    inputs = dir(fullfile(opts.myPath, '**', '*.set'));
    inputs = inputs(~contains({inputs.folder}, 'RELAXProcessed'));

    n = height(T);
    summaries = repmat(RELAX_QC_empty_summary(), n, 1);
    todo = find(T.Outcome == "cleaned")';
    files = T.File(todo);
    fileCfgs = arrayfun(@(i) cfgs(char(T.Participant(i))), todo, 'UniformOutput', false);
    startWorkerPool(numel(todo));
    results = cell(numel(todo), 1);
    parfor k = 1:numel(todo)
        results{k} = summariseOneFile(processed, inputs, files(k), fileCfgs{k}, opts);
    end
    for k = 1:numel(todo)
        summaries(todo(k)) = results{k};
        fprintf('RELAX_QC_extract %d/%d %s: %s %s\n', todo(k), n, files(k), results{k}.status, results{k}.skipReason);
    end

    unusable = [summaries.status] ~= "ok";
    scores.blinkAmplitude = arrayfun(@(s) s.scores.blinkAmplitude, summaries);
    scores.musclePower    = arrayfun(@(s) s.scores.musclePower, summaries);
    scores.popSteps       = arrayfun(@(s) s.scores.popSteps, summaries);
    scores.blinkAmplitude(unusable) = NaN; scores.musclePower(unusable) = NaN; scores.popSteps(unusable) = NaN;
    picks = RELAX_QC_pick_examples(T.File, scores, opts.overrides);

    saveCache(opts, cfgs, T, summaries, picks, emptyExamples(), true);

    examples = RELAX_QC_build_examples(picks, T, summaries, cfgs, processed, inputs, opts);

    saveCache(opts, cfgs, T, summaries, picks, examples, false);
    fprintf('RELAX_QC_extract: %d ok, %d skipped, cache %s\n', sum(~unusable), sum([summaries.status] == "skipped"), opts.cachePath);
end

function examples = emptyExamples()
    examples = struct('file', {}, 'role', {}, 'windows', {});
end

function picks = emptyPicks()
    picks = struct('file', {}, 'role', {}, 'centres', {});
end

% meta.partial marks a checkpoint: the draw stage warns and draws what the cache holds.
function saveCache(opts, cfgs, T, summaries, picks, examples, partial)
    anyCfg = cfgs(char(T.Participant(1)));
    cache.meta = struct('created', datetime('now'), 'myPath', string(opts.myPath), 'opts', opts, ...
                        'muscleSlopeThreshold', anyCfg.MuscleSlopeThreshold, ...
                        'maxChannelProportion', anyCfg.MaxProportionOfElectrodesThatCanBeDeleted, ...
                        'minBlinksForGroup', opts.minBlinksForGroup, 'partial', partial);
    cache.tables = T;
    cache.summaries = summaries;
    cache.picks = picks;
    cache.examples = examples;
    save(opts.cachePath, 'cache', '-v7.3');
end

function summary = summariseOneFile(processed, inputs, file, cfg, opts)
    try
        summary = RELAX_QC_file_summary(RELAX_QC_file_paths(processed, inputs, file), cfg, opts);
    catch err
        summary = RELAX_QC_empty_summary();
        summary.status = "skipped";
        summary.skipReason = string(err.message);
        warning('RELAX_QC:fileSkipped', '%s skipped: %s', file, err.message);
    end
end

% Each SLURM job needs its own job storage, or concurrent pools clobber each other.
function startWorkerPool(nFiles)
    if ~isempty(gcp('nocreate')) || nFiles == 0, return, end
    nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));
    if isnan(nWorkers), return, end
    nWorkers = min(nFiles, max(1, nWorkers));
    cluster = parcluster('local');
    storage = fullfile(getenv('TMPDIR'), sprintf('qc_extract_%s', getenv('SLURM_JOB_ID')));
    if ~isempty(getenv('TMPDIR')) && ~isempty(getenv('SLURM_JOB_ID'))
        mkdir(storage);
        cluster.JobStorageLocation = storage;
    end
    parpool(cluster, nWorkers);
end
