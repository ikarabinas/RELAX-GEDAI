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

    processed = fullfile(opts.myPath, 'RELAXProcessed');
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
    for i = 1:n
        if T.Outcome(i) ~= "cleaned", continue, end
        try
            summaries(i) = RELAX_QC_file_summary(filePaths(processed, inputs, T.File(i)), cfgs(char(T.Participant(i))), opts);
        catch err
            summaries(i).status = "skipped";
            summaries(i).skipReason = string(err.message);
            warning('RELAX_QC:fileSkipped', '%s skipped: %s', T.File(i), err.message);
        end
        fprintf('RELAX_QC_extract %d/%d %s: %s %s\n', i, n, T.File(i), summaries(i).status, summaries(i).skipReason);
    end

    unusable = [summaries.status] ~= "ok";
    scores.blinkAmplitude = arrayfun(@(s) s.scores.blinkAmplitude, summaries);
    scores.musclePower    = arrayfun(@(s) s.scores.musclePower, summaries);
    scores.popSteps       = arrayfun(@(s) s.scores.popSteps, summaries);
    scores.blinkAmplitude(unusable) = NaN; scores.musclePower(unusable) = NaN; scores.popSteps(unusable) = NaN;
    picks = RELAX_QC_pick_examples(T.File, scores, opts.overrides);

    examples = struct('file', {}, 'role', {}, 'windows', {});
    for k = 1:numel(picks)
        i = find(T.File == picks(k).file, 1);
        if summaries(i).status ~= "ok"
            warning('RELAX_QC:exampleUnavailable', 'Example %s has status %s; not shown.', picks(k).file, summaries(i).status);
            continue
        end
        windows = RELAX_QC_example_windows(picks(k), summaries(i), filePaths(processed, inputs, T.File(i)), ...
                                           cfgs(char(T.Participant(i))), opts);
        examples(end+1) = struct('file', picks(k).file, 'role', picks(k).role, 'windows', windows); %#ok<AGROW>
    end

    anyCfg = cfgs(char(T.Participant(1)));
    cache.meta = struct('created', datetime('now'), 'myPath', string(opts.myPath), 'opts', opts, ...
                        'muscleSlopeThreshold', anyCfg.MuscleSlopeThreshold, ...
                        'maxChannelProportion', anyCfg.MaxProportionOfElectrodesThatCanBeDeleted, ...
                        'minBlinksForGroup', opts.minBlinksForGroup);
    cache.tables = T;
    cache.summaries = summaries;
    cache.picks = picks;
    cache.examples = examples;
    save(opts.cachePath, 'cache', '-v7.3');
    fprintf('RELAX_QC_extract: %d ok, %d skipped, cache %s\n', sum(~unusable), sum([summaries.status] == "skipped"), opts.cachePath);
end

function paths = filePaths(processed, inputs, file)
    hit = find(strcmp({inputs.name}, char(file + ".set")), 1);
    if isempty(hit), error('RELAX_QC:missingInput', 'No input file %s.set under myPath', file); end
    paths.input    = fullfile(inputs(hit).folder, inputs(hit).name);
    paths.extremes = char(fullfile(processed, 'Extremes_Rejected', file + "_Extremes_Rejected.set"));
    paths.cleaned  = char(fullfile(processed, 'Cleaned_Data', file + "_RELAX.set"));
end
