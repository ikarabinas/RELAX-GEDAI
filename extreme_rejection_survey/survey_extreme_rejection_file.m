function survey_extreme_rejection_file(runDir, fileIndex, outDir, kOverride)
% Re-runs RELAX's extreme-period marking on one <runDir>/Extremes_Rejected file and
% records, per epoch, which channels tripped which criterion. Then asks how much data
% each alternative epoch-rejection rule would delete. Writes <outDir>/<file>.mat.
%
% The saved file is continuousEEG exactly as RELAX_excluding_extreme_values saw it
% (channels already removed, filtering applied), plus the blink events and RELAX fields
% added afterwards; those are stripped so the blink detector sees what it saw then.
% kOverride (optional) reruns RELAX with that ExtremeMinBadChannelsToRejectEpoch instead
% of the run's, to check RELAX's own k rule against this survey's prediction of it.

    files = dir(fullfile(runDir, 'Extremes_Rejected', '*_Extremes_Rejected.set'));
    [~, order] = sort({files.name});
    files = files(order);
    file = files(fileIndex);
    base = erase(file.name, '_Extremes_Rejected.set');
    fprintf('Surveying %d/%d: %s\n', fileIndex, numel(files), base);

    loaded = load(fullfile(runDir, 'RELAX_cfg', [base '.mat']), 'RELAX_cfg');
    cfg = loaded.RELAX_cfg;
    runK = 1;  % the run's ExtremeMinBadChannelsToRejectEpoch; absent = RELAX's any-one-channel rule
    if isfield(cfg, 'ExtremeMinBadChannelsToRejectEpoch'); runK = cfg.ExtremeMinBadChannelsToRejectEpoch; end
    if nargin > 3 && ~isempty(kOverride)
        runK = kOverride;
        cfg.ExtremeMinBadChannelsToRejectEpoch = kOverride;
    end

    EEG = pop_loadset('filename', file.name, 'filepath', file.folder);
    saved.proportionExcluded = EEG.RELAXProcessingExtremeRejections.ProportionExcludedForExtremeOutlier;
    saved.proportionDeleted = deleted_fraction(EEG.RELAX.ExtremelyBadPeriodsForDeletion, EEG.pnts);

    EEG.event(cellfun(@(t) ischar(t) && startsWith(t, 'EyeBlink'), {EEG.event.type})) = [];
    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG = rmfield(EEG, intersect(fieldnames(EEG), ...
        {'RELAX', 'RELAXProcessing', 'RELAXProcessingExtremeRejections'}));
    cfg.ms_per_sample = 1000 / EEG.srate;

    channelsAlreadyScreened = cfg;
    channelsAlreadyScreened.MaxProportionOfElectrodesThatCanBeDeleted = 0;
    if cfg.RespectBoundaries == 1
        [continuousEEG, epochedEEG] = RELAX_excluding_channels_and_epoching_boundary(EEG, channelsAlreadyScreened);
    else
        [continuousEEG, epochedEEG] = RELAX_excluding_channels_and_epoching(EEG, channelsAlreadyScreened);
    end
    clear EEG
    % RELAX_blinks_IQR_method reads RELAX.ExtremelyBadPeriodsForDeletion when present, so
    % the blink-aware shift bound must be rebuilt from the data before that is set:
    continuousBeforeMarking = continuousEEG;
    [continuousEEG, epochedEEG] = RELAX_excluding_extreme_values(continuousEEG, epochedEEG, cfg);

    flags = extreme_flags_per_channel(continuousBeforeMarking, epochedEEG, cfg);
    clear continuousBeforeMarking
    nEpochs = size(epochedEEG.data, 3);
    nChan = epochedEEG.nbchan;
    % kurtSingle is left out: RELAX_excluding_extreme_values computes it but only the
    % whole-montage kurtosis reaches AllMethodsExtremeEpochRejections (its lines 176-180),
    % so single-channel kurtosis never deletes data. It is reported separately below.
    criteria = {'shiftHigh', 'shiftLow', 'absVolt', 'jointProb', 'drift'};
    anyChannel = false(nChan, nEpochs);
    for c = 1:numel(criteria)
        anyChannel = anyChannel | flags.(criteria{c});
    end
    badChannelsPerEpoch = sum(anyChannel, 1);

    d = epochedEEG.RELAXProcessing.Details;
    relaxEpochVectors = struct( ...
        'shift', padded(d.EpochsExceedingShiftThreshold, nEpochs), ...
        'absVolt', padded(field_or_zero(d, 'AbsoluteAmplitudeExceededThreshold'), nEpochs), ...
        'kurt', padded(d.MarksForBothSingleAndAllChannelKurt, nEpochs), ...
        'jointProb', padded(d.ImprobableVoltageDistributionExceededThreshold, nEpochs), ...
        'drift', padded(field_or_zero(d, 'driftslopemaskepochs'), nEpochs), ...
        'all', padded(d.AllMethodsExtremeEpochRejections, nEpochs));
    mine = struct( ...
        'shift', any(flags.shiftHigh | flags.shiftLow, 1), ...
        'absVolt', any(flags.absVolt, 1), ...
        'kurt', any(flags.kurtSingle, 1) | flags.kurtAllEpoch, ...
        'jointProb', any(flags.jointProb, 1), ...
        'drift', any(flags.drift, 1), ...
        'all', badChannelsPerEpoch >= runK | flags.kurtAllEpoch);
    check.mismatchedEpochs = structfun(@(x) 0, mine, 'UniformOutput', false);
    for name = fieldnames(mine)'
        check.mismatchedEpochs.(name{1}) = sum(mine.(name{1}) ~= relaxEpochVectors.(name{1}));
    end
    check.rerunProportionExcluded = epochedEEG.RELAXProcessingExtremeRejections.ProportionExcludedForExtremeOutlier;
    check.savedProportionExcluded = saved.proportionExcluded;
    check.savedProportionDeleted = saved.proportionDeleted;

    rule = struct('name', {}, 'k', {}, 'rejectEpochs', {});
    for k = [1 2 3 4 5 10 20]
        rule(end+1) = struct('name', sprintf('k>=%d', k), 'k', k, ...
            'rejectEpochs', badChannelsPerEpoch >= k | flags.kurtAllEpoch); %#ok<AGROW>
    end
    for pct = [1 2 5 10]
        k = ceil(pct / 100 * nChan);
        rule(end+1) = struct('name', sprintf('>=%d%%(k=%d)', pct, k), 'k', k, ...
            'rejectEpochs', badChannelsPerEpoch >= k | flags.kurtAllEpoch); %#ok<AGROW>
    end
    segments = RELAX_get_segments(continuousEEG);
    windows = tms_windows(continuousEEG, segments);
    for r = 1:numel(rule)
        [rule(r).proportionExcluded, rule(r).proportionDeleted, mask] = data_loss(rule(r).rejectEpochs, epochedEEG, continuousEEG);
        rule(r).geometry = rejection_geometry(mask, segments, windows, epochedEEG.srate);
    end
    neverEvaluated = isnan(mask);  % the same for every rule
    windowSummary = struct('usedTrainIdx', windows.usedTrainIdx, ...
        'nPost', size(windows.post, 1), 'nPre', size(windows.pre, 1), ...
        'nDroppedPost', windows.nDroppedPost, 'nDroppedPre', windows.nDroppedPre, ...
        'nPostNeverEvaluated', count_windows_touching(neverEvaluated, windows.post), ...
        'nPreNeverEvaluated', count_windows_touching(neverEvaluated, windows.pre));
    rule = rmfield(rule, 'rejectEpochs');
    sameRule = startsWith({rule.name}, 'k') & [rule.k] == runK;
    check.ruleMatchesRerun = any(sameRule) && abs(rule(sameRule).proportionExcluded - check.rerunProportionExcluded) < 1e-9;

    evaluated = ~flags.kurtAllEpoch;  % epochs the per-channel rules can act on
    byCriterion = struct();
    for c = 1:numel(criteria)
        epochHit = any(flags.(criteria{c}), 1);
        otherHit = false(1, nEpochs);
        for o = setdiff(1:numel(criteria), c)
            otherHit = otherHit | any(flags.(criteria{o}), 1);
        end
        byCriterion.(criteria{c}) = struct( ...
            'proportionEpochsFlagged', mean(epochHit), ...
            'proportionEpochsFlaggedOnlyByThis', mean(epochHit & ~otherHit & evaluated), ...
            'channelEpochFlags', nnz(flags.(criteria{c})));
    end
    byCriterion.kurtSingleUnused = struct('proportionEpochsFlagged', mean(any(flags.kurtSingle, 1)), ...
        'proportionEpochsFlaggedOnlyByThis', NaN, 'channelEpochFlags', nnz(flags.kurtSingle));
    byCriterion.kurtAll = struct('proportionEpochsFlagged', mean(flags.kurtAllEpoch), ...
        'proportionEpochsFlaggedOnlyByThis', mean(flags.kurtAllEpoch & badChannelsPerEpoch == 0), ...
        'channelEpochFlags', NaN);

    flagsPerChannel = sum(anyChannel, 2);
    [~, worstFirst] = sort(flagsPerChannel, 'descend');
    dropWorst = struct('nChannelsDropped', {}, 'proportionExcluded', {}, 'geometry', {});
    for n = [1 2 5 10 20]
        keep = true(nChan, 1);
        keep(worstFirst(1:n)) = false;
        reject = any(anyChannel(keep, :), 1) | flags.kurtAllEpoch;
        [excluded, ~, mask] = data_loss(reject, epochedEEG, continuousEEG);
        dropWorst(end+1) = struct('nChannelsDropped', n, 'proportionExcluded', excluded, ...
            'geometry', rejection_geometry(mask, segments, windows, epochedEEG.srate)); %#ok<AGROW>
    end

    result = struct('file', base, 'nChan', nChan, 'nEpochs', nEpochs, 'srate', epochedEEG.srate, ...
        'minutes', continuousEEG.pnts / continuousEEG.srate / 60, ...
        'check', check, 'rule', rule, 'byCriterion', byCriterion, 'dropWorst', dropWorst, ...
        'windows', windowSummary, 'segments', segments, ...
        'epochOnsets', [epochedEEG.event.originallatency], ...
        'badChannelsPerEpoch', uint16(badChannelsPerEpoch), 'kurtAllEpoch', flags.kurtAllEpoch, ...
        'flagsPerChannel', flagsPerChannel, 'channelLabels', {{epochedEEG.chanlocs.labels}});
    if ~exist(outDir, 'dir'); mkdir(outDir); end
    save(fullfile(outDir, [base '.mat']), 'result');
    fprintf('Saved %s: RELAX rerun excluded %.4f (saved %.4f); k>=1 %.4f, k>=5 %.4f\n', base, ...
        check.rerunProportionExcluded, check.savedProportionExcluded, ...
        rule(1).proportionExcluded, rule([rule.k] == 5 & startsWith({rule.name}, 'k')).proportionExcluded);
    disp(check.mismatchedEpochs)
end

function [proportionExcluded, proportionDeleted, mask] = data_loss(rejectEpochs, epochedEEG, continuousEEG)
% Same sample mask RELAX_excluding_extreme_values builds (lines 243-273): evaluated
% samples start at 0, never-evaluated at NaN, rejected epochs set to 1. Everything NaN
% or 1 is later deleted by eeg_eegrej.
    mask = epochedEEG.RELAXProcessing.Details.NaNsForNonEvents;
    epochLength = round(epochedEEG.srate);
    for e = find(rejectEpochs)
        onset = epochedEEG.event(e).originallatency;
        mask(onset:onset + epochLength - 1) = 1;
    end
    segments = RELAX_get_segments(continuousEEG);
    segmentEnds = segments(segments(:,2) > segments(:,1), 2)';
    mask(segmentEnds(mask(segmentEnds - 1) == 1)) = 1;
    proportionExcluded = mean(mask, 'omitnan');
    proportionDeleted = mean(isnan(mask) | mask == 1);
end

function fraction = deleted_fraction(periods, nSamples)
    if isempty(periods)
        fraction = 0;
    else
        fraction = sum(periods(:,2) - periods(:,1) + 1) / nSamples;
    end
end

function v = padded(v, n)
% RELAX grows these vectors only as far as the last flagged epoch.
    v = logical([v(:)', zeros(1, n - numel(v))]);
end

function v = field_or_zero(s, name)
    if isfield(s, name); v = s.(name); else; v = 0; end
end
