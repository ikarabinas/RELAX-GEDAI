function [EEG, covariates] = RELAX_second_pass_channel_epochs(EEG, locked, RELAX_cfg)
% Post-cleaning re-marking (RELAX_cfg.InterpolateExtremeChannelEpochsAfterCleaning = 1):
% epochs a copy of the cleaned data exactly as the first extreme pass did, flags the
% channel-epochs that are extreme by the first pass's per-channel statistics (locked, from
% RELAX_lock_extreme_thresholds), and replaces each with a spherical spline of the epoch's
% other channels (RELAX_interpolate_channel_epochs). Epochs never straddle a boundary, so
% neither does a repair. Design:
% docs/superpowers/specs/2026-11-06-post-cleaning-channel-epoch-interpolation-design.md
%
% Off (absent or 0), EEG is returned unchanged and every count covariate is NaN, so a
% table mixing runs can tell "off" from "found nothing".

    if ~isfield(RELAX_cfg, 'InterpolateExtremeChannelEpochsAfterCleaning')
        RELAX_cfg.InterpolateExtremeChannelEpochsAfterCleaning = 0;
    end
    if ~isfield(RELAX_cfg, 'SecondPassMaxProportionFlaggedChannels')
        RELAX_cfg.SecondPassMaxProportionFlaggedChannels = 0.20;
    end
    criteria = {'ShiftHigh', 'ShiftLow', 'AbsVolt', 'JointProb', 'Drift', 'Kurtosis'};
    covariates.InterpolateExtremeChannelEpochsAfterCleaning = RELAX_cfg.InterpolateExtremeChannelEpochsAfterCleaning;
    covariates.SecondPassMaxProportionFlaggedChannels = RELAX_cfg.SecondPassMaxProportionFlaggedChannels;
    covariates.SecondPassEpochsEvaluated = NaN;
    covariates.SecondPassEpochsWithInterpolatedChannels = NaN;
    covariates.SecondPassChannelEpochsInterpolated = NaN;
    covariates.SecondPassProportionChannelTimeInterpolated = NaN;
    for c = criteria
        covariates.(['SecondPassFlaggedBy' c{1}]) = NaN;
    end
    if RELAX_cfg.InterpolateExtremeChannelEpochsAfterCleaning ~= 1
        return
    end

    [epoched, onsets] = epoch_like_first_pass(EEG, RELAX_cfg);
    labels = {epoched.chanlocs.labels};
    scores = RELAX_score_locked_extremes(locked, epoched.data, labels);
    flags = RELAX_flag_extreme_channel_epochs(scores, epoched.RELAXProcessing.Details.DriftRatioEpochsxChannels, ...
        (onsets - 1) / EEG.srate, RELAX_cfg);
    [EEG, changed] = RELAX_interpolate_channel_epochs(EEG, flags.Any, onsets, epoched.pnts, 0.25);

    covariates.SecondPassEpochsEvaluated = epoched.trials;
    covariates.SecondPassEpochsWithInterpolatedChannels = nnz(any(flags.Any, 1));
    covariates.SecondPassChannelEpochsInterpolated = nnz(flags.Any);
    covariates.SecondPassProportionChannelTimeInterpolated = nnz(changed) / numel(changed);
    for c = criteria
        covariates.(['SecondPassFlaggedBy' c{1}]) = nnz(flags.(c{1}));
    end
    fprintf('Second-pass re-marking: %d channel-epochs in %d of %d epochs interpolated (%.3f%% of channel x time)\n', ...
        covariates.SecondPassChannelEpochsInterpolated, covariates.SecondPassEpochsWithInterpolatedChannels, ...
        covariates.SecondPassEpochsEvaluated, 100 * covariates.SecondPassProportionChannelTimeInterpolated);
end

function [epoched, onsets] = epoch_like_first_pass(EEG, RELAX_cfg)
% The first pass's boundary-aware epoching on a copy, with no channel deletion, no blink
% marking (the blink-aware shift bound is not used) and no muscle screen. With deletion off
% the epoching skips its screens, drift slopes included, so RELAX_excluding_extreme_values
% computes them, as it does whenever they are missing.
    copy = EEG;
    if ~isempty(copy.event)
        copy.event(cellfun(@(t) ischar(t) && startsWith(t, 'EyeBlink'), {copy.event.type})) = [];
    end
    copy = eeg_checkset(copy, 'eventconsistency');
    copy = rmfield(copy, intersect(fieldnames(copy), {'RELAX', 'RELAXProcessing', 'RELAXProcessingExtremeRejections'}));
    RELAX_cfg.ms_per_sample = 1000 / EEG.srate;
    RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0;
    RELAX_cfg.ProbabilityDataHasNoBlinks = 2;
    RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel = 1;
    [continuousCopy, epoched] = RELAX_excluding_channels_and_epoching_boundary(copy, RELAX_cfg);
    [~, epoched] = RELAX_excluding_extreme_values(continuousCopy, epoched, RELAX_cfg);
    assert(epoched.nbchan == EEG.nbchan, 'Second-pass epoching deleted channels');

    epochEvents = epoched.event(strcmp({epoched.event.type}, 'X'));
    assert(numel(epochEvents) == epoched.trials, 'Expected one epoch-onset event per epoch');
    onsets = zeros(1, epoched.trials);
    onsets([epochEvents.epoch]) = [epochEvents.originallatency];
end
