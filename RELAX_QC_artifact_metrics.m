% RELAX_QC_artifact_metrics() - Blink amplitude ratio and muscle metrics, computed by the
%                               pipeline's own metric functions, identically for any input.
%
% Usage:
%   >> m = RELAX_QC_artifact_metrics(EEG, cfg, blinkLabels);
%
% Used on both the pre-wICA baseline and the cleaned file so "before" and "after" are
% measured the same way (this config never computes raw muscle metrics). Muscle metrics
% use only RELAX.ListOfChannelsAfterRejections, as the wrapper does before its final
% interpolation, so the cleaned result reproduces CleanedMetrics exactly.
%
% Outputs:
%   m - .blinkRatio (mean over blinkLabels), .blinkRatioPerChannel (EEG.allchan order),
%       .muscleProportion, .muscleStrength

function m = RELAX_QC_artifact_metrics(EEG, cfg, blinkLabels)

    template = RELAX_QC_empty_summary();
    m = template.metricsCleaned;

    % The metric functions recompute from scratch, and write to .Cleaned, only when this is set:
    EEG.RELAX.Data_has_been_cleaned = 1;
    EEG = removeFields(EEG, {'RELAXProcessing', 'RELAX_Metrics'});

    [continuousEEG, epochedEEG] = RELAX_epoching(EEG, cfg);
    continuousEEG = RELAX_metrics_blinks(continuousEEG, epochedEEG);
    if isfield(continuousEEG, 'RELAX_Metrics') && isfield(continuousEEG.RELAX_Metrics.Cleaned, 'BlinkAmplitudeRatio')
        % RELAX_metrics_blinks interpolates to EEG.allchan before measuring, so rows follow that montage:
        m.blinkRatioPerChannel = continuousEEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio;
        m.blinkRatio = mean(m.blinkRatioPerChannel(ismember({EEG.allchan.labels}, blinkLabels)), 'omitnan');
    end

    muscleEEG = pop_select(EEG, 'channel', EEG.RELAX.ListOfChannelsAfterRejections);
    [continuousEEG, epochedEEG] = RELAX_epoching(muscleEEG, cfg);
    continuousEEG = RELAX_metrics_muscle(continuousEEG, epochedEEG, cfg);
    if isfield(continuousEEG, 'RELAX_Metrics') && isfield(continuousEEG.RELAX_Metrics.Cleaned, 'ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel')
        m.muscleProportion = continuousEEG.RELAX_Metrics.Cleaned.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel;
        m.muscleStrength   = continuousEEG.RELAX_Metrics.Cleaned.MeanMuscleStrengthFromOnlySuperThresholdValues;
    end
end

function EEG = removeFields(EEG, names)
    present = names(isfield(EEG, names));
    if ~isempty(present), EEG = rmfield(EEG, present); end
end
