function flags = extreme_flags_per_channel(continuousEEG, epochedEEG, RELAX_cfg)
% Channel x epoch version of every criterion in RELAX_excluding_extreme_values.m.
% epochedEEG must be the output of RELAX_excluding_extreme_values, so the kurtosis,
% joint-probability and drift statistics it leaves behind are reused rather than
% recomputed. The voltage-shift bound is a local variable there, so it is rebuilt here
% line for line (lines 64-124 of that file).

    shift = range(epochedEEG.data, 2);              % chan x 1 x epoch
    upperBound = median(shift, 3) + RELAX_cfg.ExtremeVoltageShiftThreshold * mad(shift, 1, 3);
    if RELAX_cfg.ProbabilityDataHasNoBlinks < 2
        upperBound = max(blink_shift_bound(continuousEEG, epochedEEG, RELAX_cfg, shift), upperBound);
    end

    flags.shiftHigh = squeeze(shift > upperBound);
    flags.shiftLow  = squeeze(shift < 2);
    flags.absVolt   = squeeze(max(epochedEEG.data, [], 2) > RELAX_cfg.ExtremeAbsoluteVoltageThreshold) ...
                    | squeeze(min(epochedEEG.data, [], 2) < -RELAX_cfg.ExtremeAbsoluteVoltageThreshold);
    flags.kurtSingle = epochedEEG.stats.kurtE > RELAX_cfg.ExtremeSingleChannelKurtosisThreshold;
    flags.jointProb  = epochedEEG.stats.jpE > RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold;
    flags.drift      = epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels < RELAX_cfg.ExtremeDriftSlopeThreshold;
    % Whole-montage criterion, not attributable to any channel:
    flags.kurtAllEpoch = (epochedEEG.stats.kurt > RELAX_cfg.ExtremeAllChannelKurtosisThreshold)';
end

function bound = blink_shift_bound(continuousEEG, epochedEEG, RELAX_cfg, shift)
    bound = zeros(size(continuousEEG.data, 1), 1);
    [blinkMarked, ~] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, RELAX_cfg);
    if blinkMarked.RELAX.IQRmethodDetectedBlinks == 1
        eyeOnly = pop_selectevent(blinkMarked, 'type', {'EyeBlinkMax'}, 'deleteevents', 'on');
        eyeOnly = eeg_checkset(eyeOnly, 'eventconsistency');
        eyeOnly = pop_epoch(eyeOnly, {'EyeBlinkMax'}, [-0.5 0.5], 'epochinfo', 'yes');
        blinkShift = range(eyeOnly.data, 2);
        bound = median(blinkShift, 3) + RELAX_cfg.ExtremeBlinkShiftThreshold * mad(blinkShift, 1, 3);
    else
        largest = shift;
        largest(prctile(largest, 80, 3) > largest) = NaN;
        bound = median(largest, 3, 'omitnan') + RELAX_cfg.ExtremeBlinkShiftThreshold * mad(largest, 1, 3);
    end
end
