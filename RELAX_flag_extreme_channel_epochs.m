function flags = RELAX_flag_extreme_channel_epochs(scores, driftSlopes, onsetSeconds, RELAX_cfg)
% Flags the channel-epochs of cleaned data that are extreme by the first pass's standard
% (scores from RELAX_score_locked_extremes). Joint probability and kurtosis flag only the
% extreme side. Each criterion is a channels x epochs logical; Any is their OR.
%
% An epoch with more than RELAX_cfg.SecondPassMaxProportionFlaggedChannels (absent = 0.20)
% of its channels flagged should not have survived the first pass, so it is an error
% (RELAX:secondPassTooManyFlaggedChannels) rather than something to interpolate.

    flags.ShiftHigh = scores.range > scores.shiftUpper;
    flags.ShiftLow = scores.range < 2;
    flags.AbsVolt = scores.maxAbs > RELAX_cfg.ExtremeAbsoluteVoltageThreshold;
    flags.JointProb = scores.jpZ > RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold;
    flags.Drift = driftSlopes < RELAX_cfg.ExtremeDriftSlopeThreshold;
    flags.Kurtosis = scores.kurtZ > RELAX_cfg.ExtremeSingleChannelKurtosisThreshold;
    flags.Any = flags.ShiftHigh | flags.ShiftLow | flags.AbsVolt | flags.JointProb | flags.Drift | flags.Kurtosis;

    if ~isfield(RELAX_cfg, 'SecondPassMaxProportionFlaggedChannels')
        RELAX_cfg.SecondPassMaxProportionFlaggedChannels = 0.20;
    end
    nChan = size(flags.Any, 1);
    flaggedPerEpoch = sum(flags.Any, 1);
    tooMany = find(flaggedPerEpoch > RELAX_cfg.SecondPassMaxProportionFlaggedChannels * nChan, 1);
    if ~isempty(tooMany)
        error('RELAX:secondPassTooManyFlaggedChannels', ...
            'Epoch at %.1f s has %d of %d channels flagged after cleaning (limit %.0f%%); the first pass should have removed it', ...
            onsetSeconds(tooMany), flaggedPerEpoch(tooMany), nChan, 100 * RELAX_cfg.SecondPassMaxProportionFlaggedChannels);
    end
end
