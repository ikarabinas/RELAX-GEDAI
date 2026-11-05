function test_RELAX_excluding_extreme_values_min_bad_channels()
% RELAX_cfg.ExtremeMinBadChannelsToRejectEpoch = k rejects an epoch only when at least k
% channels trip a per-channel extreme criterion. On a 200-channel montage the original
% any-one-channel rule rejected most epochs on the strength of a single channel.
%
% Fixture: 300 s of noise on 8 channels with a 1-channel voltage burst at 100 s and a
% 2-channel burst at 200 s. Only the voltage-shift detector can fire.

    srate = 100;
    oneChannelBurst = 100 * srate;
    twoChannelBurst = 200 * srate;
    EEG = noise_fixture(srate, 8, 300);
    EEG.data(1, oneChannelBurst + (0:19)) = 1000;
    EEG.data([3 4], twoChannelBurst + (0:19)) = 1000;
    cfg = shift_only_cfg(srate);

    % Absent setting keeps the original rule: both bursts rejected.
    [defaultMarks, defaultOut] = run_extreme(EEG, cfg);
    assert(defaultMarks(oneChannelBurst) && defaultMarks(twoChannelBurst), ...
        'with the setting absent, a single bad channel no longer rejects its epoch');
    assert(defaultOut.RELAXProcessingExtremeRejections.ExtremeMinBadChannelsToRejectEpoch == 1, ...
        'the k in force must be recorded even when defaulted');

    cfg.ExtremeMinBadChannelsToRejectEpoch = 2;
    [marks, out] = run_extreme(EEG, cfg);
    assert(~marks(oneChannelBurst), 'an epoch with 1 bad channel was rejected under k = 2');
    assert(marks(twoChannelBurst), 'an epoch with 2 bad channels was not rejected under k = 2');

    rej = out.RELAXProcessingExtremeRejections;
    assert(rej.ExtremeMinBadChannelsToRejectEpoch == 2, 'k was not recorded');
    assert(abs(rej.ProportionExcludedIfAnyChannelRejected ...
        - defaultOut.RELAXProcessingExtremeRejections.ProportionExcludedForExtremeOutlier) < 1e-12, ...
        'the any-channel comparison proportion does not match what the original rule excludes');
    assert(rej.ProportionExcludedForExtremeOutlier < rej.ProportionExcludedIfAnyChannelRejected, ...
        'k = 2 did not exclude less than the any-channel rule');
    % 1 s epochs at a ~0.5 s step: the 0.2 s burst on channel 1 sits in 2 or 3 epochs
    % depending on where the grid falls, each kept with exactly that 1 flagged channel.
    assert(rej.EpochsKeptWithFlaggedChannels >= 2, ...
        'expected at least 2 kept epochs carrying flagged channels, got %d', rej.EpochsKeptWithFlaggedChannels);
    assert(rej.FlaggedChannelEpochsKept == rej.EpochsKeptWithFlaggedChannels, ...
        'each kept epoch should carry exactly 1 flagged channel (%d channel-epochs over %d epochs)', ...
        rej.FlaggedChannelEpochsKept, rej.EpochsKeptWithFlaggedChannels);

    % The whole-montage kurtosis test is not a per-channel criterion and still rejects
    % on its own, however large k is.
    spike = 150 * srate;
    EEG = noise_fixture(srate, 4, 300);
    EEG.data(:, spike) = 400;
    cfg = shift_only_cfg(srate);
    cfg.ExtremeVoltageShiftThreshold = 1e6;
    cfg.ExtremeAllChannelKurtosisThreshold = 10;
    cfg.ExtremeMinBadChannelsToRejectEpoch = 5;
    marks = run_extreme(EEG, cfg);
    assert(marks(spike), 'all-channel kurtosis stopped rejecting under a k above the channel count');

    fprintf('test_RELAX_excluding_extreme_values_min_bad_channels: PASS\n');
end

function EEG = noise_fixture(srate, nChan, seconds)
    rng(1);
    nSamples = seconds * srate;
    EEG = eeg_emptyset();
    EEG.data = 10 * randn(nChan, nSamples);
    EEG.srate = srate;
    EEG.nbchan = nChan;
    EEG.pnts = nSamples;
    EEG.trials = 1;
    EEG.xmin = 0;
    EEG.xmax = (nSamples-1)/srate;
    EEG.chanlocs = struct('labels', arrayfun(@(c) sprintf('E%d', c), 1:nChan, 'UniformOutput', false));
    EEG.allchan = EEG.chanlocs;
    EEG.event = struct('type',{},'latency',{});
end

function cfg = shift_only_cfg(srate)
    cfg = struct();
    cfg.ms_per_sample = 1000/srate;
    cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0;
    cfg.ProbabilityDataHasNoBlinks = 2;
    cfg.ExtremeVoltageShiftThreshold = 25;
    cfg.ExtremeAbsoluteVoltageThreshold = 1e6;
    cfg.ExtremeImprobableVoltageDistributionThreshold = 1e6;
    cfg.ExtremeSingleChannelKurtosisThreshold = 1e6;
    cfg.ExtremeAllChannelKurtosisThreshold = 1e6;
end

function [marks, continuousOut] = run_extreme(EEG, cfg)
    [continuousEEG, epochedEEG] = RELAX_excluding_channels_and_epoching_boundary(EEG, cfg);
    epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels = zeros(EEG.nbchan, epochedEEG.trials);
    [continuousOut, ~] = RELAX_excluding_extreme_values(continuousEEG, epochedEEG, cfg);
    marks = continuousOut.RELAXProcessing.Details.ExtremeOutlierPeriods == 1;
end
