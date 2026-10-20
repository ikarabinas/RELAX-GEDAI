function test_RELAX_excluding_extreme_values_kurtosis()
% An epoch whose all-channel kurtosis exceeds ExtremeAllChannelKurtosisThreshold must be
% rejected itself. RELAX upstream indexed the rejection with the stale loop variable
% `epoch` left over from the absolute-voltage loop, so every all-channel kurtosis hit
% rejected the file's last epoch instead of the one that tripped it.
%
% One 300 s segment (enough epochs for a kurtosis z-score above 10) with a single sharp
% spike on every channel mid-recording. Every other detector's threshold is set out of
% reach so only kurtosis can mark anything.

    srate = 100;
    nChan = 4;
    nSamples = 300 * srate;
    spikeSample = 15000;

    rng(1);
    EEG = eeg_emptyset();
    EEG.data = 10 * randn(nChan, nSamples);
    EEG.data(:, spikeSample) = 400;
    EEG.srate = srate;
    EEG.nbchan = nChan;
    EEG.pnts = nSamples;
    EEG.trials = 1;
    EEG.xmin = 0;
    EEG.xmax = (nSamples-1)/srate;
    EEG.chanlocs = struct('labels', {'E1','E2','E3','E4'});
    EEG.allchan = EEG.chanlocs;
    EEG.event = struct('type',{},'latency',{});

    RELAX_cfg = struct();
    RELAX_cfg.ms_per_sample = 1000/srate;
    RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0;
    RELAX_cfg.ProbabilityDataHasNoBlinks = 2;
    RELAX_cfg.ExtremeVoltageShiftThreshold = 1e6;
    RELAX_cfg.ExtremeAbsoluteVoltageThreshold = 1e6;
    RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold = 1e6;
    RELAX_cfg.ExtremeSingleChannelKurtosisThreshold = 10;
    RELAX_cfg.ExtremeAllChannelKurtosisThreshold = 10;

    [continuousEEG, epochedEEG] = RELAX_excluding_channels_and_epoching_boundary(EEG, RELAX_cfg);
    epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels = zeros(nChan, epochedEEG.trials);
    [continuousOut, ~] = RELAX_excluding_extreme_values(continuousEEG, epochedEEG, RELAX_cfg);

    marked = continuousOut.RELAXProcessing.Details.ExtremeOutlierPeriods == 1;
    assert(marked(spikeSample), 'the epoch containing the kurtosis spike was not rejected');
    markedSamples = find(marked);
    assert(all(abs(markedSamples - spikeSample) < 2*srate), ...
        'samples far from the spike were rejected (first %d, last %d) - kurtosis marks landed on the wrong epoch', ...
        markedSamples(1), markedSamples(end));

    fprintf('test_RELAX_excluding_extreme_values_kurtosis: PASS\n');
end
