function test_RELAX_epoching_muscle_deletion_drift_rows()
% When the epoching stage deletes a muscle-contaminated channel, its drift slopes must go
% with it: RELAX_excluding_extreme_values combines the drift matrix channel by channel
% with the other per-channel criteria, so a leftover row breaks it (or, worse, shifts
% every later channel's drift onto its neighbour).
%
% Fixture: 60 s of pink noise on 8 channels; E8 is white noise, whose flat spectrum the
% muscle screen reads as muscle in every epoch.

    for variant = {@RELAX_excluding_channels_and_epoching_boundary, @RELAX_excluding_channels_and_epoching}
        [EEG, cfg] = fixture();
        [continuousEEG, epochedEEG] = variant{1}(EEG, cfg);
        name = func2str(variant{1});
        assert(any(strcmp(epochedEEG.RELAXProcessingExtremeRejections.MuscleBasedElectrodesToReject, 'E8')), ...
            '%s: the white-noise channel was not deleted for muscle', name);
        assert(size(epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels, 1) == epochedEEG.nbchan, ...
            '%s: drift slopes have %d rows for %d channels after muscle deletion', name, ...
            size(epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels, 1), epochedEEG.nbchan);
        RELAX_excluding_extreme_values(continuousEEG, epochedEEG, cfg);

        % Over budget: E7 and E8 are identical white noise but only 1 channel may go, so
        % the >= pass picks 2, and the fallback restores and reruns with > (picking none).
        [EEG, cfg] = fixture();
        EEG.data(7, :) = EEG.data(8, :);
        cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0.125;
        [continuousEEG, epochedEEG] = variant{1}(EEG, cfg);
        assert(size(epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels, 1) == epochedEEG.nbchan, ...
            '%s, over-budget fallback: drift slopes have %d rows for %d channels', name, ...
            size(epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels, 1), epochedEEG.nbchan);
        RELAX_excluding_extreme_values(continuousEEG, epochedEEG, cfg);
    end
    fprintf('test_RELAX_epoching_muscle_deletion_drift_rows: PASS\n');
end

function [EEG, cfg] = fixture()
    rng(1);
    srate = 250;
    nSamples = 60 * srate;
    data = 10 * make_pink_noise(8, nSamples);
    data(8, :) = 10 * randn(1, nSamples);
    EEG = eeg_emptyset();
    EEG.data = data;
    EEG.srate = srate; EEG.nbchan = 8; EEG.pnts = nSamples; EEG.trials = 1;
    EEG.xmin = 0; EEG.xmax = (nSamples - 1) / srate;
    EEG.chanlocs = make_qc_chanlocs(8);
    EEG.allchan = EEG.chanlocs;
    EEG.event = struct('type', {}, 'latency', {});

    cfg = struct();
    cfg.ms_per_sample = 1000 / srate;
    cfg.ProbabilityDataHasNoBlinks = 2;
    cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0.25;
    cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel = 0.5;
    cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel = 0.25;
    cfg.ExtremeMinBadChannelsToRejectEpoch = 1;
end
