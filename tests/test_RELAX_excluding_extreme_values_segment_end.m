function test_RELAX_excluding_extreme_values_segment_end()
% RELAX_excluding_channels_and_epoching_boundary stops each segment's last epoch one
% sample short of the segment end (an epoch reaching it would be dropped by pop_epoch)
% and stamps that final sample clean on the last epoch's behalf. If the last epoch is
% then rejected as extreme, the final sample must be rejected with it; otherwise
% eeg_eegrej leaves a one-sample island between the deletion and the boundary, which
% GEDAI (and anything else treating segments as signals) mishandles.
%
% Seven 4 s segments so the middle segment's epochs are all 'X' (not relabelled 'Y' by
% the first/last-11-event edge exclusion). A spike in the middle segment's final 0.1 s
% lies only in its last epoch, so only that epoch is extreme.
%
% DriftRatioEpochsxChannels is pre-filled with flat slopes so the drift check, which
% would otherwise need FieldTrip, marks nothing.

    srate = 100;
    nChan = 4;
    segLen = 400;
    nSeg = 7;
    nSamples = segLen * nSeg;
    spikeSegment = 4;

    rng(1);
    EEG = eeg_emptyset();
    EEG.data = 10 * randn(nChan, nSamples);
    EEG.srate = srate;
    EEG.nbchan = nChan;
    EEG.pnts = nSamples;
    EEG.trials = 1;
    EEG.xmin = 0;
    EEG.xmax = (nSamples-1)/srate;
    EEG.chanlocs = struct('labels', {'E1','E2','E3','E4'});
    EEG.allchan = EEG.chanlocs;
    EEG.event = struct('type',{},'latency',{});
    for b = 1:nSeg-1
        EEG.event(end+1).type = 'boundary';
        EEG.event(end).latency = b*segLen + 0.5;
    end

    segments = RELAX_get_segments(EEG);
    spikeSegmentEnd = segments(spikeSegment,2);
    EEG.data(1, spikeSegmentEnd-5) = 5000;

    RELAX_cfg = struct();
    RELAX_cfg.ms_per_sample = 1000/srate;
    RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0;
    RELAX_cfg.ProbabilityDataHasNoBlinks = 2;

    [continuousEEG, epochedEEG] = RELAX_excluding_channels_and_epoching_boundary(EEG, RELAX_cfg);
    epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels = zeros(nChan, epochedEEG.trials);
    [continuousOut, ~] = RELAX_excluding_extreme_values(continuousEEG, epochedEEG, RELAX_cfg);

    periods = continuousOut.RELAX.ExtremelyBadPeriodsForDeletion;
    spikeRow = find(periods(:,1) <= spikeSegmentEnd-5 & periods(:,2) >= spikeSegmentEnd-5);
    assert(numel(spikeRow) == 1, 'the spike should fall inside exactly one deletion period');
    assert(periods(spikeRow,2) == spikeSegmentEnd, ...
        'deletion around the spike ends at sample %d, should reach the segment end %d (one-sample island left behind)', ...
        periods(spikeRow,2), spikeSegmentEnd);
    assert(continuousOut.RELAXProcessing.Details.ExtremeOutlierPeriods(spikeSegmentEnd) == 1, ...
        'segment-end sample should be counted in ExtremeOutlierPeriods');
    assert(isnan(continuousOut.RELAXProcessing.Details.NoiseMaskFullLength(spikeSegmentEnd)), ...
        'segment-end sample should be excluded from the MWF noise mask');

    % Segment ends whose last epoch was clean must stay clean:
    otherEnd = segments(spikeSegment-1,2);
    assert(continuousOut.RELAXProcessing.Details.ExtremeOutlierPeriods(otherEnd) == 0, ...
        'segment-end sample of a clean segment should not be marked extreme');

    fprintf('test_RELAX_excluding_extreme_values_segment_end: PASS\n');
end
