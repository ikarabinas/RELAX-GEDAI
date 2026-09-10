function test_RELAX_excluding_channels_and_epoching_boundary()
% RELAX_excluding_channels_and_epoching_boundary epochs each inter-boundary segment
% independently (via RELAX_get_segments) instead of on a fixed grid across the whole
% recording, so boundary-adjacent data is evaluated rather than left un-epoched (NaN)
% and later swept into ExtremelyBadPeriodsForDeletion regardless of whether it's extreme.
%
% Builds a synthetic 7-segment recording: six ~4s segments plus one 0.3s segment (too
% short for even one 1s epoch) sandwiched in the middle, and checks:
%   - no epoch straddles a boundary (pop_epoch's own discontinuity check stays a no-op)
%   - a segment too short for an epoch contributes none, and does not crash the
%     first/last-11-event edge relabelling
%   - a segment away from the file's edges is covered end-to-end with NO NaN gaps,
%     including its very last sample (regression test for the one-sample-per-segment
%     gap left by the epoch-tiling loop stopping short of end_latencies(k), which is
%     patched back in via `segmentHadEpoch`/segments(:,2) after the main mask-fill loop)
%
% RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted is set to 0 so the function takes
% its "no channels may be rejected" branch, skipping the kurtosis/jointprob/muscle-slope
% block (which needs FieldTrip and enough epochs/channels for stable statistics) -
% that block is unmodified from RELAX_excluding_channels_and_epoching_iTBS.m and is out
% of scope for this test.

    srate = 100;
    nChan = 4;
    segLens = [400 400 400 30 400 400 400]; % samples: six 4s segments, one 0.3s segment
    segEnds = cumsum(segLens);
    nSamples = segEnds(end);

    EEG = eeg_emptyset();
    EEG.data = randn(nChan, nSamples);
    EEG.srate = srate;
    EEG.nbchan = nChan;
    EEG.pnts = nSamples;
    EEG.trials = 1;
    EEG.xmin = 0;
    EEG.xmax = (nSamples-1)/srate;
    EEG.chanlocs = struct('labels', {'E1','E2','E3','E4'});
    EEG.allchan = EEG.chanlocs;
    EEG.event = struct('type',{},'latency',{});
    for b = 1:length(segEnds)-1
        EEG.event(end+1).type = 'boundary';
        EEG.event(end).latency = segEnds(b) + 0.5; % EEGLAB's between-sample convention
    end

    RELAX_cfg = struct();
    RELAX_cfg.ms_per_sample = 1000/srate;
    RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0;

    [continuousOut, epochedOut] = RELAX_excluding_channels_and_epoching_boundary(EEG, RELAX_cfg);

    segments = RELAX_get_segments(EEG);
    assert(size(segments,1) == 7, 'fixture should split into 7 segments, got %d', size(segments,1));
    shortSegmentIdx = 4;
    assert(isequal(segments(shortSegmentIdx,:), [1201 1230]), 'short segment should be samples 1201:1230');

    % No surviving epoch may contain a boundary event:
    boundaryInEpoch = 0;
    for ep = 1:numel(epochedOut.epoch)
        types = epochedOut.epoch(ep).eventtype;
        if ~iscell(types); types = {types}; end
        if any(strcmp(types, 'boundary'))
            boundaryInEpoch = boundaryInEpoch + 1;
        end
    end
    assert(boundaryInEpoch == 0, ...
        '%d epoch(s) contain a boundary event - pop_epoch should have dropped them', boundaryInEpoch);
    assert(size(epochedOut.data,3) == epochedOut.trials && epochedOut.trials > 0, ...
        'expected at least one surviving epoch');

    mask = continuousOut.RELAXProcessing.Details.NaNsForNonEvents;

    % The short segment is too short for even one epoch, so it contributes nothing and
    % stays entirely NaN - and must not have crashed the first/last-11 relabelling above:
    shortRange = segments(shortSegmentIdx,1):segments(shortSegmentIdx,2);
    assert(all(isnan(mask(shortRange))), 'short segment should be entirely un-evaluated (NaN)');

    % Segments 3 and 5 sit away from the file's edges, so none of their epochs are
    % relabelled 'Y' by the global first/last-11-event edge exclusion. Both must be
    % covered end-to-end, with NO gap at all - in particular not at the segment's very
    % last sample, which the epoch-tiling loop cannot itself reach (an epoch reaching
    % end_latencies(k) would span into the following boundary and get dropped by
    % pop_epoch) and which is patched back in separately:
    for segIdx = [3 5]
        segRange = segments(segIdx,1):segments(segIdx,2);
        nanCount = sum(isnan(mask(segRange)));
        assert(nanCount == 0, ...
            'segment %d should have zero NaN samples, found %d (last-sample patch may be missing)', ...
            segIdx, nanCount);
    end

    fprintf('test_RELAX_excluding_channels_and_epoching_boundary: PASS\n');
end
