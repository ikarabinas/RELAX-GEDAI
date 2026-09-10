function test_zapline_chunk_segmentation()
% Boundary-aligned covariance epochs for the vendored zapline-plus adaptive chunker.
% No epoch may straddle a boundary, and no inter-boundary segment short enough to
% be a single epoch may ever be subdivided.

    RELAX_add_vendored_zapline();

    srate  = 500;
    segLen = 3965;          % 7.93 s - the m208_dlpfc_day1 inter-train segment length
    nSeg   = 5;
    nSamp  = nSeg*segLen;
    bounds = (1:nSeg-1)'*segLen + 1;   % each boundary is the FIRST sample of a new segment

    % --- 1. segments tile the record exactly, one row per continuous stretch ---
    segments = zapline_boundary_segments(nSamp, bounds);
    assert(isequal(size(segments), [nSeg 2]), 'expected %d segments, got %d', nSeg, size(segments,1));
    assert(segments(1,1) == 1 && segments(end,2) == nSamp, 'segments must span the record');
    assert(all(segments(2:end,1) == segments(1:end-1,2) + 1), 'segments must be contiguous');
    assert(all(ismember(bounds, segments(:,1))), 'every boundary must start a segment');

    % --- 2. no boundaries at all degrades to a single whole-record segment ---
    assert(isequal(zapline_boundary_segments(nSamp, []), [1 nSamp]), ...
        'continuous data must produce exactly one segment');

    % --- 3. a segment below the whole-segment threshold gets exactly ONE epoch
    %        spanning it, and is never subdivided. 7.93 s < 15 s. ---
    epochSamples        = 1*srate;
    wholeSegmentSamples = 15*srate;
    minEpochSamples     = epochSamples;
    [epochs, epochSegment] = zapline_covariance_epochs(segments, ...
        epochSamples, wholeSegmentSamples, minEpochSamples);
    assert(size(epochs,1) == nSeg, 'each short segment must yield exactly one epoch (got %d)', size(epochs,1));
    assert(isequal(epochs, segments), 'a short segment''s epoch must span the whole segment');
    assert(isequal(epochSegment(:), (1:nSeg)'), 'epochSegment must map each epoch to its segment');

    % --- 4. no epoch may escape its segment, for a ragged layout ---
    raggedSamp = 100123;
    ragged = zapline_boundary_segments(raggedSamp, [3001; 3801; 40001; 40501]);
    [rEpochs, rSegIdx] = zapline_covariance_epochs(ragged, ...
        epochSamples, wholeSegmentSamples, minEpochSamples);
    for k = 1:size(rEpochs,1)
        seg = ragged(rSegIdx(k), :);
        assert(rEpochs(k,1) >= seg(1) && rEpochs(k,2) <= seg(2), ...
            'epoch %d [%d %d] escapes its segment [%d %d]', k, rEpochs(k,1), rEpochs(k,2), seg(1), seg(2));
    end

    % --- 5. a long segment is tiled from ITS OWN start, every epoch is at least
    %        segmentLength, and the last one ends exactly on the boundary ---
    longSeg = ragged(5,:);                       % 40501 .. 100123, 59623 samples
    inLong  = find(rSegIdx == 5);
    assert(rEpochs(inLong(1),1) == longSeg(1), 'tiling must start at the segment start');
    assert(rEpochs(inLong(end),2) == longSeg(2), 'the last epoch must end exactly on the boundary');
    lens = rEpochs(inLong,2) - rEpochs(inLong,1) + 1;
    assert(all(lens >= epochSamples), 'every epoch must be at least segmentLength (min was %d)', min(lens));
    assert(all(lens < 2*epochSamples), 'only the last epoch absorbs a remainder, so all stay under 2x');
    assert(lens(end) > epochSamples, 'this fixture was chosen to leave a real remainder');
    assert(isequal(rEpochs(inLong(2:end),1), rEpochs(inLong(1:end-1),2)+1), 'epochs must tile the segment');

    % --- 6. segments too short for a usable covariance emit no epoch at all,
    %        rather than contributing a rank-starved one that fakes a distance peak ---
    tiny = zapline_boundary_segments(nSamp, [201; 401]);   % two 200-sample slivers
    [tEpochs, tSegIdx] = zapline_covariance_epochs(tiny, ...
        epochSamples, wholeSegmentSamples, minEpochSamples);
    for k = 1:size(tEpochs,1)
        assert(tEpochs(k,2) - tEpochs(k,1) + 1 >= epochSamples, 'no epoch may be shorter than segmentLength');
    end
    assert(~any(ismember(tSegIdx, [1 2])), 'the two sub-500-sample segments must be dropped entirely');
end
