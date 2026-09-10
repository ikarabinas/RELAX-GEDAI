function test_zapline_chunk_indices()
% Chunk edges for the vendored zapline-plus adaptive chunker must never cut a stretch
% of continuous data into a fragment shorter than the edge margin, and minChunkLength
% must be enforced in samples rather than epoch indices.

    RELAX_add_vendored_zapline();

    srate  = 500;
    segLen = 3965;          % 7.93 s
    nSeg   = 40;
    nSamp  = nSeg*segLen;
    bounds = (1:nSeg-1)'*segLen + 1;

    epochSamples        = 1*srate;
    wholeSegmentSamples = 15*srate;
    minChunkSamples     = 30*srate;

    segments = zapline_boundary_segments(nSamp, bounds);
    [epochs, epochSegment] = zapline_covariance_epochs(segments, ...
        epochSamples, wholeSegmentSamples, epochSamples);

    % Propose an edge after every epoch, so nothing but the legality and spacing
    % rules can be what constrains the result.
    locs = (1:size(epochs,1)-1)';
    pks  = ones(size(locs));
    chunkIndices = zapline_chunk_indices(nSamp, epochs, epochSegment, segments, ...
        locs, pks, minChunkSamples, wholeSegmentSamples);

    % --- 1. shape ---
    assert(chunkIndices(1) == 1 && chunkIndices(end) == nSamp+1, ...
        'chunk indices must span the record');
    assert(all(diff(chunkIndices) > 0), 'chunk indices must be strictly increasing');

    % --- 2. On 7.93 s segments every epoch IS a whole segment, so every chunk edge
    %        must land exactly on a boundary. This is the invariant. ---
    interior = chunkIndices(2:end-1);
    assert(all(ismember(interior, bounds)), ...
        'chunk edges cut inside continuous segments: %s', mat2str(setdiff(interior, bounds)'));

    % --- 3. minChunkLength is enforced in SAMPLES: every chunk is at least 30 s even
    %        though epochs are 7.93 s each. ---
    chunkLens = diff(chunkIndices);
    assert(all(chunkLens >= minChunkSamples), ...
        'shortest chunk was %d samples, below the %d-sample minimum', min(chunkLens), minChunkSamples);

    % --- 4. ...and it must not degenerate. Interpreting minChunkLength as 30 EPOCHS
    %        would be 238 s here and would permit at most 2 chunks. ---
    assert(numel(chunkIndices)-1 >= 5, ...
        'expected many chunks from sample-domain spacing, got %d', numel(chunkIndices)-1);

    % --- 5. Inside a long segment an edge is legal only away from both its edges. ---
    longSamp = 100000;
    longSegs = zapline_boundary_segments(longSamp, 5001);
    [lEpochs, lSegIdx] = zapline_covariance_epochs(longSegs, ...
        epochSamples, wholeSegmentSamples, epochSamples);
    lLocs   = (1:size(lEpochs,1)-1)';
    lChunks = zapline_chunk_indices(longSamp, lEpochs, lSegIdx, longSegs, ...
        lLocs, ones(size(lLocs)), minChunkSamples, wholeSegmentSamples);
    for e = lChunks(2:end-1)'
        seg = longSegs(find(longSegs(:,1) <= e & longSegs(:,2) >= e, 1), :);
        assert(e == seg(1) || (e - seg(1) >= wholeSegmentSamples && seg(2) - e + 1 >= wholeSegmentSamples), ...
            'edge %d cuts segment [%d %d] closer than %d samples to an edge', ...
            e, seg(1), seg(2), wholeSegmentSamples);
    end
    assert(any(lChunks(2:end-1) > 5001), 'a long segment must still be splittable in its interior');

    % --- 6. Taller peaks win ties the way findpeaks'' MinPeakDistance does. ---
    % Both candidates sit far enough into the record that the first-chunk guard
    % cannot delete the survivor, so this isolates the spacing rule.
    twoLocs = [5; 6];
    twoPks  = [1; 5];                       % the second candidate is much taller
    twoChunks = zapline_chunk_indices(nSamp, epochs, epochSegment, segments, ...
        twoLocs, twoPks, minChunkSamples, wholeSegmentSamples);
    assert(numel(twoChunks) == 3, 'the two candidates are 3965 samples apart, so only one survives');
    assert(twoChunks(2) == epochs(7,1), 'the taller peak''s edge must be the survivor');

    % --- 7. No proposed peaks at all gives a single whole-record chunk. ---
    none = zapline_chunk_indices(nSamp, epochs, epochSegment, segments, [], [], ...
        minChunkSamples, wholeSegmentSamples);
    assert(isequal(none, [1; nSamp+1]), 'with no peaks the record must be one chunk');
end
