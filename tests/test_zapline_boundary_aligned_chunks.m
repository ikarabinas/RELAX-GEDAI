function test_zapline_boundary_aligned_chunks()
% End-to-end: on data whose continuous stretches are all shorter than the
% whole-segment threshold, every zapline chunk must be a whole number of complete
% inter-boundary segments. Fixture mirrors m208_dlpfc_day1's 7.93 s inter-train
% segments at half the sampling rate.

    RELAX_add_vendored_zapline();

    srate  = 250;
    segLen = 1983;                 % 7.93 s
    nSeg   = 30;
    nChan  = 8;
    nSamp  = nSeg*segLen;
    bounds = (1:nSeg-1)'*segLen + 1;

    rng(0);
    t    = (0:nSamp-1)'/srate;
    data = randn(nSamp, nChan) + 8*sin(2*pi*60*t)*ones(1, nChan);

    cfg = struct('plotResults', 0, 'noisefreqs', 60, 'boundaries', bounds);
    [cleaned, ~, analytics] = clean_data_with_zapline_plus(data, srate, cfg);

    % --- 1. The chunk layout is now recorded, per noise frequency. ---
    assert(isfield(analytics, 'chunkIndices'), 'analyticsResults must record chunkIndices');
    assert(iscell(analytics.chunkIndices) && ~isempty(analytics.chunkIndices), ...
        'chunkIndices must be one cell per noise frequency');
    for iFreq = 1:numel(analytics.chunkIndices)
        ci = analytics.chunkIndices{iFreq};
        assert(ci(1) == 1 && ci(end) == nSamp+1, ...
            'noise frequency %d: chunk indices must span the record', iFreq);
        assert(all(diff(ci) > 0), 'noise frequency %d: chunk indices must increase', iFreq);
        % Without this, assertions 2 and 4 below pass vacuously on a degenerate
        % single-chunk result (ci = [1; nSamp+1]): they are only checking properties
        % of interior edges, and a collapse to one whole-record chunk has none. This
        % is exactly the failure mode the numel(distances)>=3 / non-empty-proms guards
        % added by this task can produce, so it needs its own explicit gate.
        assert(numel(ci) >= 3, ...
            'noise frequency %d: expected at least one interior chunk edge, got a single whole-record chunk', iFreq);
    end

    % --- 2. THE INVARIANT: no chunk edge falls inside a continuous segment. ---
    for iFreq = 1:numel(analytics.chunkIndices)
        ci = analytics.chunkIndices{iFreq};
        stray = setdiff(ci(2:end-1), bounds);
        assert(isempty(stray), 'noise frequency %d: chunk edges cut inside continuous segments: %s', ...
            iFreq, mat2str(stray'));
    end

    % --- 3. The cleaning still works: 60 Hz must be strongly attenuated. ---
    [pRaw,   f] = pwelch_segmented(data,    [], [], 8192, srate, bounds);
    [pClean, ~] = pwelch_segmented(cleaned, [], [], 8192, srate, bounds);
    noiseBin = f > 59.5 & f < 60.5;
    attenDb  = 10*log10(mean(pRaw(noiseBin,:), 'all') / mean(pClean(noiseBin,:), 'all'));
    assert(attenDb > 10, 'expected >10 dB attenuation at 60 Hz, got %.2f dB', attenDb);

    % --- 4. Every chunk's shortest continuous stretch is a FULL segment, so the
    %        pwelch window inside that chunk can never collapse onto a fragment. ---
    for iFreq = 1:numel(analytics.chunkIndices)
        ci = analytics.chunkIndices{iFreq};
        for k = 1:numel(ci)-1
            chunkStart = ci(k);
            chunkEnd   = ci(k+1) - 1;
            inChunk    = bounds(bounds > chunkStart & bounds <= chunkEnd);
            pieces     = diff([chunkStart; inChunk; chunkEnd+1]);
            assert(all(pieces == segLen), ...
                'noise frequency %d, chunk %d contains a partial segment (piece lengths %s)', ...
                iFreq, k, mat2str(pieces'));
        end
    end

    assert(all(isfinite(analytics.proportionRemovedNoise)), ...
        'proportionRemovedNoise must be finite');

    % --- 5. SEGMENTS LONGER THAN THE WHOLE-SEGMENT THRESHOLD: exercise the
    %        "leaves room" half of the edge-legality rule (isSegmentStart ||
    %        leavesRoom in zapline_chunk_indices). The short-segment fixture above
    %        can never reach that branch: every segment there is shorter than
    %        wholeSegmentSamples, so it collapses to a single whole-segment epoch and
    %        every candidate edge is trivially a segment start. Here every segment is
    %        far longer than W, so an edge can legally land inside a segment as long
    %        as it leaves at least W samples of that segment on both sides. Uses
    %        clean_data_with_zapline_plus's defaults (minChunkLength=30s,
    %        segmentLength=1s), not overridden here, matching the fixture above -
    %        so W = max(30/2*250, 1*250) = 3750.
    W          = 3750;
    longSegLen = 20000;                 % 80 s, >> 2*W so an interior edge can exist
    nSeg2      = 2;
    nSamp2     = nSeg2*longSegLen;
    bounds2    = (1:nSeg2-1)'*longSegLen + 1;

    rng(1);
    t2   = (0:nSamp2-1)'/srate;
    data2 = randn(nSamp2, nChan);

    % Engineer a hard covariance jump exactly halfway through segment 1 (10000
    % samples from each end, comfortably >= W), by switching which channels carry
    % the big 60 Hz component partway through the segment. This is what should
    % surface as a dominant peak in the covariance-distance series and give
    % zapline_chunk_indices a genuine candidate edge that is NOT a segment start -
    % the case that exercises `leavesRoom`. The transition sample is chosen to land
    % exactly on an epoch boundary (10000 samples = 40 epochs of 250 samples each).
    half   = longSegLen/2;
    gainsA = (1:nChan);
    gainsB = (nChan:-1:1);
    sig60  = 20*sin(2*pi*60*t2);
    data2(1:half, :)            = data2(1:half, :)            + sig60(1:half, :)            .* gainsA;
    data2(half+1:longSegLen, :) = data2(half+1:longSegLen, :) + sig60(half+1:longSegLen, :) .* gainsB;
    % Segment 2: uniform gains throughout, no internal jump - just a plain segment.
    data2(longSegLen+1:end, :) = data2(longSegLen+1:end, :) + sig60(longSegLen+1:end, :) .* gainsA;

    cfg2 = struct('plotResults', 0, 'noisefreqs', 60, 'boundaries', bounds2);
    [~, ~, analytics2] = clean_data_with_zapline_plus(data2, srate, cfg2);

    assert(isfield(analytics2, 'chunkIndices') && iscell(analytics2.chunkIndices) && ~isempty(analytics2.chunkIndices), ...
        'long-segment run: analyticsResults must record chunkIndices');

    segments2 = zapline_boundary_segments(nSamp2, bounds2);

    for iFreq = 1:numel(analytics2.chunkIndices)
        ci2 = analytics2.chunkIndices{iFreq};
        assert(ci2(1) == 1 && ci2(end) == nSamp2+1, ...
            'long-segment run, noise frequency %d: chunk indices must span the record', iFreq);
        % Same non-vacuous gate as assertion 1 above, applied to this run.
        assert(numel(ci2) >= 3, ...
            'long-segment run, noise frequency %d: expected at least one interior chunk edge, got a single whole-record chunk', iFreq);

        interiorEdges = ci2(2:end-1);
        for iEdge = 1:numel(interiorEdges)
            edge = interiorEdges(iEdge);
            isBoundary = any(edge == bounds2);
            segIdx = find(segments2(:,1) <= edge & edge <= segments2(:,2), 1);
            assert(~isempty(segIdx), ...
                'long-segment run, noise frequency %d: edge %d is not inside any segment', iFreq, edge);
            seg = segments2(segIdx, :);
            leavesRoom = (edge - seg(1) >= W) && (seg(2) - edge + 1 >= W);
            assert(isBoundary || leavesRoom, ...
                ['long-segment run, noise frequency %d: edge %d is neither a boundary nor leaves >= ' ...
                 '%d samples of margin in segment [%d %d]'], ...
                iFreq, edge, W, seg(1), seg(2));
        end
    end
end
