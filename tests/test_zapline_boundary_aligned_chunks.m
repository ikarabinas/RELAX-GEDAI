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
end
