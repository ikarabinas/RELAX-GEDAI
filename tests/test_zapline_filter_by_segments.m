function test_zapline_filter_by_segments()
% The narrowband filter that drives zapline's adaptive chunk detection must not smear
% data across a discontinuity.

    RELAX_add_vendored_zapline();

    srate  = 500;
    segLen = 3965;
    nSamp  = 3*segLen;
    bounds = [segLen+1; 2*segLen+1];
    band   = [57.5 62.5];

    rng(0);
    data     = randn(nSamp, 4);
    segments = zapline_boundary_segments(nSamp, bounds);

    % bandpass() designs an order-202 linear-phase FIR here, independent of signal
    % length, and applies it single-pass with delay compensation. Reuse it as-is.
    [~, d] = bandpass(data(:,1), band, srate);
    p  = filtord(d);
    fn = @(x) bandpass(x, band, srate);

    out = zapline_filter_by_segments(data, segments, fn, p, p+1);
    assert(isequal(size(out), size(data)), 'output must keep the input size');

    % --- 1. Changing one segment must not change any other segment's output. ---
    perturbed = data;
    perturbed(1:segLen, :) = randn(segLen, 4);
    outP = zapline_filter_by_segments(perturbed, segments, fn, p, p+1);
    tail = segLen+1:nSamp;
    assert(max(abs(out(tail,:) - outP(tail,:)), [], 'all') < 1e-12, ...
        'segment 1 leaked across the boundary into segments 2-3');

    % --- 2. Each segment must equal filtering that segment on its own. ---
    for k = 1:size(segments,1)
        s = segments(k,1);
        e = segments(k,2);
        alone = zapline_filter_by_segments(data(s:e,:), [1 e-s+1], fn, p, p+1);
        assert(max(abs(out(s:e,:) - alone), [], 'all') < 1e-12, ...
            'segment %d differs from filtering it in isolation', k);
    end

    % --- 3. Regression witness: the naive whole-record call DOES leak, which is the
    %        behaviour this helper exists to replace. ---
    naive  = bandpass(data,      band, srate);
    naiveP = bandpass(perturbed, band, srate);
    leak = max(abs(naive(tail,:) - naiveP(tail,:)), [], 'all');
    assert(leak > 1e-6, 'expected the naive whole-record filter to leak (leak was %.3g)', leak);

    % --- 4. Continuous data must match a plain whole-record call away from the record
    %        edges, where the reflection padding has no influence (FIR support is p). ---
    cont = zapline_filter_by_segments(data, [1 nSamp], fn, p, p+1);
    core = 2*p+1 : nSamp-2*p;
    assert(max(abs(cont(core,:) - naive(core,:)), [], 'all') < 1e-10, ...
        'away from the record edges, one-segment filtering must match a plain call');

    % --- 5. A segment too short to filter is passed through unchanged, with a warning
    %        rather than a crash. Such segments emit no covariance epoch either, so the
    %        unfiltered data never reaches a covariance. ---
    shortSegs = zapline_boundary_segments(nSamp, [101; 201]);
    lastwarn('');
    % evalc swallows the expected warning text so the suite's output stays pristine,
    % while lastwarn still records that the warning fired.
    [~] = evalc('outShort = zapline_filter_by_segments(data, shortSegs, fn, p, p+1);');
    assert(isequal(outShort(1:200,:), data(1:200,:)), ...
        'segments below the filter minimum must be copied through unchanged');
    assert(~isempty(lastwarn), 'skipping a segment must warn');
end
