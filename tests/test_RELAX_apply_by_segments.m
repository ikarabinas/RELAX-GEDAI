function test_RELAX_apply_by_segments()
% Per-segment application with reflection padding.

    data = reshape(1:20, 1, 20);
    segments = [1 10; 11 20];

    % Identity function must round-trip the data exactly
    out = RELAX_apply_by_segments(data, segments, @(d) d, 3);
    assert(isequal(out, data), 'identity fn must round-trip exactly');

    % A pointwise function must apply within every segment
    out = RELAX_apply_by_segments(data, segments, @(d) d*2, 3);
    assert(isequal(out, data*2), 'pointwise fn must apply across all segments');

    % Output size always matches input size
    assert(isequal(size(out), size(data)), 'output size must match input size');

    % ISOLATION: scaling segment 2 by 1e6 must not change segment 1's output at all.
    % A function that depends on the whole block it is given (mean subtraction) would
    % leak across segments if the implementation concatenated them.
    dataA = [ones(1,10), ones(1,10)];
    dataB = [ones(1,10), ones(1,10)*1e6];
    fn = @(d) d - mean(d,2);
    outA = RELAX_apply_by_segments(dataA, segments, fn, 3);
    outB = RELAX_apply_by_segments(dataB, segments, fn, 3);
    assert(isequal(outA(1:10), outB(1:10)), ...
        'segment 1 output must be bitwise unchanged when segment 2 changes');

    % Multi-channel data
    data2 = [reshape(1:20,1,20); reshape(21:40,1,20)];
    out2 = RELAX_apply_by_segments(data2, segments, @(d) d*3, 2);
    assert(isequal(out2, data2*3), 'multi-channel pointwise fn');

    % A single whole-record segment behaves like a plain call
    out3 = RELAX_apply_by_segments(data, [1 20], @(d) d*2, 3);
    assert(isequal(out3, data*2), 'single segment case');

    % I1: a segment shorter than an explicitly-passed minSegmentLength is passed
    % through UNCHANGED rather than handed to fn - a 3-sample island is too short for
    % a stable filter/wavelet transform. The DEFAULT minSegmentLength is 1 (fn applies
    % to every segment, see the assertions above), so this test opts into the skip
    % path explicitly, exactly as RELAX_filtbutter.m and RELAX_targeted_wICA.m do for
    % their filtfilt/wavelet calls. Suppress the expected warning so test output stays
    % pristine, and restore the warning state afterwards even on failure.
    warnState = warning('off', 'RELAX_apply_by_segments:segmentTooShort');
    cleanupWarn = onCleanup(@() warning(warnState)); %#ok<NASGU>

    dataShort = 1:200;
    shortSegments = [1 100; 101 103; 104 200];
    fnDestructive = @(d) -d; % would change every value it is applied to
    outShort = RELAX_apply_by_segments(dataShort, shortSegments, fnDestructive, 5, 16);

    assert(isequal(outShort(101:103), dataShort(101:103)), ...
        'the 3-sample segment must be passed through byte-identically');
    assert(isequal(outShort(1:100), -dataShort(1:100)), ...
        'the 100-sample segment must still be transformed');
    assert(isequal(outShort(104:200), -dataShort(104:200)), ...
        'the 97-sample segment must still be transformed');

    fprintf('test_RELAX_apply_by_segments: PASS\n');
end
