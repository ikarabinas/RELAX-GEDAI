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

    fprintf('test_RELAX_apply_by_segments: PASS\n');
end
