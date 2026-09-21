function test_RELAX_match_intervals()
% Interval overlap with slack, both directions, inclusive ends, empty-safe.
    aLo = [100; 500; 900]; aHi = [120; 520; 920];
    bLo = [110; 600; 930]; bHi = [115; 610; 940];
    [aHasB, bHasA] = RELAX_match_intervals(aLo, aHi, bLo, bHi, 0);
    assert(islogical(aHasB) && isequal(size(aHasB), [3 1]) && isequal(size(bHasA), [3 1]));
    assert(isequal(aHasB, [true; false; false]) && isequal(bHasA, [true; false; false]));

    [aHasB, bHasA] = RELAX_match_intervals(aLo, aHi, bLo, bHi, 10);
    assert(isequal(aHasB, [true; false; true]), 'slack 10 joins 920 and 930');
    assert(isequal(bHasA, [true; false; true]));

    [aHasB, bHasA] = RELAX_match_intervals(aLo, aHi, zeros(0,1), zeros(0,1), 10);
    assert(isequal(aHasB, false(3,1)) && isequal(size(bHasA), [0 1]));

    [aHasB, ~] = RELAX_match_intervals(zeros(0,1), zeros(0,1), bLo, bHi, 10);
    assert(isequal(size(aHasB), [0 1]));

    [aHasB, ~] = RELAX_match_intervals(100, 100, 100, 100, 0);
    assert(aHasB, 'touching single-sample intervals overlap');
    fprintf('test_RELAX_match_intervals: PASS\n');
end
