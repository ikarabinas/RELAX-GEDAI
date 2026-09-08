function test_RELAX_get_segments()
% Segment derivation from boundary events.

    % No event field at all -> one whole-record segment
    EEG = struct('pnts', 200);
    assert(isequal(RELAX_get_segments(EEG), [1 200]), 'no event field should give one segment');

    % Empty event struct -> one whole-record segment
    EEG = struct('pnts', 200, 'event', struct('type', {}, 'latency', {}));
    assert(isequal(RELAX_get_segments(EEG), [1 200]), 'empty events should give one segment');

    % Events present but none are boundaries -> one whole-record segment
    EEG = struct('pnts', 200, 'event', struct('type', {'X','Y'}, 'latency', {50, 120}));
    assert(isequal(RELAX_get_segments(EEG), [1 200]), 'non-boundary events should give one segment');

    % Half-integer boundary (EEGLAB convention) splits without duplicating a sample
    EEG = struct('pnts', 200, 'event', struct('type', {'boundary'}, 'latency', {100.5}));
    assert(isequal(RELAX_get_segments(EEG), [1 100; 101 200]), 'half-integer boundary split');

    % Integer boundary splits without duplicating a sample
    EEG = struct('pnts', 200, 'event', struct('type', {'boundary'}, 'latency', {100}));
    assert(isequal(RELAX_get_segments(EEG), [1 100; 101 200]), 'integer boundary split');

    % Boundaries out of order in the event array are sorted
    EEG = struct('pnts', 300, 'event', struct('type', {'boundary','boundary'}, 'latency', {200.5, 100.5}));
    assert(isequal(RELAX_get_segments(EEG), [1 100; 101 200; 201 300]), 'boundaries should be sorted');

    % Numeric event types must not crash the type test
    EEG = struct('pnts', 200, 'event', struct('type', {1, 'boundary'}, 'latency', {10, 100.5}));
    assert(isequal(RELAX_get_segments(EEG), [1 100; 101 200]), 'numeric event types tolerated');

    % Segments always tile the record with no gaps and no overlaps
    EEG = struct('pnts', 1000, 'event', struct('type', {'boundary','boundary','boundary'}, ...
                                               'latency', {250.5, 500.5, 750.5}));
    segs = RELAX_get_segments(EEG);
    assert(segs(1,1) == 1, 'first segment starts at 1');
    assert(segs(end,2) == 1000, 'last segment ends at pnts');
    assert(all(segs(2:end,1) == segs(1:end-1,2) + 1), 'segments must tile without gap or overlap');

    fprintf('test_RELAX_get_segments: PASS\n');
end
