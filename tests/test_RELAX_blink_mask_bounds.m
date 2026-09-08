function test_RELAX_blink_mask_bounds()
% Blink mask windows, truncated to the segment containing the blink.

    segments = [1 4000; 4001 8000];
    nPoints  = 8000;
    hw       = 400; % half-width in samples

    % A blink well inside a segment is not truncated at all
    [s, e] = RELAX_blink_mask_bounds(2000, hw, segments, nPoints);
    assert(s == 1600 && e == 2400, 'interior blink must keep its full window');

    % A blink near the end of segment 1 is truncated on the right
    [s, e] = RELAX_blink_mask_bounds(3900, hw, segments, nPoints);
    assert(s == 3500 && e == 4000, 'window must be truncated at the segment end');

    % A blink near the start of segment 2 is truncated on the left
    [s, e] = RELAX_blink_mask_bounds(4100, hw, segments, nPoints);
    assert(s == 4001 && e == 4500, 'window must be truncated at the segment start');

    % Exactly on a segment start: belongs to that segment, truncated to it
    [s, e] = RELAX_blink_mask_bounds(4001, hw, segments, nPoints);
    assert(s == 4001 && e == 4401, 'blink on a segment start');

    % Exactly on a segment end: belongs to the earlier segment
    [s, e] = RELAX_blink_mask_bounds(4000, hw, segments, nPoints);
    assert(s == 3600 && e == 4000, 'blink on a segment end');

    % Truncation at the very start of the record, with one whole-record segment
    [s, e] = RELAX_blink_mask_bounds(50, hw, [1 nPoints], nPoints);
    assert(s == 1 && e == 450, 'window must be truncated at the record start');

    % Truncation at the very end of the record
    [s, e] = RELAX_blink_mask_bounds(nPoints-50, hw, [1 nPoints], nPoints);
    assert(s == nPoints-450 && e == nPoints, 'window must be truncated at the record end');

    % A latency outside every segment falls back to the whole record
    [s, e] = RELAX_blink_mask_bounds(99999, hw, segments, nPoints);
    assert(s == nPoints-hw && e == nPoints, 'out-of-range latency clamps to the record');

    % Degenerate case: a one-sample segment still yields a valid single-sample window
    [s, e] = RELAX_blink_mask_bounds(10, hw, [10 10], nPoints);
    assert(s == 10 && e == 10, 'one-sample segment yields a one-sample window');

    % Invariants that must hold for every case above
    testLatencies = [1 50 2000 3900 4000 4001 4100 7999 8000];
    for lat = testLatencies
        [s, e] = RELAX_blink_mask_bounds(lat, hw, segments, nPoints);
        assert(s <= e, 'maskStart must never exceed maskEnd (latency %d)', lat);
        assert(s >= 1 && e <= nPoints, 'window must stay inside the record (latency %d)', lat);
        assert(s == round(s) && e == round(e), 'bounds must be integers (latency %d)', lat);
        % The window must lie wholly within one segment
        inOne = any(segments(:,1) <= s & segments(:,2) >= e);
        assert(inOne, 'window must lie inside a single segment (latency %d)', lat);
    end

    % Non-integer latencies (the caller passes o+I from a max search) are rounded
    [s, e] = RELAX_blink_mask_bounds(2000.4, hw, segments, nPoints);
    assert(s == round(2000.4-hw) && e == round(2000.4+hw), 'non-integer latency handled');

    fprintf('test_RELAX_blink_mask_bounds: PASS\n');
end
