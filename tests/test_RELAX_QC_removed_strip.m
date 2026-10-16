function test_RELAX_QC_removed_strip()
% Blink/muscle strips keep every channel at full range; the pop strip is the first channel only, with robust limits.
    rng(1);
    removed = randn(4, 1000);
    removed(2, 500) = 80;
    [traces, yRange, clipped] = RELAX_QC_removed_strip(removed, false);
    assert(isequal(traces, removed) && ~clipped, 'blink/muscle keep all channels, unclipped');
    assert(yRange(1) <= min(removed, [], 'all') && yRange(2) >= 80, 'blink/muscle range spans the data');

    removed(1, 300) = 5000;
    [traces, yRange, clipped] = RELAX_QC_removed_strip(removed, true);
    assert(isequal(traces, removed(1, :)), 'pop strip is the repaired channel only');
    assert(clipped && yRange(2) < 20 && yRange(1) == -yRange(2), 'pop limits are robust and symmetric');

    [~, yRange, clipped] = RELAX_QC_removed_strip(zeros(2, 100), true);
    assert(~clipped && yRange(2) >= 1, 'near-zero pop strip keeps a minimum range');
    fprintf('test_RELAX_QC_removed_strip: PASS\n');
end
