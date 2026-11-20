function test_RELAX_QC_cleaned_on_extremes_timeline()
% Cleaned samples return at their pre-excision times, NaN across the deleted period.

    srate = 100;
    periods = [301 400];
    kept = [1:300, 401:1000];
    cleaned = struct('srate', srate, 'pnts', numel(kept), 'chanlocs', struct('labels', {'A', 'B'}), ...
                     'data', [kept; -kept]);
    seconds = [0.5 2.0 3.5 3.99 4.5];
    data = RELAX_QC_cleaned_on_extremes_timeline(cleaned, periods, 1000, {'B', 'A'}, seconds);

    assert(isequal(size(data), [2 5]));
    assert(abs(data(2, 1) - 51) < 1e-9 && abs(data(1, 1) + 51) < 1e-9, 'value at 0.5 s');
    assert(abs(data(2, 2) - 201) < 1e-9, 'value at 2.0 s');
    assert(all(isnan(data(:, 3))) && all(isnan(data(:, 4))), 'gap is NaN');
    assert(abs(data(2, 5) - 451) < 1e-9, 'value after the gap');

    threw = false;
    try, RELAX_QC_cleaned_on_extremes_timeline(cleaned, [301 399], 1000, {'A'}, seconds); catch, threw = true; end
    assert(threw, 'sample-count mismatch must error');
    fprintf('test_RELAX_QC_cleaned_on_extremes_timeline: PASS\n');
end
