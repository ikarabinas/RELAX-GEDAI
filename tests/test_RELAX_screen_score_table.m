function test_RELAX_screen_score_table()
    n = 8;
    T = table("f" + (1:n)', repmat("p1", n, 1), [repmat("cleaned", 7, 1); "excluded"], ...
              'VariableNames', {'File', 'Participant', 'Outcome'});
    mp = [0.10 0.11 0.12 0.13 0.14 0.60 0.12 99]';               % f6 is the muscle outlier; f8 excluded, huge
    metrics = struct('MuscleProportion', num2cell(mp), 'MuscleStrength', num2cell(-0.3 * ones(n, 1)), ...
                     'BlinkRatio', num2cell([1.1 1.0 1.2 NaN 1.1 1.0 1.1 NaN]'), 'NInterpolated', num2cell(5 * ones(n, 1)));
    opts = struct('cutoff', 3, 'sampleSize', 2, 'seed', 1);
    [s, meta] = RELAX_screen_score_table(T, metrics, opts);
    assert(height(s) == n && meta.nCleaned == 7);
    assert(s.MuscleProportion_median(1) == 0.12, 'excluded file does not enter the median');
    assert(s.ScreenOutcome(6) == "fail" && s.FlaggedBy(6) == "MuscleProportion");
    assert(s.ScreenOutcome(4) == "unscreened" && s.ScreenOutcome(8) == "notCleaned");
    assert(all(isnan(s.NInterpolated_z)) && any(meta.zeroScaleMetrics == "NInterpolated"));
    assert(sum(s.InReviewSample) == 2 && all(s.ScreenOutcome(s.InReviewSample) == "pass"));
    rng(0); before = rand;
    rng(0); s2 = RELAX_screen_score_table(T, metrics, opts); after = rand;
    assert(isequal(s.InReviewSample, s2.InReviewSample) && before == after, 'seeded, and leaves global rng alone');
    opts.sampleSize = 50;
    s3 = RELAX_screen_score_table(T, metrics, opts);
    assert(isequal(s3.InReviewSample, s3.ScreenOutcome == "pass"));
    fprintf('test_RELAX_screen_score_table: PASS\n');
end
