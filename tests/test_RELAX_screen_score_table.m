function test_RELAX_screen_score_table()
    n = 8;
    T = table("f" + (1:n)', repmat("p1", n, 1), [repmat("cleaned", 7, 1); "excluded"], ...
              'VariableNames', {'File', 'Participant', 'Outcome'});
    mp = [0.10 0.11 0.12 0.13 0.14 0.60 0.12 99]';               % f6 is the muscle outlier; f8 excluded, huge
    holeDepth = [1.5 1.4 1.0 1.0 1.0 1.0 1.0 99]';               % f1 is deeper than the limit, f2 exactly at it
    nInterpolated = [5 5 49 5 50 5 NaN 99]';                     % f5 reaches its limit of 50, f3 stops one short
    metrics = struct('MuscleProportion', num2cell(mp), 'MuscleStrength', num2cell(-0.3 * ones(n, 1)), ...
                     'BlinkRatio', num2cell([1.1 1.0 1.2 NaN 1.1 1.0 1.1 NaN]'), 'NInterpolated', num2cell(nInterpolated), ...
                     'HoleDepth', num2cell(holeDepth), 'HoleDeepestElectrode', num2cell(repmat("E7", n, 1)), ...
                     'HoleElectrodes', num2cell(repmat("E7 E8", n, 1)), 'Spacing', num2cell([0.02 0.02 0.02 0.02 0.02 0.03 0.02 99]'), ...
                     'NInterpolatedLimit', num2cell(50 * ones(n, 1)), 'NChannels', num2cell(203 * ones(n, 1)), ...
                     'MaxProportionDeleted', num2cell(0.25 * ones(n, 1)));
    opts = struct('cutoff', 3, 'sampleSize', 2, 'seed', 1, 'holeDepthLimit', 1.4);
    [s, meta] = RELAX_screen_score_table(T, metrics, opts);
    assert(height(s) == n && meta.nCleaned == 7);
    assert(s.MuscleProportion_median(1) == 0.12, 'excluded file does not enter the median');
    assert(s.ScreenOutcome(6) == "fail" && s.FlaggedBy(6) == "MuscleProportion");
    assert(s.ScreenOutcome(4) == "unscreened" && s.ScreenOutcome(8) == "notCleaned");
    assert(~ismember("NInterpolated_z", string(s.Properties.VariableNames)) && ~any(meta.zeroScaleMetrics == "NInterpolated"));
    assert(isequal(meta.unscreeningMetrics, ["MuscleProportion", "BlinkRatio", "HoleDepth"]));
    assert(isequal(meta.metricNames, ["MuscleProportion", "MuscleStrength", "BlinkRatio", "HoleDepth", "NInterpolated"]));
    assert(isequal(meta.zMetricNames, ["MuscleProportion", "MuscleStrength", "BlinkRatio"]));
    assert(meta.holeDepthLimit == 1.4);
    assert(meta.medianSpacing == 0.02, 'median spacing over cleaned files only');
    assert(sum(s.InReviewSample) == 2 && all(s.ScreenOutcome(s.InReviewSample) == "pass"));
    assert(s.MuscleProportion_flag(6) && ~s.MuscleProportion_flag(1), 'flag is true only for the outlier');
    assert(~s.MuscleProportion_flag(8), 'a non-cleaned row never flags');
    assert(s.ScreenOutcome(1) == "fail" && s.FlaggedBy(1) == "HoleDepth" && s.HoleDepth_flag(1), 'deeper than the limit fails');
    assert(s.ScreenOutcome(2) == "pass" && ~s.HoleDepth_flag(2), 'a depth exactly at the limit passes');
    assert(all(s.HoleDepth_limit == 1.4));
    assert(~s.HoleDepth_flag(8) && ~s.NInterpolated_flag(8), 'a non-cleaned row never flags');
    assert(s.ScreenOutcome(3) == "pass" && ~s.NInterpolated_flag(3), 'one channel short of the limit passes');
    assert(s.ScreenOutcome(5) == "fail" && s.FlaggedBy(5) == "NInterpolated" && s.NInterpolated_flag(5), 'reaching the limit fails');
    assert(s.NInterpolated_limit(5) == 50);
    assert(s.ScreenOutcome(7) == "pass" && ~s.NInterpolated_flag(7), 'a missing NInterpolated neither flags nor leaves the file unscreened');
    assert(s.HoleElectrodes(1) == "E7 E8" && s.HoleDeepestElectrode(1) == "E7");
    assert(all(s.NChannels(1:7) == 203) && all(s.MaxProportionDeleted(1:7) == 0.25));
    columns = string(s.Properties.VariableNames);
    assert(isequal(columns(1:3), ["File", "Participant", "PipelineOutcome"]));
    assert(isequal(columns(4:8), ["MuscleProportion", "MuscleProportion_median", "MuscleProportion_scale", "MuscleProportion_z", "MuscleProportion_flag"]));
    assert(isequal(columns(19:end), ["HoleDepth", "HoleDepth_limit", "HoleDepth_flag", "HoleDeepestElectrode", "HoleElectrodes", ...
        "NInterpolated", "NInterpolated_limit", "NInterpolated_flag", "NChannels", "MaxProportionDeleted", ...
        "ScreenOutcome", "FlaggedBy", "InReviewSample", "SheetStatus", "SheetMessage"]));
    rng(0); before = rand;
    rng(0); s2 = RELAX_screen_score_table(T, metrics, opts); after = rand;
    assert(isequal(s.InReviewSample, s2.InReviewSample) && before == after, 'seeded, and leaves global rng alone');
    opts.sampleSize = 50;
    s3 = RELAX_screen_score_table(T, metrics, opts);
    assert(isequal(s3.InReviewSample, s3.ScreenOutcome == "pass"));
    atCutoff = opts; atCutoff.cutoff = s.MuscleProportion_z(6);
    s4 = RELAX_screen_score_table(T, metrics, atCutoff);
    assert(~s4.MuscleProportion_flag(6), 'a z exactly at the cutoff does not flag');
    justBelowCutoff = opts; justBelowCutoff.cutoff = s.MuscleProportion_z(6) - 1e-9;
    s5 = RELAX_screen_score_table(T, metrics, justBelowCutoff);
    assert(s5.MuscleProportion_flag(6), 'a z just above the cutoff flags');
    fprintf('test_RELAX_screen_score_table: PASS\n');
end
