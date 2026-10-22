function test_RELAX_screen_outcome()
    names = ["MuscleProportion", "MuscleStrength", "BlinkRatio", "NInterpolated"];
    missingMatters = [true false true true];
    pipeline = ["cleaned"; "cleaned"; "cleaned"; "cleaned"; "excluded"; "cleaned"; "cleaned"];
    values = [0.1 -0.3 1.1 3; 0.9 -0.3 1.1 3; 0.1 -0.3 NaN 3; 0.9 -0.3 NaN 3; NaN NaN NaN NaN; 0.1 NaN 1.1 3; 0.9 -0.1 1.1 3];
    z      = [0   0    0   0; 4   0    0   0; 0   0    NaN 0; 5   0    NaN 0; NaN NaN NaN NaN; 0   NaN 0   0; 3   3.5  0   0];
    [o, by] = RELAX_screen_outcome(pipeline, values, z, 3, names, missingMatters);
    assert(isequal(o, ["pass"; "fail"; "unscreened"; "fail"; "notCleaned"; "pass"; "fail"]));
    assert(by(2) == "MuscleProportion" && by(4) == "MuscleProportion" && by(1) == "");
    assert(by(7) == "MuscleStrength", 'z exactly at the cutoff does not flag');
    z(7, 1) = 3.2;
    [~, by] = RELAX_screen_outcome(pipeline, values, z, 3, names, missingMatters);
    assert(by(7) == "MuscleProportion; MuscleStrength");
    fprintf('test_RELAX_screen_outcome: PASS\n');
end
