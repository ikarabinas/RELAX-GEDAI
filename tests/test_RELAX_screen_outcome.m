function test_RELAX_screen_outcome()
    names = ["MuscleProportion", "MuscleStrength", "BlinkRatio", "HoleDepth"];
    missingMatters = [true false true true];
    pipeline = ["cleaned"; "cleaned"; "cleaned"; "cleaned"; "excluded"; "cleaned"; "cleaned"; "cleaned"; "cleaned"];
    values = [0.1 -0.3 1.1 1.0; 0.9 -0.3 1.1 1.0; 0.1 -0.3 NaN 1.0; 0.9 -0.3 NaN 1.0; NaN NaN NaN NaN; 0.1 NaN 1.1 1.0; 0.9 -0.1 1.1 1.0; 0.1 -0.3 1.1 1.7; 0.1 -0.3 1.1 NaN];
    flagged = logical([0 0 0 0; 1 0 0 0; 0 0 0 0; 1 0 0 0; 1 1 1 1; 0 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 0 0]);
    [o, by] = RELAX_screen_outcome(pipeline, values, flagged, names, missingMatters);
    assert(isequal(o, ["pass"; "fail"; "unscreened"; "fail"; "notCleaned"; "pass"; "fail"; "fail"; "unscreened"]));
    assert(by(2) == "MuscleProportion" && by(4) == "MuscleProportion" && by(1) == "");
    assert(by(5) == "", 'a not-cleaned row is never flagged');
    assert(by(7) == "MuscleStrength");
    assert(by(8) == "HoleDepth", 'a row flagged only by the fourth metric fails');
    flagged(7, 1) = true;
    [~, by] = RELAX_screen_outcome(pipeline, values, flagged, names, missingMatters);
    assert(by(7) == "MuscleProportion; MuscleStrength");
    fprintf('test_RELAX_screen_outcome: PASS\n');
end
