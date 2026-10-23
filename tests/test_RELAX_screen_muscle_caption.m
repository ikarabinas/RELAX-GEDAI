function test_RELAX_screen_muscle_caption()
% The caption names the epoch, its start, its highest slope and where, and says when the slope is under the threshold.
    above = RELAX_screen_muscle_caption(1, 340.6, -0.12, "E90", -0.31);
    assert(above == "Muscle epoch 1, starts at 340.6 s: slope -0.12 on E90 (threshold -0.31)");
    below = RELAX_screen_muscle_caption(3, 12, -0.45, "E7", -0.31);
    assert(below == "Muscle epoch 3, starts at 12.0 s: slope -0.45 on E7 (threshold -0.31), below threshold");
    atThreshold = RELAX_screen_muscle_caption(2, 5.5, -0.31, "E7", -0.31);
    assert(endsWith(atThreshold, ", below threshold"), 'a slope equal to the threshold is not above it');
    fprintf('test_RELAX_screen_muscle_caption: PASS\n');
end
