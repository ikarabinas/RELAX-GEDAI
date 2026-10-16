function test_RELAX_QC_pick_pops()
% popLargest is the largest own step; popKept the largest own step whose transition clears
% every later-excised period, never the same step as popLargest.

    labels = {'E1', 'E2', 'E3'};
    steps = struct( ...
        'label',       {{'E1'; 'E2'; 'E3'; 'E1'; 'E9'; 'E2'}}, ...
        'amplitude',   [900; 700; 600; -500; 4000; 300], ...
        'sample',      [1000; 5000; 9000; 13000; 15000; 17000], ...
        'onset',       [990; 4990; 8990; 12990; 14990; 16990], ...
        'offset',      [1010; 5010; 9010; 13010; 15010; 17010], ...
        'commonMode',  [false; false; false; false; false; true], ...
        'eyeMovement', [false; false; false; false; false; false]);
    inputSrate = 1000; extremesSrate = 250;

    picks = RELAX_QC_pick_pops(steps, labels, zeros(0, 2), inputSrate, extremesSrate);
    assert(picks.largest.label == "E1" && picks.largest.amplitude == 900, 'largest own step, E9 not in montage');
    assert(abs(picks.largest.centreSeconds - 0.999) < 1e-9 && all(abs(picks.largest.shadeSeconds - [0.989 1.009]) < 1e-9), 'times in seconds at the input rate');
    assert(picks.kept.label == "E2" && picks.kept.amplitude == 700, 'nothing excised: kept is the next-largest, never the largest');

    excisedAroundLargest = [1 + 0.9 * extremesSrate, 1 + 1.1 * extremesSrate];
    picks = RELAX_QC_pick_pops(steps, labels, excisedAroundLargest, inputSrate, extremesSrate);
    assert(picks.largest.amplitude == 900 && picks.kept.amplitude == 700, 'largest ignores excision; kept is the largest clear one');

    excisedAroundSecond = [1 + 4.95 * extremesSrate, 1 + 5.05 * extremesSrate];
    picks = RELAX_QC_pick_pops(steps, labels, excisedAroundSecond, inputSrate, extremesSrate);
    assert(picks.largest.amplitude == 900 && picks.kept.amplitude == 600, 'overlapping step skipped, with differing sample rates');

    touchingOnset = [1 + 8.9 * extremesSrate, 1 + 8.99 * extremesSrate];
    picks = RELAX_QC_pick_pops(steps, labels, [excisedAroundSecond; touchingOnset], inputSrate, extremesSrate);
    assert(picks.kept.amplitude == 500, 'the transition interval, not just the centre, must clear the excision');

    everything = [1, 1 + 20 * extremesSrate];
    picks = RELAX_QC_pick_pops(steps, labels, everything, inputSrate, extremesSrate);
    assert(~isempty(picks.largest) && isempty(picks.kept), 'no qualifying step leaves kept empty');

    onlyOne = structfun(@(f) f(1), steps, 'UniformOutput', false);
    picks = RELAX_QC_pick_pops(onlyOne, labels, zeros(0, 2), inputSrate, extremesSrate);
    assert(~isempty(picks.largest) && isempty(picks.kept), 'a single step cannot fill both columns');

    picks = RELAX_QC_pick_pops(table2struct(struct2table(steps), 'ToScalar', true), labels, [], inputSrate, extremesSrate);
    assert(picks.largest.amplitude == 900, 'a steps table works and empty periods mean none excised');

    picks = RELAX_QC_pick_pops([], labels, [], inputSrate, extremesSrate);
    assert(isempty(picks.largest) && isempty(picks.kept), 'no steps saved');

    fprintf('test_RELAX_QC_pick_pops: PASS\n');
end
