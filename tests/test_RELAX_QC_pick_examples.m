function test_RELAX_QC_pick_examples()
% Automatic examples: the most typical file plus the worst by blink, muscle and pops,
% with a file filling several roles shown once. Overrides replace the choice.

    files = ["a"; "b"; "c"; "d"; "e"];
    scores = struct('blinkAmplitude', [1; 5; 2; 3; NaN], 'musclePower', [2; 1; 5; 3; 4], 'popSteps', [1; 2; 3; 9; 4]);
    % mean ranks a 1.33, b 2.33, c 3.33, d 3.67, e 4 -> median 3.33 is c
    picks = RELAX_QC_pick_examples(files, scores, []);
    assert(isequal([picks.file], ["c", "b", "d"]), 'got %s', strjoin([picks.file], ','));
    assert(picks(1).role == "median, worst muscle", 'roles merge: %s', picks(1).role);
    assert(picks(2).role == "worst blink" && picks(3).role == "most pops", 'single roles');
    assert(all(isnan([picks(1).centres.blink, picks(1).centres.muscle, picks(1).centres.pop])), 'automatic centres');

    overrides = struct('file', {"b", "b", "a", "zz"}, 'kind', {"blink", "pop", "muscle", "blink"}, ...
                       'centreSeconds', {12.5, 30, NaN, 1});
    lastwarn('');
    picks = RELAX_QC_pick_examples(files, scores, overrides);
    [~, warnId] = lastwarn();
    assert(strcmp(warnId, 'RELAX_QC:unknownOverride'), 'an override naming no file warns');
    assert(isequal([picks.file], ["b", "a"]), 'override order, unknown file dropped');
    assert(all([picks.role] == "chosen"), 'override role');
    assert(picks(1).centres.blink == 12.5 && picks(1).centres.pop == 30 && isnan(picks(1).centres.muscle), 'centres from overrides');

    allMissing = struct('blinkAmplitude', nan(2, 1), 'musclePower', nan(2, 1), 'popSteps', nan(2, 1));
    assert(isempty(RELAX_QC_pick_examples(["x"; "y"], allMissing, [])), 'nothing scored, nothing picked');

    fprintf('test_RELAX_QC_pick_examples: PASS\n');
end
