function test_RELAX_screen_top_separated()
    times = 0:0.1:10; scores = zeros(size(times));
    at = @(s) abs(times - s) < 1e-9;
    scores(at(2)) = 10; scores(at(2.5)) = 9; scores(at(6)) = 8; scores(at(9.9)) = 7; scores(at(5)) = NaN;
    c = RELAX_screen_top_separated(times, scores, 3, 3);
    assert(isequal(round(c, 6), [2 6 9.9]), 'the 2.5 s peak is suppressed by the 2 s one');

    c = RELAX_screen_top_separated([1 1.5 8], [5 4 3], 3, 3);
    assert(isequal(c, [1 8]), 'fewer than k separated candidates returns fewer');

    c = RELAX_screen_top_separated([1 5 9 13], [4 3 2 1], 2, 3);
    assert(isequal(c, [1 5]), 'more separated candidates than k returns the k highest');

    assert(isempty(RELAX_screen_top_separated([], [], 3, 3)));
    assert(isempty(RELAX_screen_top_separated([1 2], [NaN NaN], 3, 3)));
    fprintf('test_RELAX_screen_top_separated: PASS\n');
end
