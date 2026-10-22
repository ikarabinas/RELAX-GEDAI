function test_RELAX_screen_loudest_channels()
    data = zeros(4, 1000); data(1, 1:100) = 50; data(3, 500:600) = 10; data(4, 500:600) = 5; data(2, 500:600) = 99;
    labels = RELAX_screen_loudest_channels(data, {'E1', 'E2', 'E3', 'E4'}, 450:650, {'E2'}, 2);
    assert(isequal(labels, {'E3', 'E4'}), 'only power inside the span counts; excluded labels skipped');
    labels = RELAX_screen_loudest_channels(data, {'E1', 'E2', 'E3', 'E4'}, 450:650, {'E2'}, 5);
    assert(numel(labels) == 3, 'no more labels than channels left after exclusion');
    fprintf('test_RELAX_screen_loudest_channels: PASS\n');
end
