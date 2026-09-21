function test_RELAX_flat_channels()
% A channel is flat when it has no variation inside ANY usable segment: globally constant,
% railed at a different constant per segment, or all-NaN. One varying segment is enough to
% keep it.

    n = 1000;
    segments = [1 500; 501 n];
    data = zeros(5, n);
    data(1, :) = randn(1, n);              % ordinary
    data(2, :) = 7;                        % globally constant
    data(3, 501:end) = 333333;             % railed: 0 then a different constant
    data(4, :) = NaN;                      % non-finite everywhere
    data(5, 1:500) = 3; data(5, 501:end) = randn(1, 500);   % constant in one segment only

    flat = RELAX_flat_channels(data, segments);
    assert(islogical(flat) && isequal(size(flat), [5 1]), 'nbchan x 1 logical');
    assert(isequal(flat, [false; true; true; true; false]), 'got %s', mat2str(flat'));

    flat1 = RELAX_flat_channels(data(3, :), [1 n]);
    assert(~flat1, 'over one segment the railed channel does vary');

    flat0 = RELAX_flat_channels(zeros(0, n), segments);
    assert(islogical(flat0) && isequal(size(flat0), [0 1]));

    fprintf('test_RELAX_flat_channels: PASS\n');
end
