function test_RELAX_reflect_pad()
% Odd reflection padding.

    % A linear ramp must be continued exactly by odd reflection
    x = [1 2 3 4 5];
    padded = RELAX_reflect_pad(x, 2);
    assert(isequal(padded, [-1 0 1 2 3 4 5 6 7]), 'odd reflection must continue a linear ramp');

    % The original data must be preserved in the middle
    assert(isequal(padded(3:7), x), 'interior must be untouched');

    % p is clamped to nSamples-1
    x = [1 2 3];
    padded = RELAX_reflect_pad(x, 10);
    assert(size(padded,2) == 3 + 2*2, 'p must clamp to nSamples-1');

    % p of zero returns the input unchanged
    x = [1 2 3];
    assert(isequal(RELAX_reflect_pad(x, 0), x), 'p=0 returns input unchanged');

    % Multi-channel shape is preserved
    x = [1 2 3 4 5; 10 20 30 40 50];
    padded = RELAX_reflect_pad(x, 2);
    assert(isequal(size(padded), [2 9]), 'multi-channel output shape');
    assert(isequal(padded(2,:), [-10 0 10 20 30 40 50 60 70]), 'second channel reflected independently');

    fprintf('test_RELAX_reflect_pad: PASS\n');
end
