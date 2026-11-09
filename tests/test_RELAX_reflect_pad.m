function test_RELAX_reflect_pad()
% Even (mirror) reflection padding.

    % Mirror about the edge sample, which is not repeated
    x = [1 2 3 4 5];
    padded = RELAX_reflect_pad(x, 2);
    assert(isequal(padded, [3 2 1 2 3 4 5 4 3]), 'even reflection must mirror about the edge sample');

    % The original data must be preserved in the middle
    assert(isequal(padded(3:7), x), 'interior must be untouched');

    % An outlying edge sample must not shift the pad level: with odd reflection the pad
    % would sit near 2*x(1) - level, and a high-pass would ring on the step for ~1/fc
    x = [100 0 0 0 0 0];
    padded = RELAX_reflect_pad(x, 3);
    assert(isequal(padded(1:3), [0 0 0]), 'pad level must not follow an outlying edge sample');

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
    assert(isequal(padded(2,:), [30 20 10 20 30 40 50 40 30]), 'second channel reflected independently');

    fprintf('test_RELAX_reflect_pad: PASS\n');
end
