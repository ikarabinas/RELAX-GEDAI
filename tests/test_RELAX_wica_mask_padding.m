function test_RELAX_wica_mask_padding()
% Blink mask dilation and smoothing, applied to one block at a time.

    % A single flagged sample must dilate symmetrically by moving_mean_length
    mask = zeros(1,101);
    mask(51) = 1;
    out = RELAX_pad_blink_mask_block(mask, 10);
    assert(isequal(size(out), size(mask)), 'output size must match input size');
    assert(all(out(41:61) > 0), 'the dilated region must be non-zero');
    assert(out(51) > 0, 'the flagged sample itself must stay non-zero');
    assert(out(1) == 0 && out(101) == 0, 'far-away samples must stay zero');

    % An all-zero mask stays all zero
    assert(all(RELAX_pad_blink_mask_block(zeros(1,101), 10) == 0), 'zero mask stays zero');

    % An all-ones mask stays all ones
    assert(all(abs(RELAX_pad_blink_mask_block(ones(1,101), 10) - 1) < 1e-12), 'ones mask stays ones');

    % ISOLATION through RELAX_apply_by_segments: a blink at the end of segment 1
    % must not bleed into segment 2.
    segments = [1 100; 101 200];
    mask = zeros(1,200);
    mask(98) = 1; % 2 samples from the end of segment 1
    out = RELAX_apply_by_segments(mask, segments, @(d) RELAX_pad_blink_mask_block(d, 10), 0);
    assert(all(out(101:200) == 0), 'dilation must not cross a segment boundary');
    assert(any(out(88:100) > 0), 'dilation must still occur inside segment 1');

    fprintf('test_RELAX_wica_mask_padding: PASS\n');
end
