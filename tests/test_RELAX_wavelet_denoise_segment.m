function test_RELAX_wavelet_denoise_segment()
% Level-5 coif5 stationary wavelet thresholding of one block of data.

    rng(42);
    n = 4001; % 4001 = 125*32 + 1, so the internal zero-padding branch is exercised
    t = (0:n-1)/500;
    x = sin(2*pi*10*t) + 0.1*randn(1,n);

    % Length is always preserved, including for non-multiples of 32
    y = RELAX_wavelet_denoise_segment(x, 0.5, 's');
    assert(isequal(size(y), size(x)), 'output length must equal input length');
    assert(~any(isnan(y)), 'output must contain no NaNs');
    assert(mod(numel(x), 32) ~= 0, 'the primary test length must NOT be a multiple of 32');

    % A length that IS an exact multiple of 32 takes the no-padding path
    x32 = x(1:3968); % 3968 = 124*32 exactly, so padding is empty
    y32 = RELAX_wavelet_denoise_segment(x32, 0.5, 's');
    assert(isequal(size(y32), size(x32)), 'output length preserved for multiple-of-32 input');

    % A zero threshold is a near-perfect reconstruction (swt/iswt round-trip)
    % Note: SWT/ISWT round-trip has floating-point precision limit ~1.5-1.8e-8
    y0 = RELAX_wavelet_denoise_segment(x, 0, 's');
    assert(max(abs(y0 - x)) < 2e-8, 'zero threshold must reconstruct the input');

    % A larger threshold removes more energy than a smaller one
    y_small = RELAX_wavelet_denoise_segment(x, 0.1, 's');
    y_large = RELAX_wavelet_denoise_segment(x, 1.0, 's');
    assert(sum(y_large.^2) < sum(y_small.^2), 'a larger threshold must remove more energy');

    % Hard thresholding is also accepted
    yh = RELAX_wavelet_denoise_segment(x, 0.5, 'h');
    assert(isequal(size(yh), size(x)), 'hard threshold type accepted');

    % Determinism
    assert(isequal(RELAX_wavelet_denoise_segment(x, 0.5, 's'), y), 'must be deterministic');

    fprintf('test_RELAX_wavelet_denoise_segment: PASS\n');
end
