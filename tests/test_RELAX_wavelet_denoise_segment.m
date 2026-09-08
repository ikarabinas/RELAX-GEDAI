function test_RELAX_wavelet_denoise_segment()
% Level-5 coif5 stationary wavelet thresholding of one block of data.

    rng(42);
    n = 4000; % not a multiple of 32, so internal padding is exercised
    t = (0:n-1)/500;
    x = sin(2*pi*10*t) + 0.1*randn(1,n);

    % Length is always preserved, including for non-multiples of 32
    y = RELAX_wavelet_denoise_segment(x, 0.5, 's');
    assert(isequal(size(y), size(x)), 'output length must equal input length');
    assert(~any(isnan(y)), 'output must contain no NaNs');

    % A length that IS a multiple of 32 also works
    x32 = x(1:3968); % 3968 = 124*32
    y32 = RELAX_wavelet_denoise_segment(x32, 0.5, 's');
    assert(isequal(size(y32), size(x32)), 'output length preserved for multiple-of-32 input');

    % A zero threshold is a near-perfect reconstruction (swt/iswt round-trip)
    y0 = RELAX_wavelet_denoise_segment(x, 0, 's');
    assert(max(abs(y0 - x)) < 1e-8, 'zero threshold must reconstruct the input');

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
