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

    % A zero threshold is a near-perfect reconstruction (swt/iswt round-trip).
    % MATLAB's level-5 coif5 swt/iswt round-trip has a precision floor around
    % sqrt(eps) (measured: 1.4e-8 for a bare iswt(swt(x)) with no padding at all,
    % and 1.6e-8 for an exact-multiple-of-32 length), so the tolerance is keyed to
    % that floor with an order of magnitude of headroom rather than to a bare
    % constant that happens to pass on one machine.
    y0 = RELAX_wavelet_denoise_segment(x, 0, 's');
    reconstruction_tolerance = 10 * sqrt(eps) * max(abs(x));
    assert(max(abs(y0 - x)) < reconstruction_tolerance, ...
        'zero threshold must reconstruct the input (max err %.3e, tolerance %.3e)', ...
        max(abs(y0 - x)), reconstruction_tolerance);

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
