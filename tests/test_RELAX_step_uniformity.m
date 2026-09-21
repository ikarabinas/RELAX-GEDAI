function test_RELAX_step_uniformity()
% Per-channel level change across one interval and how uniform it is across channels.
% A reference pop gives the same amplitude everywhere; a movement event does not.

    rng(11); n = 4000; W = 150; onset = 2000; offset = 2030;
    noise = 5; N = 32;

    %% 1. the same +100 uV step on every channel
    X = noise * randn(N, n); X(:, offset:end) = X(:, offset:end) + 100;
    [amp, s] = RELAX_step_uniformity(X, onset, offset, W, 1, n);
    assert(isequal(size(amp), [N 1]));
    assert(all(abs(amp - 100) < 3 * noise), 'each channel measures ~100 (max dev %.1f)', max(abs(amp - 100)));
    assert(abs(s.medianAmplitude - 100) < noise && s.nChannels == N);
    assert(s.signAgreement == 1 && s.fractionNear == 1);
    assert(s.amplitudeSpread < 0.15, 'spread %.3f', s.amplitudeSpread);

    %% 2. mixed-sign event: half the channels go up, half down
    X2 = noise * randn(N, n);
    X2(1:N/2, offset:end) = X2(1:N/2, offset:end) + 100;
    X2(N/2+1:end, offset:end) = X2(N/2+1:end, offset:end) - 100;
    [~, s2] = RELAX_step_uniformity(X2, onset, offset, W, 1, n);
    assert(s2.signAgreement <= 0.6, 'sign agreement %.2f', s2.signAgreement);

    %% 3. same sign, variable amplitude (ground-pop-like): agreement high, spread large
    X3 = noise * randn(N, n);
    X3(:, offset:end) = X3(:, offset:end) + linspace(20, 400, N)';
    [~, s3] = RELAX_step_uniformity(X3, onset, offset, W, 1, n);
    assert(s3.signAgreement == 1);
    assert(s3.amplitudeSpread > 0.5, 'spread %.3f', s3.amplitudeSpread);
    assert(s3.fractionNear < 0.5, 'fractionNear %.2f', s3.fractionNear);

    %% 4. windows are clipped to the segment without error
    [amp4, s4] = RELAX_step_uniformity(X, 20, 40, W, 1, n);
    assert(all(isfinite(amp4)) && isfinite(s4.medianAmplitude));
    [amp5, ~] = RELAX_step_uniformity(X, n - 40, n - 20, W, 1, n);
    assert(all(isfinite(amp5)));

    %% 5. a NaN channel is excluded from the summary, NaN in amp
    X6 = X; X6(3, :) = NaN;
    [amp6, s6] = RELAX_step_uniformity(X6, onset, offset, W, 1, n);
    assert(isnan(amp6(3)) && s6.nChannels == N - 1 && s6.signAgreement == 1);

    %% 6. single input, one channel
    [amp7, s7] = RELAX_step_uniformity(single(X(1, :)), onset, offset, W, 1, n);
    assert(isa(amp7, 'double') && s7.nChannels == 1 && s7.signAgreement == 1 && s7.amplitudeSpread == 0);

    fprintf('test_RELAX_step_uniformity: PASS\n');
end
