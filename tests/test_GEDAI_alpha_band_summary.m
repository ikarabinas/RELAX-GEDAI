function test_GEDAI_alpha_band_summary()
% GEDAI_alpha_band_summary on a synthetic file where the answer is known. A 10 Hz component
% with one shared time course and distinct per-channel amplitudes sits on channels 1-20 of 32;
% the removed arm is half of it, so the on arm keeps a quarter of its alpha power.

    opts = GEDAI_alpha_figures('opts');
    srate = 250;
    nChannels = 32;
    nSamples = 120 * srate;
    nAlphaChannels = 20;
    noiseSigma = 0.1;

    stream = RandStream('mt19937ar', 'Seed', 4);
    time = (0:nSamples - 1) / srate;
    alphaAmplitudes = zeros(nChannels, 1);
    alphaAmplitudes(1:nAlphaChannels) = 1 + 0.1 * (1:nAlphaChannels)';
    alphaComponent = alphaAmplitudes * sin(2 * pi * 10 * time);
    noise = noiseSigma * randn(stream, nChannels, nSamples);

    A.off = noise + alphaComponent;
    A.removed = 0.5 * alphaComponent;
    A.on = A.off - A.removed;
    A.threshold = 0;
    A.enova = 0;
    A.seconds = 0;

    EEG = eeg_emptyset();
    EEG.srate = srate;
    EEG.nbchan = nChannels;
    EEG.pnts = nSamples;
    EEG.xmax = (nSamples - 1) / srate;
    EEG.event = [];
    M = GEDAI_alpha_masks(EEG, opts);
    assert(M.restIsWholeFile, 'fixture has no TMS_train events, so rest is the whole file');

    B = GEDAI_alpha_band_summary(A, M, srate, opts);

    assert(abs(B.iaf - 10) <= 0.5, 'iaf is %g Hz, expected 10 +/- 0.5', B.iaf);
    assert(abs(B.alphaRatioTop - 0.25) <= 0.03, ...
        'alphaRatioTop is %g, expected 0.25 +/- 0.03', B.alphaRatioTop);
    assert(isequal(sort(B.topAlpha(:)), (1:nAlphaChannels)'), ...
        'topAlpha must be channels 1-20, got %s', mat2str(sort(B.topAlpha(:))'));
    vectorCorrelation = abs(corr(B.vecRemoved(:), B.vecInput(:)));
    assert(vectorCorrelation > 0.95, ...
        '|corr(vecRemoved, vecInput)| is %g, expected > 0.95', vectorCorrelation);

    assert(size(B.band.off, 1) == nChannels && size(B.band.off, 2) == size(opts.bands, 1), ...
        'band.off must be nch x 6');
    assert(size(B.P.removed, 1) == nChannels && size(B.P.removed, 2) == numel(B.f), ...
        'P.removed must be nch x nf');

    fprintf('test_GEDAI_alpha_band_summary: PASS\n');
end
