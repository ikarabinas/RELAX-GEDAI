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

    peakVersusRawPowerRois(opts, M, srate);

    fprintf('test_GEDAI_alpha_band_summary: PASS\n');
end

function peakVersusRawPowerRois(opts, M, srate)
% Channels 1-20 carry a 10 Hz alpha peak. Channels 21-40 carry steep 1/f noise scaled to
% four times their raw 8-13 Hz power and no peak. The raw-power ROI (topAlpha) must pick
% 21-40 and sit on the 7 Hz edge of iafRange; the prominence ROI (topAlphaPeak) must pick
% 1-20 and find 10 Hz. Centroids need chanlocs, passed as the optional 5th argument.
    nChannels = 40;
    nAlphaChannels = 20;
    nSamples = size(M.rest, 2);
    stream = RandStream('mt19937ar', 'Seed', 5);
    time = (0:nSamples - 1) / srate;

    alphaComponent = zeros(nChannels, nSamples);
    alphaComponent(1:nAlphaChannels, :) = repmat(sin(2 * pi * 10 * time), nAlphaChannels, 1);
    alphaBandPower = 0.5;
    rawPowerRatio = 4;

    pinkNoise = steepPinkNoise(stream, nChannels - nAlphaChannels, nSamples, srate, [8 13], rawPowerRatio * alphaBandPower);
    whiteNoise = 0.05 * randn(stream, nChannels, nSamples);
    A.off = whiteNoise + alphaComponent;
    A.off(nAlphaChannels + 1:end, :) = A.off(nAlphaChannels + 1:end, :) + pinkNoise;
    A.removed = 0.5 * alphaComponent;
    A.on = A.off - A.removed;
    A.threshold = 0;
    A.enova = 0;
    A.seconds = 0;

    B = GEDAI_alpha_band_summary(A, M, srate, opts);

    peakRoi = 1:nAlphaChannels;
    rawRoi = nAlphaChannels + 1:nChannels;
    assert(isfield(B, 'topAlphaPeak'), 'B.topAlphaPeak is missing');
    assert(numel(B.topAlphaPeak) == opts.nTopAlpha, 'topAlphaPeak must hold opts.nTopAlpha indices');
    assert(isequal(sort(B.topAlphaPeak(:)), peakRoi'), ...
        'topAlphaPeak must be channels 1-20, got %s', mat2str(sort(B.topAlphaPeak(:))'));
    assert(all(B.topAlpha >= rawRoi(1) & B.topAlpha <= rawRoi(end)), ...
        'topAlpha must lie within 21-40 (raw power), got %s', mat2str(sort(B.topAlpha(:))'));

    assert(abs(B.iafPeak - 10) <= 0.5, 'iafPeak is %g Hz, expected 10 +/- 0.5', B.iafPeak);
    assert(islogical(B.iafPeakAtEdge) && ~B.iafPeakAtEdge, 'iafPeakAtEdge must be a false logical (iafPeak %g Hz)', B.iafPeak);
    assert(islogical(B.iafAtEdge) && B.iafAtEdge, ...
        'iafAtEdge must be true: iaf over the 1/f ROI is %g Hz, iafRange %s', B.iaf, mat2str(opts.iafRange));

    assert(abs(B.alphaRatioPeak - 0.25) <= 0.03, 'alphaRatioPeak is %g, expected 0.25 +/- 0.03', B.alphaRatioPeak);
    assert(abs(B.alphaRatioTop - 1) <= 1e-9, ...
        'removed is zero on channels 21-40, so alphaRatioTop must be 1, got %g', B.alphaRatioTop);
    assert(B.alphaPowerPeak > 0 && B.alphaPowerTop > 0, 'alpha powers over the ROIs must be positive');

    assert(isequal(size(B.topAlphaCentroid), [1 3]) && all(isnan(B.topAlphaCentroid)), ...
        'topAlphaCentroid must be 1x3 NaN when no chanlocs are given');
    assert(isequal(size(B.topAlphaPeakCentroid), [1 3]) && all(isnan(B.topAlphaPeakCentroid)), ...
        'topAlphaPeakCentroid must be 1x3 NaN when no chanlocs are given');

    chanlocs = struct('X', num2cell(1:nChannels), ...
        'Y', num2cell(mod((1:nChannels) * 7, 13)), ...
        'Z', num2cell(0.5 * (1:nChannels) + 3));
    Bloc = GEDAI_alpha_band_summary(A, M, srate, opts, chanlocs);
    assert(isequal(Bloc.topAlpha, B.topAlpha) && isequal(Bloc.topAlphaPeak, B.topAlphaPeak), ...
        'passing chanlocs must not change the ROIs');
    assert(max(abs(Bloc.topAlphaPeakCentroid - roiCentroid(chanlocs, peakRoi))) < 1e-9, ...
        'topAlphaPeakCentroid is %s, expected %s', mat2str(Bloc.topAlphaPeakCentroid, 6), ...
        mat2str(roiCentroid(chanlocs, peakRoi), 6));
    assert(max(abs(Bloc.topAlphaCentroid - roiCentroid(chanlocs, rawRoi))) < 1e-9, ...
        'topAlphaCentroid is %s, expected %s', mat2str(Bloc.topAlphaCentroid, 6), ...
        mat2str(roiCentroid(chanlocs, rawRoi), 6));
end

function centroid = roiCentroid(chanlocs, roi)
    centroid = [mean([chanlocs(roi).X]), mean([chanlocs(roi).Y]), mean([chanlocs(roi).Z])];
end

function noise = steepPinkNoise(stream, nChannels, nSamples, srate, band, bandPower)
% Gaussian noise with amplitude spectrum ~ f^-1.5 (power ~ f^-3, monotonically falling, so
% no spectral peak anywhere), each channel scaled so its power inside band is bandPower.
    frequencies = (0:nSamples - 1)' * srate / nSamples;
    frequencies(frequencies > srate / 2) = srate - frequencies(frequencies > srate / 2);
    amplitude = zeros(nSamples, 1);
    amplitude(frequencies > 0.5) = frequencies(frequencies > 0.5) .^ -1.5;
    noise = zeros(nChannels, nSamples);
    inBand = frequencies >= band(1) & frequencies < band(2);
    for c = 1:nChannels
        shaped = real(ifft(fft(randn(stream, nSamples, 1)) .* amplitude));
        spectrum = fft(shaped);
        measuredBandPower = sum(abs(spectrum(inBand)) .^ 2) / nSamples ^ 2;
        noise(c, :) = (shaped * sqrt(bandPower / measuredBandPower))';
    end
end
