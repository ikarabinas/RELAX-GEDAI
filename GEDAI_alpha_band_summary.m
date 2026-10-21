function B = GEDAI_alpha_band_summary(A, M, srate, opts, varargin)
% GEDAI_ALPHA_BAND_SUMMARY  Rest-block spectra, band powers, alpha statistics of the three GEDAI arms.
%
%   B = GEDAI_alpha_band_summary(A, M, srate, opts, chanlocs)
%
%   A    struct with off, on, removed (nch x pnts) from GEDAI_alpha_run
%   M    struct from GEDAI_alpha_masks; only M.rest is used
%
%   B.f, B.P.(off|on|removed) nch x nf, B.band.(...) nch x nbands (power summed over bins),
%   B.topAlpha (the nTopAlpha channels with the highest off 8-13 Hz power), B.iaf,
%   B.alphaRatioTop (median over topAlpha of on/off alpha power), B.alphaPowerTop,
%   B.alphaRelTop (alpha power / 4-45 Hz power, off arm, mean over topAlpha), and
%   B.vecInput / B.vecRemoved, the leading eigenvectors of the 8-13 Hz band-passed
%   off and removed covariances over the rest samples.
%
%   Second ROI: B.topAlphaPeak, the nTopAlpha channels with the highest off-arm alpha
%   prominence (mean 8-13 Hz power / mean of 4-7 and 14-20 Hz power), with B.iafPeak,
%   B.alphaRatioPeak, B.alphaPowerPeak. B.iafAtEdge / B.iafPeakAtEdge are true when the IAF lies
%   within one bin of either end of opts.iafRange. B.topAlphaCentroid / B.topAlphaPeakCentroid
%   are the 1x3 mean [X Y Z] of the optional chanlocs over each ROI (NaN without chanlocs).

    p = inputParser;
    addRequired(p, 'A', @isstruct);
    addRequired(p, 'M', @isstruct);
    addRequired(p, 'srate', @isnumeric);
    addRequired(p, 'opts', @isstruct);
    addOptional(p, 'chanlocs', [], @(c) isempty(c) || isstruct(c));
    parse(p, A, M, srate, opts, varargin{:});
    chanlocs = p.Results.chanlocs;

    arms = {'off', 'on', 'removed'};
    for a = 1:numel(arms)
        [B.f, B.P.(arms{a})] = GEDAI_alpha_mask_psd(A.(arms{a}), srate, M.rest, opts.winSec, opts.maxHz);
    end

    binWidth = B.f(2) - B.f(1);
    for a = 1:numel(arms)
        for b = 1:size(opts.bands, 1)
            inBand = B.f >= opts.bands(b, 1) & B.f < opts.bands(b, 2);
            B.band.(arms{a})(:, b) = sum(B.P.(arms{a})(:, inBand), 2) * binWidth;
        end
    end

    alphaBins = B.f >= opts.alphaBand(1) & B.f < opts.alphaBand(2);
    offAlpha = sum(B.P.off(:, alphaBins), 2) * binWidth;
    onAlpha = sum(B.P.on(:, alphaBins), 2) * binWidth;
    broadBins = B.f >= 4 & B.f < 45;
    offBroad = sum(B.P.off(:, broadBins), 2) * binWidth;
    shoulderBins = (B.f >= 4 & B.f < 8) | (B.f >= 14 & B.f < 20);
    prominence = mean(B.P.off(:, alphaBins), 2) ./ mean(B.P.off(:, shoulderBins), 2);

    [~, byPower] = sort(offAlpha, 'descend');
    [~, byProminence] = sort(prominence, 'descend');
    top = alphaRoi(byPower(1:opts.nTopAlpha), B, offAlpha, onAlpha, offBroad, opts, binWidth, chanlocs);
    peak = alphaRoi(byProminence(1:opts.nTopAlpha), B, offAlpha, onAlpha, offBroad, opts, binWidth, chanlocs);

    B.topAlpha = top.channels;
    B.iaf = top.iaf;
    B.iafAtEdge = top.iafAtEdge;
    B.alphaRatioTop = top.ratio;
    B.alphaPowerTop = top.power;
    B.alphaRelTop = top.relative;
    B.topAlphaCentroid = top.centroid;
    B.topAlphaPeak = peak.channels;
    B.iafPeak = peak.iaf;
    B.iafPeakAtEdge = peak.iafAtEdge;
    B.alphaRatioPeak = peak.ratio;
    B.alphaPowerPeak = peak.power;
    B.alphaRelPeak = peak.relative;
    B.topAlphaPeakCentroid = peak.centroid;

    B.vecInput = leadingAlphaVector(A.off, M.rest, srate, opts.alphaBand);
    B.vecRemoved = leadingAlphaVector(A.removed, M.rest, srate, opts.alphaBand);
end

function vector = leadingAlphaVector(data, mask, srate, band)
    [b, a] = butter(2, band / (srate / 2), 'bandpass');
    paddedMask = [false, mask(:)', false];
    runs = [find(diff(paddedMask) == 1)', find(diff(paddedMask) == -1)' - 1];
    minRun = 3 * max(numel(a), numel(b));
    nch = size(data, 1);
    covariance = zeros(nch);
    for r = 1:size(runs, 1)
        if runs(r, 2) - runs(r, 1) + 1 <= minRun, continue, end
        x = filtfilt(b, a, double(data(:, runs(r, 1):runs(r, 2)))')';
        covariance = covariance + x * x';
    end
    [vectors, values] = eig((covariance + covariance') / 2);
    [~, last] = max(diag(values));
    vector = vectors(:, last);
end

function roi = alphaRoi(channels, B, offAlpha, onAlpha, offBroad, opts, binWidth, chanlocs)
    roi.channels = channels;
    iafBins = find(B.f >= opts.iafRange(1) & B.f <= opts.iafRange(2));
    [~, peak] = max(mean(B.P.off(channels, iafBins), 1));
    roi.iaf = B.f(iafBins(peak));
    roi.iafAtEdge = roi.iaf <= opts.iafRange(1) + binWidth || roi.iaf >= opts.iafRange(2) - binWidth;
    roi.ratio = median(onAlpha(channels) ./ offAlpha(channels));
    roi.power = mean(offAlpha(channels));
    roi.relative = roi.power / mean(offBroad(channels));
    roi.centroid = nan(1, 3);
    if ~isempty(chanlocs)
        roi.centroid = [mean([chanlocs(channels).X]), mean([chanlocs(channels).Y]), mean([chanlocs(channels).Z])];
    end
end
