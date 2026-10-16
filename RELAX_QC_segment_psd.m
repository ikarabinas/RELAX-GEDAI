% RELAX_QC_segment_psd() - Welch power spectral density from windows that lie wholly
%                          inside one boundary segment.
%
% Usage:
%   >> [freqs, power] = RELAX_QC_segment_psd(data, srate, segments, windowSeconds, maxHz);
%
% Hann windows with 50 % overlap; every qualifying window across all segments counts
% equally. A segment shorter than one window contributes nothing, so excised data,
% boundaries and slivers never leak into the spectrum.
%
% Outputs:
%   freqs - [1 x F], 0..maxHz in steps of 1/windowSeconds
%   power - [nChannels x F] one-sided density; all NaN when no window qualifies

function [freqs, power] = RELAX_QC_segment_psd(data, srate, segments, windowSeconds, maxHz)

    windowLength = round(windowSeconds * srate);
    hop = floor(windowLength / 2);
    taper = hann(windowLength, 'periodic')';
    allFreqs = (0:floor(windowLength / 2)) * srate / windowLength;
    keep = allFreqs <= maxHz;
    freqs = allFreqs(keep);
    nBins = numel(allFreqs);

    summed = zeros(size(data, 1), nBins);
    nWindows = 0;
    for s = 1:size(segments, 1)
        for first = segments(s, 1) : hop : segments(s, 2) - windowLength + 1
            x = double(data(:, first:first + windowLength - 1));
            spectrum = fft((x - mean(x, 2)) .* taper, [], 2);
            summed = summed + abs(spectrum(:, 1:nBins)).^2;
            nWindows = nWindows + 1;
        end
    end

    if nWindows == 0
        power = nan(size(data, 1), numel(freqs));
        return
    end
    density = summed / nWindows / (srate * sum(taper.^2));
    density(:, 2:end-1) = 2 * density(:, 2:end-1);   % fold negative frequencies; DC and Nyquist have no mirror
    power = density(:, keep);
end
