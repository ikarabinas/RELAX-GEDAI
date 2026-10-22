% RELAX_QC_muscle_trace() - Where in a recording muscle activity is strongest.
%
% Usage:
%   >> [trace, channelPower, filtered] = RELAX_QC_muscle_trace(EEG, windowSeconds);
%
% Outputs:
%   trace        - [1 x pnts] channel-mean 20–80 Hz power, moving-averaged over
%                  windowSeconds within each boundary segment (never across one);
%                  NaN in segments too short to filter or smooth
%   channelPower - [nChannels x 1] each channel's mean 20–80 Hz power
%   filtered     - the 20–80 Hz band-passed EEG the two outputs above are computed from

function [trace, channelPower, filtered] = RELAX_QC_muscle_trace(EEG, windowSeconds)

    filtered = RELAX_filtbutter(EEG, 20, 80, 4, 'bandpass', 'acausal', true);

    width = round(windowSeconds * EEG.srate);
    segments = RELAX_get_segments(EEG);
    % RELAX_filtbutter passes segments under its 16-sample minimum through unfiltered, and a
    % segment shorter than the smoothing window has no meaningful smoothed value either:
    usable = (segments(:, 2) - segments(:, 1) + 1) >= max(width, 16);
    isUsable = false(1, EEG.pnts);
    for s = find(usable)'
        isUsable(segments(s, 1):segments(s, 2)) = true;
    end

    meanPower = zeros(1, EEG.pnts);
    channelPower = zeros(EEG.nbchan, 1);
    for c = 1:EEG.nbchan                                     % per channel: a full double copy of 256-channel data would not fit comfortably
        p = double(filtered.data(c, :)).^2;
        meanPower = meanPower + p;
        channelPower(c) = mean(p(isUsable));
    end
    meanPower = meanPower / EEG.nbchan;

    trace = nan(1, EEG.pnts);
    for s = find(usable)'
        span = segments(s, 1):segments(s, 2);
        trace(span) = movmean(meanPower(span), width, 'Endpoints', 'shrink');
    end
end
