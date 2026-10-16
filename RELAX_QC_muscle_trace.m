% RELAX_QC_muscle_trace() - Where in a recording muscle activity is strongest.
%
% Usage:
%   >> [trace, channelPower] = RELAX_QC_muscle_trace(EEG, windowSeconds);
%
% Outputs:
%   trace        - [1 x pnts] channel-mean 20–80 Hz power, moving-averaged over
%                  windowSeconds within each boundary segment (never across one)
%   channelPower - [nChannels x 1] each channel's mean 20–80 Hz power

function [trace, channelPower] = RELAX_QC_muscle_trace(EEG, windowSeconds)

    filtered = RELAX_filtbutter(EEG, 20, 80, 4, 'bandpass', 'acausal', true);

    meanPower = zeros(1, EEG.pnts);
    channelPower = zeros(EEG.nbchan, 1);
    for c = 1:EEG.nbchan                                     % per channel: a full double copy of 256-channel data would not fit comfortably
        p = double(filtered.data(c, :)).^2;
        meanPower = meanPower + p;
        channelPower(c) = mean(p);
    end
    meanPower = meanPower / EEG.nbchan;

    trace = nan(1, EEG.pnts);
    width = round(windowSeconds * EEG.srate);
    segments = RELAX_get_segments(EEG);
    for s = 1:size(segments, 1)
        span = segments(s, 1):segments(s, 2);
        trace(span) = movmean(meanPower(span), width, 'Endpoints', 'shrink');
    end
end
