function [EEG, changed] = RELAX_interpolate_channel_epochs(EEG, flagged, onsets, epochSamples, rampSeconds)
% Replaces each flagged channel over each flagged epoch with a spherical spline of that
% epoch's unflagged channels, blended into the real signal by a raised cosine over the
% first and last rampSeconds inside the epoch.
%
% flagged: channels x epochs logical. onsets: first sample of each epoch on EEG's
% continuous timeline. changed: channels x samples, true wherever data was altered.
%
% Epochs are applied in onset order on the progressively repaired data, so where two
% flagged epochs overlap, the later one's fully replaced middle covers the earlier one's
% closing ramp.

    ramp = round(rampSeconds * EEG.srate);
    up = 0.5 * (1 - cos(pi * (1:ramp) / (ramp + 1)));
    weight = [up, ones(1, epochSamples - 2 * ramp), fliplr(up)];
    dataClass = class(EEG.data);
    changed = false(EEG.nbchan, EEG.pnts);

    [~, order] = sort(onsets);
    for e = order(any(flagged(:, order), 1))
        span = onsets(e):onsets(e) + epochSamples - 1;
        bad = find(flagged(:, e));
        current = double(EEG.data(:, span));
        interpolated = spline_interpolate(current, bad, EEG.chanlocs, EEG.srate);
        EEG.data(bad, span) = cast(weight .* interpolated(bad, :) + (1 - weight) .* current(bad, :), dataClass);
        changed(bad, span) = true;
    end
end

function data = spline_interpolate(data, bad, chanlocs, srate)
    epoch = eeg_emptyset;
    epoch.data = data;
    epoch.nbchan = size(data, 1); epoch.pnts = size(data, 2); epoch.trials = 1;
    epoch.srate = srate; epoch.xmin = 0; epoch.xmax = (epoch.pnts - 1) / srate;
    epoch.chanlocs = chanlocs;
    epoch = pop_interp(epoch, bad, 'spherical');
    data = double(epoch.data);
end
