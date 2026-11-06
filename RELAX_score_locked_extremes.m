function scores = RELAX_score_locked_extremes(locked, epochedData, labels)
% Scores epochs against first-pass statistics from RELAX_lock_extreme_thresholds, matching
% each channel to its own first-pass row by label. Nothing is re-fitted, so an epoch scores
% the same whatever other epochs it is scored with.
%
% Returned, channels x epochs unless noted:
%   range, maxAbs  voltage range and largest |voltage| within the epoch
%   jpRaw, jpZ     joint-probability score -sum(log p) and its z against the first pass
%   kurtZ          kurtosis z against the first pass
%   shiftUpper     channels x 1, the matched first-pass voltage-shift bounds

    [known, rows] = ismember(labels, locked.labels);
    if ~all(known)
        error('RELAX:secondPassUnknownChannel', 'Channels the first extreme pass never scored: %s', ...
            strjoin(labels(~known), ', '));
    end
    data = double(epochedData);
    [nChan, ~, nEpochs] = size(data);

    scores.range = reshape(range(data, 2), nChan, nEpochs);
    scores.maxAbs = reshape(max(abs(data), [], 2), nChan, nEpochs);
    scores.shiftUpper = locked.shiftUpper(rows);

    bins = size(locked.jpHistogram, 2);
    scores.jpRaw = zeros(nChan, nEpochs);
    for c = 1:nChan
        r = rows(c);
        b = RELAX_jointprob_bins(reshape(data(c, :, :), [], nEpochs), locked.jpLo(r), locked.jpHi(r), bins);
        p = repmat(locked.jpFloor, size(b));
        inRange = b >= 1 & b <= bins;
        p(inRange) = max(locked.jpHistogram(r, b(inRange)), locked.jpFloor);
        scores.jpRaw(c, :) = -sum(log(p), 1);
    end
    scores.jpZ = (scores.jpRaw - locked.jpScoreMean(rows)) ./ locked.jpScoreStd(rows);

    kurt = reshape(kurtosis(data, 1, 2), nChan, nEpochs);
    scores.kurtZ = (kurt - locked.kurtMean(rows)) ./ locked.kurtStd(rows);
end
