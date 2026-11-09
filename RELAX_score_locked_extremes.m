function scores = RELAX_score_locked_extremes(locked, epochedData, labels)
% Scores epochs against first-pass statistics from RELAX_lock_extreme_thresholds, matching
% each channel to its own first-pass row by label. Nothing is re-fitted, so an epoch scores
% the same whatever other epochs it is scored with. Channels are taken one at a time in
% double, so the input is never copied whole.
%
% Returned, channels x epochs unless noted:
%   range, maxAbs  voltage range and largest |voltage| within the epoch
%   jpRaw, jpZ     joint-probability score -sum(log p) and its z against the first pass
%   kurtZ          kurtosis z against the first pass
%   shiftUpper     channels x 1, the matched first-pass voltage-shift bounds
% A channel with no first-pass spread in a statistic scores z = 0 on it.

    [known, rows] = ismember(labels, locked.labels);
    if ~all(known)
        error('RELAX:secondPassUnknownChannel', 'Channels the first extreme pass never scored: %s', ...
            strjoin(labels(~known), ', '));
    end
    [nChan, nPnts, nEpochs] = size(epochedData);
    bins = size(locked.jpHistogram, 2);

    scores.range = zeros(nChan, nEpochs);
    scores.maxAbs = zeros(nChan, nEpochs);
    scores.jpRaw = zeros(nChan, nEpochs);
    kurt = zeros(nChan, nEpochs);
    for c = 1:nChan
        r = rows(c);
        channel = double(reshape(epochedData(c, :, :), nPnts, nEpochs));
        scores.range(c, :) = range(channel, 1);
        scores.maxAbs(c, :) = max(abs(channel), [], 1);
        kurt(c, :) = kurtosis(channel, 1, 1);
        b = RELAX_jointprob_bins(channel, locked.jpLo(r), locked.jpHi(r), bins);
        p = repmat(locked.jpFloor, size(b));
        inRange = b >= 1 & b <= bins;
        p(inRange) = max(locked.jpHistogram(r, b(inRange)), locked.jpFloor);
        scores.jpRaw(c, :) = -sum(log(p), 1);
    end
    scores.shiftUpper = locked.shiftUpper(rows);
    scores.jpZ = z_against_first_pass(scores.jpRaw, locked.jpScoreMean(rows), locked.jpScoreStd(rows));
    scores.kurtZ = z_against_first_pass(kurt, locked.kurtMean(rows), locked.kurtStd(rows));
end

function z = z_against_first_pass(x, firstPassMean, firstPassStd)
    z = (x - firstPassMean) ./ firstPassStd;
    z(~(firstPassStd > 0), :) = 0;
end
