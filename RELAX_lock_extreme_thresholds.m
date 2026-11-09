function locked = RELAX_lock_extreme_thresholds(epochedData, labels, RELAX_cfg)
% Per-channel statistics of the first extreme-marking pass, fixed so that cleaned data can
% later be judged by the same standard (RELAX_score_locked_extremes) rather than by
% statistics re-fitted on cleaner data.
%
% epochedData: channels x samples x epochs, fitted exactly as given (the caller chooses
% the reference). labels: 1 x channels cellstr, used to match channels later. Channels are
% taken one at a time in double, so the input is never copied whole.
%
% Returned, per channel unless noted:
%   shiftUpper              median + ExtremeVoltageShiftThreshold * MAD of the epoch
%                           voltage range (no blink-aware bound)
%   jpLo, jpHi, jpHistogram the value range and 1000-bin histogram EEGLAB's jointprob
%                           builds, and jpFloor (scalar) for values it never saw
%   jpScoreMean, jpScoreStd mean and SD over epochs of the epoch score -sum(log p)
%   kurtMean, kurtStd       mean and SD over epochs of the epoch kurtosis
% A channel with no spread in the first pass (SD 0 or NaN) is not judged by joint
% probability or kurtosis later.

    [nChan, nPnts, nEpochs] = size(epochedData);
    bins = 1000;

    locked.labels = labels;
    voltageRange = zeros(nChan, nEpochs);
    kurt = zeros(nChan, nEpochs);
    locked.jpLo = zeros(nChan, 1);
    locked.jpHi = zeros(nChan, 1);
    locked.jpHistogram = zeros(nChan, bins);
    for c = 1:nChan
        channel = double(reshape(epochedData(c, :, :), nPnts, nEpochs));
        voltageRange(c, :) = range(channel, 1);
        kurt(c, :) = kurtosis(channel, 1, 1);
        locked.jpLo(c) = min(channel, [], 'all');
        locked.jpHi(c) = max(channel, [], 'all');
        b = RELAX_jointprob_bins(channel(:), locked.jpLo(c), locked.jpHi(c), bins);
        locked.jpHistogram(c, :) = accumarray(b, 1, [bins 1])' / numel(channel);
    end
    locked.shiftUpper = median(voltageRange, 2) + RELAX_cfg.ExtremeVoltageShiftThreshold * mad(voltageRange, 1, 2);
    locked.kurtMean = mean(kurt, 2);
    locked.kurtStd = std(kurt, 0, 2);
    locked.jpFloor = 1 / (nPnts * nEpochs);
    locked.jpScoreMean = zeros(nChan, 1);
    locked.jpScoreStd = ones(nChan, 1);

    rawScores = RELAX_score_locked_extremes(locked, epochedData, labels).jpRaw;
    locked.jpScoreMean = mean(rawScores, 2);
    locked.jpScoreStd = std(rawScores, 0, 2);
end
