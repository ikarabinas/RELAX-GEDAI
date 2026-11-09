function test_RELAX_flag_extreme_channel_epochs()
% Each per-channel criterion of the second extreme pass fires on its own, joint
% probability and kurtosis only on the extreme side, and an epoch with more than the
% allowed proportion of its channels flagged stops the file.

    cfg = struct('ExtremeAbsoluteVoltageThreshold', 1000, ...
        'ExtremeImprobableVoltageDistributionThreshold', 10, ...
        'ExtremeSingleChannelKurtosisThreshold', 10, ...
        'ExtremeDriftSlopeThreshold', -4);
    onsetSeconds = [0 0.5 1 1.5];

    [s, drift] = quiet_scores();
    s.range(1, 1) = 150;
    s.range(2, 1) = 1.5;
    s.maxAbs(3, 1) = 1001;
    s.jpZ(4, 1) = 10.5;
    drift(5, 1) = -4.5;
    s.kurtZ(6, 1) = 10.5;
    s.jpZ(7, 2) = -50;
    s.kurtZ(8, 2) = -50;
    criteriaCfg = cfg;
    criteriaCfg.SecondPassMaxProportionFlaggedChannels = 1;  % six criteria share epoch 1
    flags = RELAX_flag_extreme_channel_epochs(s, drift, onsetSeconds, criteriaCfg);

    criteria = {'ShiftHigh', 'ShiftLow', 'AbsVolt', 'JointProb', 'Drift', 'Kurtosis'};
    for k = 1:numel(criteria)
        expected = false(10, 4);
        expected(k, 1) = true;
        assert(isequal(flags.(criteria{k}), expected), '%s fired somewhere other than channel %d epoch 1', criteria{k}, k);
    end
    expectedAny = false(10, 4);
    expectedAny(1:6, 1) = true;
    assert(isequal(flags.Any, expectedAny), 'Any is not the OR of the criteria');
    assert(~flags.Any(7, 2), 'a too-probable epoch (jp z = -50) was flagged');
    assert(~flags.Any(8, 2), 'a very negative kurtosis z was flagged');

    % Exactly 20 % of channels flagged does not raise
    [s, drift] = quiet_scores();
    s.kurtZ(1:2, 3) = 11;
    RELAX_flag_extreme_channel_epochs(s, drift, onsetSeconds, cfg);

    % More than 20 % raises, naming the epoch and the counts
    s.kurtZ(1:3, 4) = 11;
    threw = false;
    try
        RELAX_flag_extreme_channel_epochs(s, drift, onsetSeconds, cfg);
    catch err
        threw = strcmp(err.identifier, 'RELAX:secondPassTooManyFlaggedChannels');
        assert(contains(err.message, sprintf('%.1f', onsetSeconds(4))) && contains(err.message, '3 of 10') ...
            && contains(err.message, 's into the cleaned data'), ...
            'error message does not name the epoch and counts: %s', err.message);
    end
    assert(threw, '3 of 10 channels flagged did not raise RELAX:secondPassTooManyFlaggedChannels');

    % The limit is a setting
    cfg.SecondPassMaxProportionFlaggedChannels = 0.5;
    RELAX_flag_extreme_channel_epochs(s, drift, onsetSeconds, cfg);

    fprintf('test_RELAX_flag_extreme_channel_epochs: PASS\n');
end

function [s, drift] = quiet_scores()
    s.range = 50 * ones(10, 4);
    s.maxAbs = 30 * ones(10, 4);
    s.jpZ = zeros(10, 4);
    s.kurtZ = zeros(10, 4);
    s.shiftUpper = 100 * ones(10, 1);
    drift = -1 * ones(10, 4);
end
