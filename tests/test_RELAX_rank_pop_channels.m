function test_RELAX_rank_pop_channels()
% Re-ranking bad contacts after eye-movement flagging: a channel whose "steps" are all
% saccades must not be rejected for its step rate. The sample-level instability criterion
% this once also checked was retired 2026-09-22.

    srate = 1000; n60 = 60 * srate; rng(11);
    function EEG = mkEEG(data)
        EEG = struct('data', data, 'nbchan', size(data,1), 'pnts', size(data,2), ...
                     'srate', srate, 'trials', 1, 'event', struct('type', {}, 'latency', {}));
        EEG.chanlocs = struct('labels', arrayfun(@(k) sprintf('E%d', k), 1:size(data,1), 'UniformOutput', false));
    end

    % E1: telegraph every 2 s (29 toggles/min); E2: quiet
    x = randn(1, n60) * 5; lvl = 0;
    for a = 2001:2000:n60; lvl = 300 - lvl; x(a:end) = x(a:end) + (lvl - (300 - lvl)); end
    [steps, info] = RELAX_detect_steps(mkEEG([x; randn(1, n60) * 5]));
    assert(isequal(info.badChannels, {'E1'}), 'precondition: detector rejects the telegraph channel');

    %% 1. counting the same steps the detector counted reproduces its ranking exactly
    info2 = RELAX_rank_pop_channels(steps, info, ~steps.commonMode);
    assert(isequal(info2.badChannels, info.badChannels) && isequal(info2.badChannelReason, info.badChannelReason));
    assert(isequal(info2.stepRatePerMinute, info.stepRatePerMinute));
    assert(isequal(info2.transitionFraction, info.transitionFraction));

    %% 2. when every step on E1 is an eye movement, E1 is not a bad contact
    steps.eyeMovement = steps.channel == 1;
    info3 = RELAX_rank_pop_channels(steps, info, ~steps.commonMode & ~steps.eyeMovement);
    assert(isempty(info3.badChannels), 'eye-movement steps must not count against the channel');
    assert(iscell(info3.badChannels) && isequal(size(info3.badChannels), [0 1]));
    assert(iscell(info3.badChannelReason) && isequal(size(info3.badChannelReason), [0 1]));
    assert(info3.stepRatePerMinute(1) == 0 && info3.transitionFraction(1) == 0);
    assert(isequal(size(info3.stepRatePerMinute), [2 1]));

    %% 3. the default keep mask is ~commonMode & ~eyeMovement
    info4 = RELAX_rank_pop_channels(steps, info);
    assert(isequal(info4.badChannels, info3.badChannels) && isequal(info4.stepRatePerMinute, info3.stepRatePerMinute));

    %% 4. an empty keep mask leaves the other info fields untouched
    assert(info3.nSteps == info.nSteps);

    %% 5. every criterion now answers to the keep mask: retiring the sample-level instability
    %% criterion (2026-09-22) left nothing that rejects a channel whose steps were all excluded
    tt = 0:n60-1; sq = 300 * (mod(floor(tt / 250), 2) == 0);
    [steps5, info5] = RELAX_detect_steps(mkEEG(randn(1, n60) * 5 + sq));
    assert(isequal(info5.badChannels, {'E1'}) && contains(info5.badChannelReason{1}, 'steps/min'), ...
        'the fixture must reject on rate while its steps count: %s', strjoin(info5.badChannelReason, '; '));
    steps5.eyeMovement = true(size(steps5.sample));
    info6 = RELAX_rank_pop_channels(steps5, info5);
    assert(isempty(info6.badChannels), ...
        'a channel whose every step is an eye movement must survive even at a telegraph rate');
    assert(info6.stepRatePerMinute(1) == 0 && info6.transitionFraction(1) == 0);

    %% 6. worst-first ordering survives re-ranking
    rates = [108 72 36]; X = randn(3, n60) * 5;
    for c = 1:3
        at = round(linspace(2000, n60 - 2000, rates(c)));
        for j = 1:numel(at); X(c, at(j):end) = X(c, at(j):end) + 300 * (-1)^j; end
    end
    [steps7, info7] = RELAX_detect_steps(mkEEG(X([2 1 3], :)));
    steps7.eyeMovement = false(size(steps7.sample));
    info8 = RELAX_rank_pop_channels(steps7, info7);
    assert(isequal(info8.badChannels, {'E2'; 'E1'; 'E3'}));

    %% 7. degenerate (flat / railed) channels are rejected first, reason 'flat'
    X = [zeros(1, n60); x; randn(1, n60) * 5];        % E1 flat, E2 telegraph (2 s), E3 quiet
    [steps9, info9] = RELAX_detect_steps(mkEEG(X));
    assert(isequal(info9.degenerateChannels, {'E1'}));
    assert(isequal(info9.badChannels, {'E1'; 'E2'}), 'flat channel first, then the telegraph channel (got %s)', strjoin(info9.badChannels', ' '));
    assert(contains(info9.badChannelReason{1}, 'flat') && contains(info9.badChannelReason{2}, 'steps/min'));
    steps9.eyeMovement = false(size(steps9.sample));
    info10 = RELAX_rank_pop_channels(steps9, info9);
    assert(isequal(info10.badChannels, {'E1'; 'E2'}) && contains(info10.badChannelReason{1}, 'flat'), 're-ranking keeps the flat channel first');

    fprintf('test_RELAX_rank_pop_channels: PASS\n');
end
