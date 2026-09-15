function test_RELAX_detect_steps()
% Step detection: recovers real steps, rejects ramps and noise, degrades safely.

    srate = 1000;
    n     = 10000;          % 10 s
    rng(7);

    % ---- helper: build a 2-channel EEG with no boundaries ----
    function EEG = mkEEG(data)
        EEG = struct();
        EEG.data     = data;
        EEG.nbchan   = size(data,1);
        EEG.pnts     = size(data,2);
        EEG.srate    = srate;
        EEG.trials   = 1;
        EEG.event    = struct('type', {}, 'latency', {});
        EEG.chanlocs = struct('labels', ...
            arrayfun(@(k) sprintf('E%d', k), 1:size(data,1), 'UniformOutput', false));
    end

    % ---- 1. a known step is found at the right sample with the right amplitude ----
    x = randn(1, n) * 5;                 % 5 uV rms background
    x(4001:end) = x(4001:end) + 500;     % +500 uV step, new level begins at 4001
    EEG = mkEEG([x; randn(1, n) * 5]);

    [steps, info] = RELAX_detect_steps(EEG);

    assert(info.available, 'detector should be available on ordinary data: %s', info.status);
    onCh1 = steps.sample(steps.channel == 1);
    assert(numel(onCh1) == 1, 'expected exactly 1 step on channel 1, got %d', numel(onCh1));
    assert(abs(onCh1 - 4001) <= 2, 'step sample should be ~4001, got %d', onCh1);
    amp = steps.amplitude(steps.channel == 1);
    assert(abs(amp - 500) < 25, 'step amplitude should be ~+500 uV, got %g', amp);
    assert(steps.z(steps.channel == 1) > 20, 'z should clear the threshold');
    assert(isempty(steps.sample(steps.channel == 2)), 'clean channel 2 must yield no steps');
    assert(strcmp(steps.label{steps.channel == 1}, 'E1'), 'label should accompany the index');

    % sign is preserved for a downward step
    x2 = randn(1, n) * 5;
    x2(4001:end) = x2(4001:end) - 500;
    [steps2, ~] = RELAX_detect_steps(mkEEG(x2));
    assert(steps2.amplitude(1) < 0, 'a downward step must have negative amplitude');

    % ---- 2. NO detection on a smooth ramp (the exploratory detector's failure mode) ----
    ramp = randn(1, n) * 5 + linspace(0, 500, n);   % same 500 uV excursion, spread out
    [stepsRamp, ~] = RELAX_detect_steps(mkEEG(ramp));
    assert(isempty(stepsRamp.sample), ...
        'a smooth 500 uV ramp must not be detected as a step, got %d detections', ...
        numel(stepsRamp.sample));

    % ---- 3. no detection on noise alone ----
    [stepsNoise, ~] = RELAX_detect_steps(mkEEG(randn(2, n) * 20));
    assert(isempty(stepsNoise.sample), ...
        'plain noise must not trigger detections, got %d', numel(stepsNoise.sample));

    % ---- 4. refractory: two steps 5 ms apart collapse to one ----
    x3 = randn(1, n) * 5;
    x3(4001:end) = x3(4001:end) + 400;
    x3(4006:end) = x3(4006:end) + 400;
    [steps3, ~] = RELAX_detect_steps(mkEEG(x3), 'RefractoryMs', 20);
    assert(numel(steps3.sample) == 1, ...
        'two steps inside the refractory window must give 1 detection, got %d', numel(steps3.sample));

    % ...but 500 ms apart they stay separate
    x4 = randn(1, n) * 5;
    x4(3001:end) = x4(3001:end) + 400;
    x4(3501:end) = x4(3501:end) - 400;
    [steps4, ~] = RELAX_detect_steps(mkEEG(x4), 'RefractoryMs', 20);
    assert(numel(steps4.sample) == 2, ...
        'two well-separated steps must give 2 detections, got %d', numel(steps4.sample));

    % ---- 5. flat channel: no steps, no NaN, listed as degenerate ----
    flat = zeros(1, n);
    [steps5, info5] = RELAX_detect_steps(mkEEG([flat; randn(1,n)*5]));
    assert(isempty(steps5.sample(steps5.channel == 1)), 'a flat channel must yield no steps');
    assert(ismember('E1', info5.degenerateChannels), 'flat channel must be listed as degenerate');
    assert(size(info5.degenerateChannels,1) >= 1 && size(info5.degenerateChannels,2) == 1, ...
        'info.degenerateChannels must be Nx1');
    assert(~any(isnan(steps5.z)), 'no NaN z-scores');

    % ---- 6. absolute amplitude floor suppresses a tiny step on a very quiet channel ----
    quiet = randn(1, n) * 0.01;
    quiet(4001:end) = quiet(4001:end) + 2;    % 2 uV: huge in z, below a 10 uV floor
    [steps6, ~] = RELAX_detect_steps(mkEEG(quiet), 'MinAmplitudeMicrovolts', 10);
    assert(isempty(steps6.sample), 'a 2 uV step must be suppressed by a 10 uV floor');
    [steps6b, ~] = RELAX_detect_steps(mkEEG(quiet), 'MinAmplitudeMicrovolts', 1);
    assert(numel(steps6b.sample) == 1, 'the same step must survive a 1 uV floor');

    % ---- 7. boundaries: no detection manufactured at the discontinuity itself ----
    xb = randn(1, n) * 5;
    xb(5001:end) = xb(5001:end) + 800;     % a big jump that IS the boundary, not a pop
    EEGb = mkEEG(xb);
    EEGb.event = struct('type', {'boundary'}, 'latency', {5000.5});
    [stepsB, ~] = RELAX_detect_steps(EEGb);
    assert(isempty(stepsB.sample), ...
        'a level change coinciding with a boundary is a splice, not a pop; got %d detections', ...
        numel(stepsB.sample));

    % a genuine pop well inside the second segment is still found
    xb2 = randn(1, n) * 5;
    xb2(5001:end) = xb2(5001:end) + 800;
    xb2(7001:end) = xb2(7001:end) + 400;
    EEGb2 = mkEEG(xb2);
    EEGb2.event = struct('type', {'boundary'}, 'latency', {5000.5});
    [stepsB2, ~] = RELAX_detect_steps(EEGb2);
    assert(numel(stepsB2.sample) == 1 && abs(stepsB2.sample - 7001) <= 2, ...
        'the pop inside the segment must still be found');

    % ---- 8. Population A: a chronically popping channel is flagged for rejection ----
    xa = randn(1, n) * 5;
    for s = 500:500:9500                    % 19 steps in 10 s = 114/min
        xa(s:end) = xa(s:end) + 300 * (-1)^(s/500);
    end
    [~, infoA] = RELAX_detect_steps(mkEEG([xa; randn(1,n)*5]), ...
        'StepsPerMinuteToRejectChannel', 30);
    assert(ismember('E1', infoA.badChannels), 'chronically popping channel must be flagged');
    assert(~ismember('E2', infoA.badChannels), 'the clean channel must not be flagged');
    assert(infoA.stepRatePerMinute(1) > 30, 'rate should exceed the threshold');
    assert(numel(infoA.badChannelReason) == numel(infoA.badChannels), 'one reason per bad channel');
    assert(size(infoA.badChannels,1) >= 1 && size(infoA.badChannels,2) == 1, 'info.badChannels must be Nx1');
    assert(size(infoA.badChannelReason,1) >= 1 && size(infoA.badChannelReason,2) == 1, ...
        'info.badChannelReason must be Nx1');

    % ---- 8b. Population A ordering: badChannels must be worst-first for N>=2 ----
    % Three channels with distinct step rates (~198 / ~90 / ~60 steps/min). Channel-index
    % order is E1,E2,E3, but rate order is E2 (highest) > E1 > E3 (lowest) - deliberately
    % different, so a regression to channel-index order would fail this assertion.
    countHigh = 33;  % 33 steps in 10 s = 198 steps/min -> E2
    countMid  = 15;  % 15 steps in 10 s =  90 steps/min -> E1
    countLow  = 10;  % 10 steps in 10 s =  60 steps/min -> E3

    xOrdE1 = randn(1, n) * 5;
    posMid = round(linspace(300, 9700, countMid));
    for k = 1:numel(posMid)
        xOrdE1(posMid(k):end) = xOrdE1(posMid(k):end) + 300 * (-1)^k;
    end

    xOrdE2 = randn(1, n) * 5;
    posHigh = round(linspace(300, 9700, countHigh));
    for k = 1:numel(posHigh)
        xOrdE2(posHigh(k):end) = xOrdE2(posHigh(k):end) + 300 * (-1)^k;
    end

    xOrdE3 = randn(1, n) * 5;
    posLow = round(linspace(300, 9700, countLow));
    for k = 1:numel(posLow)
        xOrdE3(posLow(k):end) = xOrdE3(posLow(k):end) + 300 * (-1)^k;
    end

    [~, infoOrd] = RELAX_detect_steps(mkEEG([xOrdE1; xOrdE2; xOrdE3]), ...
        'StepsPerMinuteToRejectChannel', 30);
    assert(numel(infoOrd.badChannels) == 3, ...
        'expected all 3 channels flagged bad, got %d', numel(infoOrd.badChannels));
    assert(isequal(infoOrd.badChannels, {'E2'; 'E1'; 'E3'}), ...
        'badChannels must be worst-first by step rate (E2>E1>E3), got %s', ...
        strjoin(infoOrd.badChannels, ','));
    assert(size(infoOrd.badChannels,1) == 3 && size(infoOrd.badChannels,2) == 1, ...
        'badChannels must be Nx1');
    assert(size(infoOrd.badChannelReason,1) == 3 && size(infoOrd.badChannelReason,2) == 1, ...
        'badChannelReason must be Nx1');

    % ---- 9. empty result has the right shape, not a bare [] ----
    [stepsE, infoE] = RELAX_detect_steps(mkEEG(zeros(2, n) + randn(2, n)*1e-9));
    assert(isnumeric(stepsE.channel) && size(stepsE.channel,1)==0 && size(stepsE.channel,2)==1, ...
        'empty .channel must be exactly 0x1 numeric, not a bare []');
    assert(isnumeric(stepsE.sample) && size(stepsE.sample,1)==0 && size(stepsE.sample,2)==1, ...
        'empty .sample must be exactly 0x1 numeric, not a bare []');
    assert(isnumeric(stepsE.amplitude) && size(stepsE.amplitude,1)==0 && size(stepsE.amplitude,2)==1, ...
        'empty .amplitude must be exactly 0x1 numeric, not a bare []');
    assert(isnumeric(stepsE.z) && size(stepsE.z,1)==0 && size(stepsE.z,2)==1, ...
        'empty .z must be exactly 0x1 numeric, not a bare []');
    assert(iscell(stepsE.label) && size(stepsE.label,1)==0 && size(stepsE.label,2)==1, ...
        'empty .label must be exactly 0x1 cell, not a bare {}');
    assert(infoE.nSteps == 0 && infoE.nChannelsAffected == 0, 'zero counts on an empty result');
    assert(iscell(infoE.badChannels) && size(infoE.badChannels,1)==0 && size(infoE.badChannels,2)==1, ...
        'empty info.badChannels must be exactly 0x1 cell, not a bare {}');
    assert(iscell(infoE.badChannelReason) && size(infoE.badChannelReason,1)==0 && size(infoE.badChannelReason,2)==1, ...
        'empty info.badChannelReason must be exactly 0x1 cell, not a bare {}');
    assert(iscell(infoE.degenerateChannels) && size(infoE.degenerateChannels,1)==0 && size(infoE.degenerateChannels,2)==1, ...
        'empty info.degenerateChannels must be exactly 0x1 cell, not a bare {}');

    % ---- 10. params are echoed so the saved struct is self-describing ----
    [~, infoP] = RELAX_detect_steps(mkEEG(randn(2,n)*5), 'ZThreshold', 33);
    assert(infoP.params.ZThreshold == 33, 'params must echo the effective settings');

    fprintf('test_RELAX_detect_steps: PASS\n');
end
