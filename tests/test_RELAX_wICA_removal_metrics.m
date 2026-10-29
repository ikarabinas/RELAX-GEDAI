function test_RELAX_wICA_removal_metrics()
% Removed-power fractions from wICA input/output pairs, split by muscle epochs.

    rng(7);
    srate = 500;
    nChan = 4;
    nSamples = 20 * srate;
    % A slope fitted to one 1 s periodogram has an sd of about 0.25, so the two classes sit far
    % either side of -0.31: brain power falls as 1/f^2 (slope -2), muscle rises as f^2 (slope +2)
    bin = 0:nSamples - 1;
    cyclesPerRecording = max(1, min(bin, nSamples - bin));
    brain = real(ifft(fft(randn(nChan, nSamples), [], 2) ./ cyclesPerRecording, [], 2));
    brain = 10 * brain ./ std(brain, 0, 2);
    muscle = zeros(nChan, nSamples);
    muscle(:, nSamples/2+1:end) = 50 * diff(randn(nChan, nSamples/2 + 1), 1, 2);

    EEGbefore = struct('srate', srate, 'nbchan', nChan, 'pnts', nSamples, ...
        'data', single(brain + muscle), ...
        'chanlocs', struct('labels', {'E1', 'E2', 'E3', 'E4'}), ...
        'event', struct('type', {'boundary'}, 'latency', {5250.5}));
    cfg = struct('MuscleSlopeThreshold', -0.31);

    % A 10.5 s and a 9.5 s segment give 10 + 9 whole 1 s epochs; the sub-second slivers are dropped
    [metrics, perChannel] = RELAX_wICA_removal_metrics(EEGbefore, EEGbefore, cfg);
    assert(metrics.NumberOfEpochs == 19, 'epochs should tile each segment, not span the boundary');
    assert(abs(metrics.ProportionChannelEpochsMuscle - 9/19) < 1e-12, 'the noisy second segment should be muscle');
    assert(metrics.ProportionEpochsBlink == 0 && isnan(metrics.RemovedFraction_Blink_Delta), 'no blink mask, no blink epochs');
    assert(metrics.RemovedFraction_Clean_Alpha == 0 && metrics.RemovedFraction_Muscle_Broadband == 0, ...
        'nothing removed should give zero fractions');
    assert(isequal(size(perChannel.removedPowerClean), [nChan numel(perChannel.bandNames)]), ...
        'per-channel output should be channels x bands');

    EEGzero = EEGbefore;
    EEGzero.data(:) = 0;
    metrics = RELAX_wICA_removal_metrics(EEGbefore, EEGzero, cfg);
    assert(abs(metrics.RemovedFraction_Clean_Beta - 1) < 1e-6 && abs(metrics.RemovedFraction_Muscle_Gamma - 1) < 1e-6, ...
        'removing everything should give fractions of one');

    EEGbrainOnly = EEGbefore;
    EEGbrainOnly.data = single(brain);
    metrics = RELAX_wICA_removal_metrics(EEGbefore, EEGbrainOnly, cfg);
    assert(metrics.RemovedFraction_Clean_Broadband < 1e-6, 'removing only muscle should leave clean epochs untouched');
    assert(metrics.RemovedFraction_Muscle_HighGamma > 0.9, 'the muscle epochs'' high gamma is almost all muscle');

    EEGnoMuscle = EEGbefore;
    EEGnoMuscle.data = single(brain);
    metrics = RELAX_wICA_removal_metrics(EEGnoMuscle, EEGnoMuscle, cfg);
    assert(isnan(metrics.RemovedFraction_Muscle_Alpha), 'no muscle epochs should give NaN muscle fractions');

    % Epochs touching a blink leave the clean set whatever their slope
    blinkMask = false(1, nSamples);
    blinkMask(1.2*srate : 1.4*srate) = true;
    [metrics, perChannel] = RELAX_wICA_removal_metrics(EEGbefore, EEGbrainOnly, cfg, blinkMask);
    assert(abs(metrics.ProportionEpochsBlink - 1/19) < 1e-12, 'one epoch touches the blink');
    assert(all(perChannel.channelEpochsClean == 9) && all(perChannel.channelEpochsBlink == 1), ...
        'the blink epoch should move from clean to blink on every channel');

    % Muscle is judged per channel: one noisy channel does not make the others muscle
    oneNoisy = brain;
    oneNoisy(1, :) = oneNoisy(1, :) + 50 * diff(randn(1, nSamples + 1));
    EEGone = EEGbefore;
    EEGone.data = single(oneNoisy);
    [~, perChannel] = RELAX_wICA_removal_metrics(EEGone, EEGone, cfg);
    assert(perChannel.channelEpochsMuscle(1) == 19 && all(perChannel.channelEpochsClean(2:end) == 19), ...
        'only the noisy channel should count as muscle');

    fprintf('test_RELAX_wICA_removal_metrics: PASS\n');
end
