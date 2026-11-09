function test_RELAX_epoch_drift_slopes()
% The log-frequency log-power slope (1-75 Hz) of each channel in each epoch, which the
% extreme pass compares against ExtremeDriftSlopeThreshold: about 0 for white noise, about
% -1 for pink noise, and steep for a channel carrying only slow drift.

    rng(1);
    srate = 250;
    nEpochs = 20;
    t = (0:srate - 1) / srate;
    data = zeros(3, srate, nEpochs);
    for e = 1:nEpochs
        data(1, :, e) = randn(1, srate);
        data(2, :, e) = make_pink_noise(1, srate);
        data(3, :, e) = 100 * sin(2 * pi * 0.3 * t + e) + 0.001 * randn(1, srate);
    end
    EEG = eeg_emptyset();
    EEG.data = data;
    EEG.nbchan = 3; EEG.pnts = srate; EEG.trials = nEpochs; EEG.srate = srate;
    EEG.xmin = 0; EEG.xmax = (srate - 1) / srate;
    EEG.chanlocs = make_qc_chanlocs(3);
    EEG = eeg_checkset(EEG);

    slopes = RELAX_epoch_drift_slopes(EEG);
    assert(isequal(size(slopes), [3 nEpochs]), 'slopes should be channels x epochs, got %s', mat2str(size(slopes)));
    assert(all(abs(slopes(1, :)) < 0.5), 'white noise slope is not near 0: %s', mat2str(slopes(1, :), 2));
    assert(all(slopes(2, :) < -0.5 & slopes(2, :) > -2), 'pink noise slope is not near -1: %s', mat2str(slopes(2, :), 2));
    assert(all(slopes(3, :) < -4), 'a slow-drift channel does not fall below the -4 drift threshold: %s', mat2str(slopes(3, :), 2));
    fprintf('test_RELAX_epoch_drift_slopes: PASS\n');
end
