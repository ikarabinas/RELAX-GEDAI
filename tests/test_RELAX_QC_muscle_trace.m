function test_RELAX_QC_muscle_trace()
% The smoothed 20–80 Hz power peaks inside the strongest muscle burst, and smoothing never
% crosses a boundary.

    rng(3);
    srate = 500; pnts = 10000;
    EEG = eeg_emptyset();
    EEG.srate = srate; EEG.nbchan = 4; EEG.pnts = pnts; EEG.trials = 1; EEG.xmin = 0;
    EEG.data = randn(4, pnts);
    EEG.data(:, pnts) = 1e4;                                   % lone final sample, its own 1-sample segment
    burst = 7001:7500;                                         % 14–15 s, segment 2
    EEG.data(2, burst) = EEG.data(2, burst) + 20 * sin(2 * pi * 50 * (burst / srate));
    edgeBurst = 4751:5000;                                     % last 0.5 s of segment 1
    EEG.data(3, edgeBurst) = EEG.data(3, edgeBurst) + 5 * sin(2 * pi * 50 * (edgeBurst / srate));
    EEG.event = struct('type', {'boundary', 'boundary'}, 'latency', {5000.5, pnts - 0.5}, 'duration', {0, 0});
    EEG = eeg_checkset(EEG);

    % A 0.5 s window is shorter than the 1 s burst, so the smoothed peak must fall inside it
    [trace, channelPower, filtered] = RELAX_QC_muscle_trace(EEG, 0.5);
    assert(numel(trace) == pnts && ~any(isnan(trace(1:pnts-1))), 'one finite value per sample outside the 1-sample segment');
    assert(isnan(trace(pnts)), 'a 1-sample segment has no trace');
    [~, peak] = max(trace);
    assert(peak >= burst(1) && peak <= burst(end), 'peak at sample %d is outside the burst', peak);
    assert(trace(5001) < trace(5000), 'segment 1''s edge burst must not smooth into segment 2');
    [~, loudest] = max(channelPower);
    assert(loudest == 2 && max(channelPower) < 1e3, 'channel 2 carries the strongest burst, and the lone huge sample does not enter channelPower');

    expected = RELAX_filtbutter(EEG, 20, 80, 4, 'bandpass', 'acausal', true);
    assert(isequal(size(filtered.data), size(EEG.data)) && isequal(filtered.data, expected.data), 'filtered is the 20-80 Hz band-passed EEG');

    fprintf('test_RELAX_QC_muscle_trace: PASS\n');
end
