function test_RELAX_screen_muscle_slopes()
% A white-noise burst on one real channel is the epoch with the highest slope, found at the burst's time.
    rng(11);
    srate = 500; seconds = 40; nChan = 20;
    burstChannel = 'E7'; burstStart = 20;
    eeg = pinkNoiseEeg(srate, seconds, nChan);
    burstSpan = round(burstStart * srate) + 1 : round((burstStart + 1) * srate);
    burstRow = find(strcmp({eeg.chanlocs.labels}, burstChannel));
    eeg.data(burstRow, burstSpan) = eeg.data(burstRow, burstSpan) + 10 * std(eeg.data(burstRow, :)) * randn(1, numel(burstSpan));
    cfg = struct('MuscleSlopeThreshold', -0.31, 'LowPassFilter', 80);

    slopes = RELAX_screen_muscle_slopes(eeg, cfg);

    assert(~any(slopes.labels == "E20"), 'a channel the pipeline interpolated is not measured');
    assert(numel(slopes.labels) == nChan - 1 && size(slopes.values, 1) == numel(slopes.labels));
    assert(size(slopes.values, 2) == numel(slopes.startSeconds), 'one start time per epoch');
    assert(slopes.threshold == -0.31 && slopes.note == "");
    [highest, epoch] = max(max(slopes.values, [], 1));
    [~, channel] = max(slopes.values(:, epoch));
    assert(slopes.labels(channel) == burstChannel, 'the highest slope is on the channel with the burst');
    assert(abs(slopes.startSeconds(epoch) - burstStart) <= 0.5, 'the highest slope is at the burst');
    assert(highest > -0.31, 'the burst is above the muscle threshold');

    unavailable = RELAX_screen_muscle_slopes(eeg, struct('MuscleSlopeThreshold', -0.31, 'LowPassFilter', 60));
    assert(isempty(unavailable.values) && strlength(unavailable.note) > 0, 'a 60 Hz low-pass cannot yield 7-75 Hz slopes');
    fprintf('test_RELAX_screen_muscle_slopes: PASS\n');
end

function eeg = pinkNoiseEeg(srate, seconds, nChan)
    chanlocs = make_qc_chanlocs(nChan);
    eeg = eeg_emptyset(); eeg.srate = srate; eeg.nbchan = nChan; eeg.pnts = seconds * srate;
    eeg.trials = 1; eeg.xmin = 0; eeg.chanlocs = chanlocs;
    eeg.data = make_pink_noise(nChan, eeg.pnts);
    eeg.data(nChan, :) = NaN;
    eeg = eeg_checkset(eeg);
    eeg.allchan = chanlocs;
    eeg.RELAX.Data_has_been_cleaned = 1;
    eeg.RELAX.ListOfChannelsAfterRejections = {chanlocs(1:nChan - 1).labels};
end
