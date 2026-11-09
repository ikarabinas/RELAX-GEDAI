function test_RELAX_second_pass_channel_epochs()
% The post-cleaning stage: off, it changes nothing and records only the setting; on, it
% epochs the cleaned data as the first pass did, repairs the channel-epochs that are
% extreme by the first pass's standard, never across a boundary, and counts what it did.
%
% Fixture: 40 s of average-referenced pink noise on 32 channels with one boundary between
% samples 5000 and 5001. The first-pass statistics are fitted on the same data before any
% spike is added.

    [EEG, cfg] = fixture();
    locked = lock_from(EEG, cfg);

    % Off: absent, then 0
    for setting = {[], 0}
        offCfg = cfg;
        if ~isempty(setting{1}); offCfg.InterpolateExtremeChannelEpochsAfterCleaning = setting{1}; end
        [out, cov] = RELAX_second_pass_channel_epochs(EEG, [], offCfg);
        assert(isequal(out.data, EEG.data), 'the stage changed data while off');
        assert(cov.InterpolateExtremeChannelEpochsAfterCleaning == 0, 'the off setting was not recorded as 0');
        assert(isnan(cov.SecondPassChannelEpochsInterpolated) && isnan(cov.SecondPassEpochsEvaluated) ...
            && isnan(cov.SecondPassFlaggedByKurtosis) && isnan(cov.SecondPassChannelEpochsNotRepaired), ...
            'counts should be NaN while off');
    end

    cfg.InterpolateExtremeChannelEpochsAfterCleaning = 1;

    % spike
    spiked = EEG;
    spikeSamples = 12 * EEG.srate + (1:20);
    spiked.data(7, spikeSamples) = spiked.data(7, spikeSamples) + 400;
    [out, cov] = RELAX_second_pass_channel_epochs(spiked, locked, cfg);
    assert(max(abs(out.data(7, spikeSamples))) < 100, 'the spike survived: peak %.0f uV', max(abs(out.data(7, spikeSamples))));
    differs = out.data ~= spiked.data;
    assert(~any(differs([1:6 8:32], :), 'all'), 'a channel other than the spiked one changed');
    assert(cov.SecondPassChannelEpochsInterpolated >= 1 && cov.SecondPassFlaggedByKurtosis >= 1, ...
        'the spike was not flagged by kurtosis and interpolated');
    epochs = cov.SecondPassEpochsEvaluated;
    assert(epochs > 50, 'expected the 40 s fixture to give over 50 epochs, got %d', epochs);
    proportion = cov.SecondPassProportionChannelTimeInterpolated;
    assert(proportion >= nnz(differs) / (32 * EEG.pnts), 'proportion interpolated is below the share of samples that changed');
    assert(proportion <= cov.SecondPassChannelEpochsInterpolated * EEG.srate / (32 * EEG.pnts), ...
        'proportion interpolated exceeds what the interpolated channel-epochs could cover');
    assert(cov.SecondPassEpochsWithInterpolatedChannels >= 1 ...
        && cov.SecondPassEpochsWithInterpolatedChannels <= cov.SecondPassChannelEpochsInterpolated, ...
        'epochs with interpolation inconsistent with channel-epochs interpolated');

    % boundary: a spike just before the boundary is repaired without reaching across it
    spiked = EEG;
    spiked.data(7, 4980:4995) = spiked.data(7, 4980:4995) + 400;
    [out, cov] = RELAX_second_pass_channel_epochs(spiked, locked, cfg);
    changed = out.data ~= spiked.data;
    assert(any(changed(7, 4980:4995)), 'the spike before the boundary was not repaired');
    assert(~any(changed(:, 5000) & changed(:, 5001)), 'a repair spans the boundary');
    assert(~any(changed(:, 5001:end), 'all'), 'the repair before the boundary reached past it');

    % label_subset: a channel removed after the first pass
    subset = pop_select(EEG, 'nochannel', {'E3'});
    subset.data(6, spikeSamples) = subset.data(6, spikeSamples) + 400;  % E7 is now row 6
    [out, cov] = RELAX_second_pass_channel_epochs(subset, locked, cfg);
    assert(size(out.data, 1) == 31, 'the stage changed the channel count');
    assert(max(abs(out.data(6, spikeSamples))) < 100, 'with E3 removed, the spike on E7 was not repaired');
    assert(cov.SecondPassEpochsEvaluated == epochs, 'epoch count changed with a channel removed');

    % blink: a blink-sized deflection is held to the plain first-pass shift bound even where
    % the cleaned data carry blink events (no blink exemption)
    spiked = EEG;
    blinkSamples = 30 * EEG.srate + (1:100);
    spiked.data(7, blinkSamples) = spiked.data(7, blinkSamples) + 600 * sin(pi * (1:100) / 100);
    spiked.event(end + 1).type = 'EyeBlinkMax';
    spiked.event(end).latency = 30 * EEG.srate + 50;
    spiked.event(end).duration = 0;
    spiked = eeg_checkset(spiked, 'eventconsistency');
    [out, cov] = RELAX_second_pass_channel_epochs(spiked, locked, cfg);
    assert(cov.SecondPassFlaggedByShiftHigh >= 1, 'the blink-sized shift was not flagged by the shift bound');
    assert(max(abs(out.data(7, blinkSamples))) < 150, 'the blink-sized shift survived: peak %.0f uV', ...
        max(abs(out.data(7, blinkSamples))));

    % plain_epoching: with RespectBoundaries = 0 the first pass used the fixed 0.5 s grid,
    % which drops the epoch straddling a boundary at 20.4 s and so never sees the 0.4 s
    % before it; the second pass must epoch the same way, while boundary-aware tiling
    % covers that stretch
    offGrid = EEG;
    offGrid.event(1).latency = 5100.5;
    offGrid.data(7, 5060:5075) = offGrid.data(7, 5060:5075) + 400;
    out = RELAX_second_pass_channel_epochs(offGrid, locked, cfg);
    assert(any(out.data(7, 5060:5075) ~= offGrid.data(7, 5060:5075)), ...
        'boundary-aware epoching did not repair the spike just before the boundary');
    plainCfg = cfg;
    plainCfg.RespectBoundaries = 0;
    [~, plainEpoched] = RELAX_excluding_channels_and_epoching(offGrid, plainCfg);
    [out, cov] = RELAX_second_pass_channel_epochs(offGrid, locked, plainCfg);
    assert(cov.SecondPassEpochsEvaluated == plainEpoched.trials, ...
        'with RespectBoundaries = 0 the second pass evaluated %d epochs, the fixed-grid epoching %d', ...
        cov.SecondPassEpochsEvaluated, plainEpoched.trials);
    assert(isequal(out.data(7, 5060:5075), offGrid.data(7, 5060:5075)), ...
        'with RespectBoundaries = 0 the second pass repaired data the fixed grid never epochs');

    % no_coordinates: a flagged channel that cannot be interpolated is counted apart
    noLocs = EEG;
    noLocs.data(7, spikeSamples) = noLocs.data(7, spikeSamples) + 400;
    [noLocs.chanlocs(7).theta, noLocs.chanlocs(7).radius, noLocs.chanlocs(7).X, noLocs.chanlocs(7).Y, noLocs.chanlocs(7).Z] = deal([]);
    [out, cov] = RELAX_second_pass_channel_epochs(noLocs, locked, cfg);
    assert(isequal(out.data, noLocs.data), 'data changed although the only flagged channel has no location');
    assert(cov.SecondPassChannelEpochsInterpolated == 0 && cov.SecondPassEpochsWithInterpolatedChannels == 0 ...
        && cov.SecondPassProportionChannelTimeInterpolated == 0, 'an unrepaired channel-epoch was counted as interpolated');
    assert(cov.SecondPassChannelEpochsNotRepaired >= 1, 'the flagged channel without a location was not counted as not repaired');

    fprintf('test_RELAX_second_pass_channel_epochs: PASS\n');
end

function [EEG, cfg] = fixture()
    rng(1);
    srate = 250;
    nSamples = 40 * srate;
    EEG = eeg_emptyset();
    data = 10 * make_pink_noise(32, nSamples);
    EEG.data = single(data - mean(data, 1));
    EEG.srate = srate; EEG.nbchan = 32; EEG.pnts = nSamples; EEG.trials = 1;
    EEG.xmin = 0; EEG.xmax = (nSamples - 1) / srate;
    EEG.chanlocs = make_qc_chanlocs(32);
    EEG.allchan = EEG.chanlocs;
    EEG.event = struct('type', 'boundary', 'latency', 5000.5, 'duration', 0);
    EEG = eeg_checkset(EEG, 'eventconsistency');

    cfg = struct();
    cfg.ms_per_sample = 1000 / srate;
    cfg.ExtremeVoltageShiftThreshold = 25;
    cfg.ExtremeAbsoluteVoltageThreshold = 1000;
    cfg.ExtremeImprobableVoltageDistributionThreshold = 10;
    cfg.ExtremeSingleChannelKurtosisThreshold = 10;
    cfg.ExtremeAllChannelKurtosisThreshold = 10;
    cfg.ExtremeDriftSlopeThreshold = -4;
    cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel = 1;
    cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0;
    cfg.ProbabilityDataHasNoBlinks = 2;
    cfg.RespectBoundaries = 1;
end

function locked = lock_from(EEG, cfg)
    [~, epoched] = RELAX_excluding_channels_and_epoching_boundary(EEG, cfg);
    data = double(epoched.data);
    locked = RELAX_lock_extreme_thresholds(data - mean(data, 1), {epoched.chanlocs.labels}, cfg);
end
