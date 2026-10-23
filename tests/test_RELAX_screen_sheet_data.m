function test_RELAX_screen_sheet_data()
% Real fixture files: windows are 3 s, centres separated, missing windows carry notes.
    fixtures = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
    processed = fullfile(fixtures, 'RELAXProcessed');
    inputs = dir(fullfile(fixtures, '**', '*.set')); inputs = inputs(~contains({inputs.folder}, 'RELAXProcessed'));
    paths = RELAX_QC_file_paths(processed, inputs, "boundary_fixture");
    loaded = load(fullfile(processed, 'PerParticipantStats', 'RELAX_cfg_smoke.mat'), 'RELAX_cfg');
    opts = struct('windowSeconds', 3, 'nWindows', 3, 'blinkHalfSeconds', 0.5, ...
                  'topoHalfSeconds', 0.05, 'blinkElectrodeFallback', NaN);
    sheet = RELAX_screen_sheet_data(paths, loaded.RELAX_cfg, opts);
    assert(numel(sheet.muscle) == 3 && numel(sheet.blink) == 3 && iscell(sheet.interpolatedLabels));
    assert(numel(sheet.muscleSlope) == 3 && numel(sheet.muscleSlopeChannel) == 3 && numel(sheet.blinkResidual) == 3);
    assert(sheet.muscleSlopeThreshold == loaded.RELAX_cfg.MuscleSlopeThreshold);
    assert(isstruct(sheet.interpolatedReasons) && isequal(sheet.interpolatedReasons.labels, string(sheet.interpolatedLabels)), ...
           'the reasons are for the interpolated electrodes');
    shown = [sheet.muscle([sheet.muscle.available])];
    for w = shown, assert(abs(w.t(end) - w.t(1) - 3) < 2 / 500); end
    starts = arrayfun(@(w) w.t(1), shown);
    assert(all(diff(sort(starts)) >= 3 - 1e-9) || numel(starts) < 2, 'muscle windows do not overlap');
    for w = [sheet.muscle sheet.blink]
        assert(w.available || strlength(w.note) > 0, 'every missing window explains itself');
        assert(~w.available || strlength(sheet.baselineNote) > 0 || ~isempty(w.before));
    end
    unavailable = ~[sheet.muscle.available];
    assert(all(isnan(sheet.muscleSlope(unavailable))) && all(sheet.muscleSlopeChannel(unavailable) == ""));
    assert(all(isfinite(sheet.muscleSlope(~unavailable))) && all(strlength(sheet.muscleSlopeChannel(~unavailable)) > 0));
    for k = find(~unavailable)
        assert(ismember(sheet.muscleSlopeChannel(k), string(sheet.muscle(k).channels)), 'a muscle window shows the channel with the epoch''s highest slope');
        assert(numel(sheet.muscle(k).channels) <= 6);
    end
    assert(all(isnan(sheet.blinkResidual(~[sheet.blink.available]))));
    opts.nWindows = 30;
    crowded = RELAX_screen_sheet_data(paths, loaded.RELAX_cfg, opts);
    missingBlinks = crowded.blink(~[crowded.blink.available]);
    assert(~isempty(missingBlinks) && all(startsWith([missingBlinks.note], "Only")), 'too few blinks: notes say how many');
    assert(all(contains([missingBlinks.note], 's apart')) || crowded.blinkN == nnz([crowded.blink.available]));
    missingMuscle = crowded.muscle(~[crowded.muscle.available]);
    assert(~isempty(missingMuscle) && all(contains([missingMuscle.note], 'muscle epoch')) && all(contains([missingMuscle.note], 'at least 3 s apart')));
    opts.nWindows = 3;
    lowPassedCfg = loaded.RELAX_cfg; lowPassedCfg.LowPassFilter = 60;
    withoutSlopes = RELAX_screen_sheet_data(paths, lowPassedCfg, opts);
    assert(~any([withoutSlopes.muscle.available]) && all(contains([withoutSlopes.muscle.note], "60 Hz")), 'without slopes every muscle window says why');
    assert(all(isnan(withoutSlopes.muscleSlope)) && all(withoutSlopes.muscleSlopeChannel == "") && any([withoutSlopes.blink.available]), 'blink windows are unaffected');
    refusedCfg = loaded.RELAX_cfg; refusedCfg.Run_GEDAI = 1; refusedCfg.GEDAI_BroadbandOnly = false;
    assertConfigRefused(refusedCfg);
    sheetWithoutBaseline = RELAX_screen_sheet_data(paths, refusedCfg, opts);
    assertNoBefore(sheetWithoutBaseline);
    assert(all(isnan(sheetWithoutBaseline.blinkAverage.erpBaseline)), 'no baseline: blink average has no baseline trace');

    missingExtremes = paths; missingExtremes.extremes = fullfile(tempname, 'nope_Extremes_Rejected.set');
    assertNoBefore(RELAX_screen_sheet_data(missingExtremes, loaded.RELAX_cfg, opts));
    fprintf('sheet fixture: %d muscle windows, %d blink windows, blinkN %d, baselineNote "%s"\n', ...
        nnz([sheet.muscle.available]), nnz([sheet.blink.available]), sheet.blinkN, sheet.baselineNote);
    test_muscle_windows_stay_off_segment_edges(loaded.RELAX_cfg);
    fprintf('test_RELAX_screen_sheet_data: PASS\n');
end

function test_muscle_windows_stay_off_segment_edges(loadedCfg)
% The window shown for a burst near the start of a segment slides to stay inside that segment, and
% bursts 3 s apart that overlap once slid inside their segment are not both chosen.
    rng(7);
    srate = 500; seconds = 20; nChan = 32;
    chanlocs = make_qc_chanlocs(nChan);
    eeg = syntheticEeg(srate, seconds, chanlocs);
    eeg = addBurst(eeg, 11.2, 20);
    folder = tempname; mkdir(folder); cleanup = onCleanup(@() rmdir(folder, 's'));
    pop_saveset(eeg, 'filename', 'synthetic_RELAX.set', 'filepath', folder);
    paths = struct('cleaned', fullfile(folder, 'synthetic_RELAX.set'), 'extremes', fullfile(folder, 'absent.set'));
    opts = struct('windowSeconds', 3, 'nWindows', 3, 'blinkHalfSeconds', 0.5, ...
                  'topoHalfSeconds', 0.05, 'blinkElectrodeFallback', NaN);

    sheet = RELAX_screen_sheet_data(paths, loadedCfg, opts);
    strongest = sheet.muscle(1);
    assert(strongest.available && abs(mean(strongest.t) - 11.2) < 0.5, 'the strongest window is the burst');
    assert(strongest.t(1) >= 10 - 1 / srate && isempty(strongest.boundaries), 'the burst window starts inside its segment and crosses no boundary');

    crowded = syntheticEeg(srate, seconds, chanlocs);
    crowded = addBurst(crowded, 1.2, 20);
    crowded = addBurst(crowded, 4.2, 18);
    pop_saveset(crowded, 'filename', 'crowded_RELAX.set', 'filepath', folder);
    twoBursts = RELAX_screen_sheet_data(struct('cleaned', fullfile(folder, 'crowded_RELAX.set'), 'extremes', paths.extremes), loadedCfg, opts);
    shownCentres = arrayfun(@(w) mean(w.t), twoBursts.muscle([twoBursts.muscle.available]));
    assert(numel(shownCentres) >= 2 && min(pdist2Seconds(shownCentres)) >= 3 - 2 / srate, ...
           'bursts 3 s apart that overlap once slid inside the segment are not both chosen');

    blinkLabelsPresent = intersect(RELAX_screen_blink_labels(loadedCfg, NaN), {chanlocs.labels}, 'stable');
    assert(~isempty(blinkLabelsPresent), 'the synthetic montage must contain a blink electrode');
    eeg.RELAX.ListOfChannelsAfterRejections = setdiff({chanlocs.labels}, blinkLabelsPresent, 'stable');
    pop_saveset(eeg, 'filename', 'interpolated_RELAX.set', 'filepath', folder);
    paths.cleaned = fullfile(folder, 'interpolated_RELAX.set');
    interpolated = RELAX_screen_sheet_data(paths, loadedCfg, opts);
    assert(isequal(sort(interpolated.interpolatedLabels(:)), sort(blinkLabelsPresent(:))) && ~any([interpolated.blink.available]), ...
           'blinks are not ranked on interpolated electrodes');
    assert(all(contains([interpolated.blink.note], "interpolated")));
end

function eeg = syntheticEeg(srate, seconds, chanlocs)
% Two segments of seconds / 2 each, 1/f noise, every channel kept.
    eeg = eeg_emptyset(); eeg.srate = srate; eeg.nbchan = numel(chanlocs); eeg.pnts = seconds * srate;
    eeg.trials = 1; eeg.xmin = 0; eeg.chanlocs = chanlocs;
    eeg.data = make_pink_noise(eeg.nbchan, eeg.pnts);
    eeg.event = struct('type', 'boundary', 'latency', eeg.pnts / 2 + 0.5, 'duration', 0);
    eeg = eeg_checkset(eeg);
    eeg.allchan = chanlocs;
    eeg.RELAX.Data_has_been_cleaned = 1;
    eeg.RELAX.ListOfChannelsAfterRejections = {chanlocs.labels};
end

function eeg = addBurst(eeg, centreSeconds, amplitude)
    span = max(1, round((centreSeconds - 1.5) * eeg.srate) + 1) : round((centreSeconds + 1.5) * eeg.srate);
    seconds = (span - 1) / eeg.srate;
    envelope = exp(-((seconds - centreSeconds) / 0.6).^2);
    eeg.data(:, span) = eeg.data(:, span) + amplitude * envelope .* randn(size(eeg.data(:, span)));
end

function separations = pdist2Seconds(centres)
    differences = abs(centres(:) - centres(:)');
    separations = differences(triu(true(numel(centres)), 1));
end

function assertConfigRefused(cfg)
    refused = false;
    try
        RELAX_QC_check_config(cfg);
    catch
        refused = true;
    end
    assert(refused, 'the test config must be one RELAX_QC_check_config refuses');
end

function assertNoBefore(sheet)
    assert(strlength(sheet.baselineNote) > 0, 'the reason for the missing baseline is reported');
    available = [sheet.muscle([sheet.muscle.available]) sheet.blink([sheet.blink.available])];
    assert(~isempty(available) && all(arrayfun(@(w) isempty(w.before), available)));
end
