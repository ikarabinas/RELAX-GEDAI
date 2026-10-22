function test_RELAX_screen_sheet_data()
% Real fixture files: windows are 3 s, centres separated, missing windows carry notes.
    fixtures = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
    processed = fullfile(fixtures, 'RELAXProcessed');
    inputs = dir(fullfile(fixtures, '**', '*.set')); inputs = inputs(~contains({inputs.folder}, 'RELAXProcessed'));
    paths = RELAX_QC_file_paths(processed, inputs, "boundary_fixture");
    loaded = load(fullfile(processed, 'PerParticipantStats', 'RELAX_cfg_smoke.mat'), 'RELAX_cfg');
    opts = struct('windowSeconds', 3, 'nWindows', 3, 'muscleSmoothSeconds', 2, 'blinkHalfSeconds', 0.5, ...
                  'topoHalfSeconds', 0.05, 'blinkElectrodeFallback', NaN);
    sheet = RELAX_screen_sheet_data(paths, loaded.RELAX_cfg, opts);
    assert(numel(sheet.muscle) == 3 && numel(sheet.blink) == 3 && iscell(sheet.interpolatedLabels));
    assert(numel(sheet.muscleRelativePower) == 3 && numel(sheet.blinkResidual) == 3);
    shown = [sheet.muscle([sheet.muscle.available])];
    for w = shown, assert(abs(w.t(end) - w.t(1) - 3) < 2 / 500); end
    starts = arrayfun(@(w) w.t(1), shown);
    assert(all(diff(sort(starts)) >= 3 - 1e-9) || numel(starts) < 2, 'muscle windows do not overlap');
    for w = [sheet.muscle sheet.blink]
        assert(w.available || strlength(w.note) > 0, 'every missing window explains itself');
        assert(~w.available || strlength(sheet.baselineNote) > 0 || ~isempty(w.before));
    end
    unavailable = ~[sheet.muscle.available];
    assert(all(isnan(sheet.muscleRelativePower(unavailable))));
    assert(all(isnan(sheet.blinkResidual(~[sheet.blink.available]))));
    opts.nWindows = 30;
    crowded = RELAX_screen_sheet_data(paths, loaded.RELAX_cfg, opts);
    missingBlinks = crowded.blink(~[crowded.blink.available]);
    assert(~isempty(missingBlinks) && all(startsWith([missingBlinks.note], "Only")), 'too few blinks: notes say how many');
    assert(all(contains([missingBlinks.note], 's apart')) || crowded.blinkN == nnz([crowded.blink.available]));
    missingMuscle = crowded.muscle(~[crowded.muscle.available]);
    assert(~isempty(missingMuscle) && all(contains([missingMuscle.note], 'muscle burst')) && all(contains([missingMuscle.note], 'at least 3 s apart')));
    opts.nWindows = 3;
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
% A short transient at a segment's end must not outrank a longer burst further in, and the
% window shown for the burst slides to stay inside its segment.
    rng(7);
    srate = 500; seconds = 20; nChan = 32;
    chanlocs = make_qc_chanlocs(nChan);
    eeg = eeg_emptyset(); eeg.srate = srate; eeg.nbchan = nChan; eeg.pnts = seconds * srate;
    eeg.trials = 1; eeg.xmin = 0; eeg.chanlocs = chanlocs;
    eeg.data = randn(nChan, eeg.pnts);
    secondsToSamples = @(t0, t1) round(t0 * srate) + 1 : round(t1 * srate);
    carrier = @(span, hz) sin(2 * pi * hz * (span - 1) / srate);
    burstSpan = secondsToSamples(10, 13.5);
    envelope = exp(-(((burstSpan - 1) / srate - 11.2) / 0.6).^2);
    eeg.data(:, burstSpan) = eeg.data(:, burstSpan) + 20 * envelope .* carrier(burstSpan, 50);
    transientSpan = secondsToSamples(9.8, 10);
    eeg.data(:, transientSpan) = eeg.data(:, transientSpan) + 35 * carrier(transientSpan, 50);
    eeg.event = struct('type', 'boundary', 'latency', 10 * srate + 0.5, 'duration', 0);
    eeg = eeg_checkset(eeg);
    eeg.allchan = chanlocs;
    eeg.RELAX.ListOfChannelsAfterRejections = {chanlocs.labels};

    folder = tempname; mkdir(folder); cleanup = onCleanup(@() rmdir(folder, 's'));
    pop_saveset(eeg, 'filename', 'synthetic_RELAX.set', 'filepath', folder);
    paths = struct('cleaned', fullfile(folder, 'synthetic_RELAX.set'), 'extremes', fullfile(folder, 'absent.set'));
    opts = struct('windowSeconds', 3, 'nWindows', 3, 'muscleSmoothSeconds', 2, 'blinkHalfSeconds', 0.5, ...
                  'topoHalfSeconds', 0.05, 'blinkElectrodeFallback', NaN);

    sheet = RELAX_screen_sheet_data(paths, loadedCfg, opts);
    strongest = sheet.muscle(1);
    assert(strongest.available && abs(mean(strongest.t) - 11.2) < 0.5, 'the strongest window is the mid-segment burst, not the transient at 9.8-10 s');
    assert(strongest.t(1) >= 10 - 1 / srate && isempty(strongest.boundaries), 'the burst window starts inside its segment and crosses no boundary');

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
