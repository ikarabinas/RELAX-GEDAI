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
    fprintf('test_RELAX_screen_sheet_data: PASS\n');
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
