function test_RELAX_screen_draw_sheet()
% The reason text names flagging or missing metrics (missing slope never listed, extra note
% appended, notCleaned handled); a sheet is saved for a failing file, a file without a
% baseline, a file without blinks, a screen with a zero-scale metric, and a file whose z is
% far beyond the strip plots' range.

    reviewDir = tempname; cleanup = onCleanup(@() removeIfPresent(reviewDir));
    figuresBefore = numel(findall(groot, 'Type', 'figure'));
    [screen, meta, sheet] = syntheticSheet();

    assert(RELAX_screen_reason(screen(2, :), meta) == "MuscleProportion is 4.2 robust SD above the dataset median (cutoff 3.0)");
    assert(RELAX_screen_reason(screen(3, :), meta) == "Not screened: BlinkRatio missing");
    assert(RELAX_screen_reason(screen(1, :), meta, "Baseline not aligned") == "All metrics within cutoff 3.0 | Baseline not aligned");
    assert(RELAX_screen_reason(screen(1, :), meta, "") == "All metrics within cutoff 3.0");

    twoMissing = screen(3, :); twoMissing.NInterpolated = NaN;
    assert(RELAX_screen_reason(twoMissing, meta) == "Not screened: BlinkRatio, NInterpolated missing");
    slopeAndRatioMissing = screen(3, :); slopeAndRatioMissing.MuscleStrength = NaN;
    assert(RELAX_screen_reason(slopeAndRatioMissing, meta) == "Not screened: BlinkRatio missing", 'a missing slope is never listed');
    twoFlags = screen(2, :); twoFlags.BlinkRatio_flag = true; twoFlags.BlinkRatio_z = 5.04;
    assert(RELAX_screen_reason(twoFlags, meta) == ...
        "MuscleProportion is 4.2 robust SD above the dataset median (cutoff 3.0); BlinkRatio is 5.0 robust SD above the dataset median (cutoff 3.0)");
    notCleaned = screen(1, :); notCleaned.ScreenOutcome = "notCleaned";
    assert(RELAX_screen_reason(notCleaned, meta) == "Not cleaned by the pipeline");

    RELAX_screen_draw_sheet(sheet, screen(2, :), screen, meta, reviewDir);

    noBaseline = sheet; noBaseline.baselineNote = "Baseline not aligned";
    for k = 1:3
        [noBaseline.muscle(k).before, noBaseline.muscle(k).removed] = deal([]);
        [noBaseline.blink(k).before, noBaseline.blink(k).removed] = deal([]);
    end
    noBaseline.blinkAverage = withoutBaseline(sheet.blinkAverage);
    RELAX_screen_draw_sheet(noBaseline, screen(1, :), screen, meta, reviewDir);

    noBlinks = sheet; noBlinks.blinkN = 0; noBlinks.blinkAverage = [];
    noBlinks.blink = repmat(RELAX_QC_empty_window("blink"), 1, 3); [noBlinks.blink.note] = deal("No blinks in this file");
    noBlinks.blinkResidual = nan(1, 3);
    RELAX_screen_draw_sheet(noBlinks, screen(3, :), screen, meta, reviewDir);

    offScale = screen; offScale.MuscleProportion_z(2:3) = [40; 55];
    RELAX_screen_draw_sheet(sheet, offScale(2, :), offScale, meta, fullfile(reviewDir, 'offScale'));
    assert(isfile(fullfile(reviewDir, 'offScale', 'f2.png')));

    allZero = screen; allZero.MuscleStrength_z(:) = NaN;
    meta.zeroScaleMetrics = "MuscleStrength";
    RELAX_screen_draw_sheet(sheet, allZero(1, :), allZero, meta, fullfile(reviewDir, 'zeroScale'));

    for f = ["f2", "f1", "f3"], for ext = [".png", ".pdf"]
        assert(isfile(fullfile(reviewDir, f + ext)), 'missing %s%s', f, ext);
    end, end
    assert(isfile(fullfile(reviewDir, 'zeroScale', 'f1.png')));
    assert(numel(findall(groot, 'Type', 'figure')) == figuresBefore, 'the drawer leaves no figure open');
    fprintf('test_RELAX_screen_draw_sheet: PASS\n');
end

function [screen, meta, sheet] = syntheticSheet()
    nChan = 32; srate = 100; seconds = 60;
    chanlocs = make_qc_chanlocs(nChan);
    cleaned = eeg_emptyset(); cleaned.srate = srate; cleaned.nbchan = nChan; cleaned.pnts = seconds * srate;
    cleaned.trials = 1; cleaned.xmin = 0; cleaned.chanlocs = chanlocs;
    cleaned.data = randn(nChan, cleaned.pnts);
    cleaned = eeg_checkset(cleaned);
    before = cleaned; before.data = cleaned.data + 5 * randn(nChan, cleaned.pnts);

    channelsShown = {'E1', 'E2', 'E5', 'E9', 'E12', 'E20'};
    window = @(kind, centre) RELAX_QC_fill_window(RELAX_QC_empty_window(kind), centre, before, cleaned, channelsShown, 3);
    sheet.chanlocs = chanlocs;
    sheet.interpolatedLabels = {'E5'};
    sheet.baselineNote = "";
    sheet.muscle = [window("muscle", 10), window("muscle", 25), window("muscle", 40)];
    sheet.blink = [window("blink", 12), window("blink", 30), window("blink", 50)];
    sheet.muscleRelativePower = [4.2 3.1 2.5];
    sheet.blinkResidual = [31 22 18];
    peaks = [1500 3000 4500];
    sheet.blinkN = numel(peaks);
    sheet.blinkAverage = RELAX_QC_blink_average(before, cleaned, peaks, {'E1', 'E2'}, 50, 0.05);

    T = table(["f1"; "f2"; "f3"], ["p1"; "p1"; "p2"], ["cleaned"; "cleaned"; "cleaned"], 'VariableNames', {'File', 'Participant', 'Outcome'});
    metrics = struct('MuscleProportion', {0.10; 0.30; 0.12}, 'MuscleStrength', {1; 2; 1.5}, ...
                     'BlinkRatio', {0.2; 0.3; NaN}, 'NInterpolated', {2; 3; 1});
    [screen, meta] = RELAX_screen_score_table(T, metrics, struct('cutoff', 3, 'sampleSize', 10, 'seed', 1));
    screen.MuscleProportion_z(2) = 4.2;
    screen.MuscleProportion_flag = [false; true; false];
    screen.ScreenOutcome = ["pass"; "fail"; "unscreened"];
    screen.FlaggedBy = ["", "MuscleProportion", ""]';
end

function avg = withoutBaseline(avg)
    avg.erpBaseline(:) = NaN; avg.erpRemoved(:) = NaN;
    avg.topoBaseline(:) = NaN; avg.topoRemoved(:) = NaN;
end

function removeIfPresent(folder)
    if isfolder(folder), rmdir(folder, 's'); end
end
