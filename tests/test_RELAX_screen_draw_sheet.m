function test_RELAX_screen_draw_sheet()
% The reason text names flagging or missing metrics (missing slope never listed, extra note
% appended, notCleaned handled); a sheet is saved for a failing file, a file without a
% baseline, a file without blinks, a screen with a zero-scale metric, a file whose z is
% far beyond the strip plots' range, a file whose hole depth is missing and a depth beyond its strip.

    reviewDir = tempname; cleanup = onCleanup(@() removeIfPresent(reviewDir));
    figuresBefore = numel(findall(groot, 'Type', 'figure'));
    [screen, meta, sheet] = syntheticSheet();

    muscleSentence = "MuscleProportion is 4.2 robust SD above the dataset median (cutoff 3.0)";
    holeSentence = "Interpolation hole 1.70× electrode spacing deep at E5 (limit 1.4)";
    assert(RELAX_screen_reason(screen(2, :), meta) == muscleSentence + "; " + holeSentence);
    holeOnly = screen(2, :); holeOnly.MuscleProportion_flag = false;
    assert(RELAX_screen_reason(holeOnly, meta) == holeSentence);
    channelLimit = screen(1, :); channelLimit.ScreenOutcome = "fail"; channelLimit.NInterpolated = 50; channelLimit.NInterpolated_flag = true;
    channelLimit.NInterpolated_limit = 50; channelLimit.NChannels = 203; channelLimit.MaxProportionDeleted = 0.25;
    assert(RELAX_screen_reason(channelLimit, meta) == "50 channels interpolated, the pipeline's limit of 50 (25% of 203): bad channels may remain");
    assert(RELAX_screen_reason(screen(3, :), meta) == "Not screened: BlinkRatio missing");
    passLine = "All metrics within cutoff 3.0, hole depth within limit 1.4";
    assert(RELAX_screen_reason(screen(1, :), meta, "Baseline not aligned") == passLine + " | No before-wICA trace: Baseline not aligned");
    assert(RELAX_screen_reason(screen(1, :), meta, "") == passLine);

    twoMissing = screen(3, :); twoMissing.HoleDepth = NaN;
    assert(RELAX_screen_reason(twoMissing, meta) == "Not screened: BlinkRatio, HoleDepth missing");
    slopeAndRatioMissing = screen(3, :); slopeAndRatioMissing.MuscleStrength = NaN;
    assert(RELAX_screen_reason(slopeAndRatioMissing, meta) == "Not screened: BlinkRatio missing", 'a missing slope is never listed');
    twoFlags = screen(2, :); twoFlags.BlinkRatio_flag = true; twoFlags.BlinkRatio_z = 5.04;
    assert(RELAX_screen_reason(twoFlags, meta) == muscleSentence + "; BlinkRatio is 5.0 robust SD above the dataset median (cutoff 3.0); " + holeSentence);
    notCleaned = screen(1, :); notCleaned.ScreenOutcome = "notCleaned";
    assert(RELAX_screen_reason(notCleaned, meta) == "Not cleaned by the pipeline");

    drawWithHoleMarkerWarningFatal(sheet, screen(2, :), screen, meta, reviewDir);

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

    missingDepth = twoMissing;
    RELAX_screen_draw_sheet(noBlinks, missingDepth, screen, meta, fullfile(reviewDir, 'missingDepth'));
    assert(isfile(fullfile(reviewDir, 'missingDepth', 'f3.png')));

    deepHoles = screen; deepHoles.HoleDepth(1:2) = [9; 1.7];
    RELAX_screen_draw_sheet(sheet, deepHoles(2, :), deepHoles, meta, fullfile(reviewDir, 'deepHoles'));
    assert(isfile(fullfile(reviewDir, 'deepHoles', 'f2.png')));

    manyLabels = arrayfun(@(k) sprintf('E%d', k), 1:24, 'UniformOutput', false);
    crowded = sheet; crowded.interpolatedLabels = manyLabels;
    crowded.interpolatedReasons = struct('labels', string(manyLabels), 'reasons', repmat("muscle", 1, 24), ...
                                         'groups', struct('reason', "muscle", 'labels', string(manyLabels)));
    crowdedRow = screen(2, :); crowdedRow.HoleElectrodes = strjoin(string(manyLabels), " ");
    drawWithHoleMarkerWarningFatal(crowded, crowdedRow, screen, meta, fullfile(reviewDir, 'crowded'));
    crowdedNoBlinks = noBlinks; crowdedNoBlinks.interpolatedLabels = manyLabels;
    crowdedNoBlinks.interpolatedReasons = crowded.interpolatedReasons;
    drawWithHoleMarkerWarningFatal(crowdedNoBlinks, crowdedRow, screen, meta, fullfile(reviewDir, 'crowdedNoBlinks'));
    assert(isfile(fullfile(reviewDir, 'crowded', 'f2.png')) && isfile(fullfile(reviewDir, 'crowdedNoBlinks', 'f2.png')));

    allZero = screen; allZero.BlinkRatio_z(:) = NaN;
    meta.zeroScaleMetrics = "BlinkRatio";
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
    sheet.interpolatedLabels = {'E5', 'E9'};
    sheet.interpolatedReasons = struct('labels', ["E5" "E9"], 'reasons', ["muscle" "impedance"], ...
                                       'groups', struct('reason', {"impedance", "muscle"}, 'labels', {"E9", "E5"}));
    sheet.baselineNote = "";
    sheet.muscle = [window("muscle", 10), window("muscle", 25), window("muscle", 40)];
    sheet.blink = [window("blink", 12), window("blink", 30), window("blink", 50)];
    sheet.muscleSlope = [-0.1 -0.2 -0.4];
    sheet.muscleSlopeChannel = ["E1" "E5" "E12"];
    sheet.muscleSlopeThreshold = -0.31;
    sheet.muscleEpochStart = [9.5 24.5 39.5];
    sheet.blinkResidual = [31 22 18];
    peaks = [1500 3000 4500];
    sheet.blinkN = numel(peaks);
    sheet.blinkAverage = RELAX_QC_blink_average(before, cleaned, peaks, {'E1', 'E2'}, 50, 0.05);

    T = table(["f1"; "f2"; "f3"], ["p1"; "p1"; "p2"], ["cleaned"; "cleaned"; "cleaned"], 'VariableNames', {'File', 'Participant', 'Outcome'});
    metrics = struct('MuscleProportion', {0.10; 0.30; 0.12}, 'MuscleStrength', {1; 2; 1.5}, ...
                     'BlinkRatio', {0.2; 0.3; NaN}, 'NInterpolated', {2; 3; 1}, 'HoleDepth', {1.0; 1.7; 1.0}, ...
                     'HoleDeepestElectrode', {""; "E5"; ""}, 'HoleElectrodes', {""; "E5"; ""}, 'Spacing', {0.02; 0.02; 0.02}, ...
                     'NInterpolatedLimit', {50; 50; 50}, 'NChannels', {203; 203; 203}, 'MaxProportionDeleted', {0.25; 0.25; 0.25}, ...
                     'InterpolatedReasons', {""; "unknown: E5"; ""}, 'NInterpolated_impedance', {0; 0; 0}, 'NInterpolated_pop', {0; 0; 0}, ...
                     'NInterpolated_bridged', {0; 0; 0}, 'NInterpolated_PREP', {0; 0; 0}, 'NInterpolated_extreme', {0; 0; 0}, ...
                     'NInterpolated_muscle', {0; 0; 0}, 'NInterpolated_unknown', {0; 1; 0});
    [screen, meta] = RELAX_screen_score_table(T, metrics, struct('cutoff', 3, 'sampleSize', 10, 'seed', 1, 'holeDepthLimit', 1.4));
    screen.MuscleProportion_z(2) = 4.2;
    screen.MuscleProportion_flag = [false; true; false];
    screen.ScreenOutcome = ["pass"; "fail"; "unscreened"];
    screen.FlaggedBy = ["", "MuscleProportion; HoleDepth", ""]';
end

function avg = withoutBaseline(avg)
    avg.erpBaseline(:) = NaN; avg.erpRemoved(:) = NaN;
    avg.topoBaseline(:) = NaN; avg.topoRemoved(:) = NaN;
end

function drawWithHoleMarkerWarningFatal(sheet, row, screen, meta, reviewDir)
    previous = warning('error', 'RELAX_screen:holeMarkers');
    restore = onCleanup(@() warning(previous));
    RELAX_screen_draw_sheet(sheet, row, screen, meta, reviewDir);
end

function removeIfPresent(folder)
    if isfolder(folder), rmdir(folder, 's'); end
end
