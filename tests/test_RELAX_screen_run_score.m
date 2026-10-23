function test_RELAX_screen_run_score()
% End-to-end score stage on the one-file smoke fixture: every metric has zero spread.
    outDir = tempname; cleanup = onCleanup(@() rmdir(outDir, 's'));
    fixtures = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
    opts = struct('myPath', [fixtures filesep], 'outputFolderName', 'RELAXProcessed', 'outDir', outDir, ...
                  'cutoff', 3, 'sampleSize', 10, 'seed', 20261022, 'holeDepthLimit', 1.4);
    [screen, meta] = RELAX_screen_run_score(opts);
    assert(any(screen.File == "boundary_fixture"));
    assert(screen.NInterpolated(screen.File == "boundary_fixture") == 5, 'interpolated = allchan minus kept channels');
    fixtureRow = screen.File == "boundary_fixture";
    reasonNames = ["impedance" "pop" "bridged" "PREP" "extreme" "muscle" "unknown"];
    fixtureReasonCounts = screen{fixtureRow, "NInterpolated_" + reasonNames};
    assert(screen.NInterpolated(fixtureRow) == sum(fixtureReasonCounts), 'the reason counts add up to the interpolated channels');
    assert(isequal(fixtureReasonCounts, [0 0 0 0 1 4 0]), 'the fixture header lists one extreme and four muscle channels');
    assert(screen.InterpolatedReasons(fixtureRow) == "extreme: E226; muscle: E81 E90 E205 E234");
    assert(isfinite(screen.HoleDepth(screen.File == "boundary_fixture")), 'the hole is measured from the electrode coordinates');
    for f = ["residual_screen.mat", "residual_screen.csv", "README.md"]
        assert(isfile(fullfile(outDir, f)), 'missing %s', f);
    end
    readme = fileread(fullfile(outDir, 'README.md'));
    assert(contains(readme, "cutoff 3") && contains(readme, "cannot flag any file in this run"));
    assert(contains(readme, "1.4826 * MAD"), 'README explains the scale columns');
    assert(contains(readme, "No before-wICA trace:") && contains(readme, "separate line for each metric") && contains(readme, "Clear that folder"), 'README guide explains the header lines and stale sheets');
    assert(contains(readme, "Open rings mark electrodes that were interpolated") && contains(readme, "3 s each"), 'README guide describes the drawn sheet');
    assert(contains(readme, "Amber filled electrodes mark the hole") && contains(readme, "NInterpolated_limit"), 'README guide describes the hole rule and the channel limit');
    assert(contains(readme, "InterpolatedReasons") && contains(readme, "NInterpolated_<reason>"), 'README explains the removal reasons');
    assert(contains(readme, "highest muscle slope") && contains(readme, "its channel and the threshold") && ~contains(readme, "20-80 Hz power as a multiple"), 'README guide describes the slope-ranked muscle windows');
    assert(contains(readme, "'(interp., <reason>)'") && contains(readme, "Interpolated for <reason>"), 'README guide describes the reasons drawn on the sheet');
    assert(contains(readme, "BlinkRatio or HoleDepth is missing"), 'README names the unscreening metrics');
    assert(~isnan(screen.BlinkRatio) || strlength(screen.SheetMessage) > 0, 'a missing BlinkRatio is explained');
    assert(~isempty(meta.zeroScaleMetrics));
    csv = readtable(fullfile(outDir, 'residual_screen.csv'));
    assert(isnumeric(csv.InReviewSample) && all(ismember(csv.InReviewSample, [0 1])));
    assert(isnumeric(csv.MuscleProportion_flag) && all(ismember(csv.MuscleProportion_flag, [0 1])));
    assertUnscreenedWithHeader(fixtures, screen, @writeTextFile, "");
    assertUnscreenedWithHeader(fixtures, screen, @writeHeaderWithoutAnyChannelInfo, "EEG.allchan");
    assertUnscreenedWithHeader(fixtures, screen, @writeHeaderWithoutRelaxField, "ListOfChannelsAfterRejections");
    fprintf('test_RELAX_screen_run_score: PASS\n');
end

function assertUnscreenedWithHeader(fixtures, healthyScreen, writeCleanedFile, expectedMessagePart)
% A copy of the fixture's stats beside a cleaned file whose header is unusable: the table's
% muscle values survive, the file is unscreened, and the reason is reported.
    root = tempname; cleanup = onCleanup(@() rmdir(root, 's'));
    processed = fullfile(root, 'RELAXProcessed');
    cleanedDir = fullfile(processed, 'Cleaned_Data');
    mkdir(cleanedDir);
    copyfile(fullfile(fixtures, 'RELAXProcessed', 'PerParticipantStats'), fullfile(processed, 'PerParticipantStats'));
    writeCleanedFile(cleanedDir, 'boundary_fixture_RELAX.set');
    warningState = warning('off', 'RELAX_screen:headerUnreadable'); restoreWarning = onCleanup(@() warning(warningState));
    opts = struct('myPath', [root filesep], 'outputFolderName', 'RELAXProcessed', 'outDir', fullfile(root, 'out'), ...
                  'cutoff', 3, 'sampleSize', 10, 'seed', 1, 'holeDepthLimit', 1.4);
    screen = RELAX_screen_run_score(opts);
    assert(screen.ScreenOutcome == "unscreened" && strlength(screen.SheetMessage) > 0);
    assert(contains(screen.SheetMessage, expectedMessagePart));
    assert(isnan(screen.NInterpolated) && isnan(screen.BlinkRatio) && isnan(screen.HoleDepth));
    assert(all(isnan(screen{:, "NInterpolated_" + ["impedance" "pop" "bridged" "PREP" "extreme" "muscle" "unknown"]})) && screen.InterpolatedReasons == "");
    assert(screen.MuscleProportion == healthyScreen.MuscleProportion && screen.MuscleStrength == healthyScreen.MuscleStrength);
end

function writeTextFile(cleanedDir, name)
    writelines("not a set file", fullfile(cleanedDir, name));
end

function writeHeaderWithoutAnyChannelInfo(cleanedDir, name)
    pop_saveset(eeg_emptyset, 'filename', name, 'filepath', cleanedDir, 'savemode', 'onefile');
end

function writeHeaderWithoutRelaxField(cleanedDir, name)
    EEG = eeg_emptyset;
    EEG.allchan = struct('labels', {'E1', 'E2'});
    pop_saveset(EEG, 'filename', name, 'filepath', cleanedDir, 'savemode', 'onefile');
end
