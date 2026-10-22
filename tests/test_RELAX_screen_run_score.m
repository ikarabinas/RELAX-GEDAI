function test_RELAX_screen_run_score()
% End-to-end score stage on the one-file smoke fixture: every metric has zero spread.
    outDir = tempname; cleanup = onCleanup(@() rmdir(outDir, 's'));
    fixtures = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
    opts = struct('myPath', [fixtures filesep], 'outputFolderName', 'RELAXProcessed', 'outDir', outDir, ...
                  'cutoff', 3, 'sampleSize', 10, 'seed', 20261022);
    [screen, meta] = RELAX_screen_run_score(opts);
    assert(any(screen.File == "boundary_fixture"));
    assert(screen.NInterpolated(screen.File == "boundary_fixture") == 5, 'interpolated = allchan minus kept channels');
    for f = ["residual_screen.mat", "residual_screen.csv", "README.md"]
        assert(isfile(fullfile(outDir, f)), 'missing %s', f);
    end
    readme = fileread(fullfile(outDir, 'README.md'));
    assert(contains(readme, "cutoff 3") && contains(readme, "cannot flag any file in this run"));
    assert(contains(readme, "1.4826 * MAD"), 'README explains the scale columns');
    assert(contains(readme, "Open rings mark electrodes that were interpolated") && contains(readme, "3 s each"), 'README guide describes the drawn sheet');
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
                  'cutoff', 3, 'sampleSize', 10, 'seed', 1);
    screen = RELAX_screen_run_score(opts);
    assert(screen.ScreenOutcome == "unscreened" && strlength(screen.SheetMessage) > 0);
    assert(contains(screen.SheetMessage, expectedMessagePart));
    assert(isnan(screen.NInterpolated) && isnan(screen.BlinkRatio));
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
