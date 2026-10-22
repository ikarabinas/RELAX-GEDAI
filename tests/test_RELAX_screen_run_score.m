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
    assert(~isempty(meta.zeroScaleMetrics));
    csv = readtable(fullfile(outDir, 'residual_screen.csv'));
    assert(isnumeric(csv.InReviewSample) && all(ismember(csv.InReviewSample, [0 1])));
    assert(isnumeric(csv.MuscleProportion_flag) && all(ismember(csv.MuscleProportion_flag, [0 1])));
    test_unreadable_header(fixtures, screen);
    fprintf('test_RELAX_screen_run_score: PASS\n');
end

function test_unreadable_header(fixtures, healthyScreen)
% A copy of the fixture's stats beside an unreadable cleaned file: the table's muscle values
% survive, the file is unscreened, and the load error is reported.
    root = tempname; cleanup = onCleanup(@() rmdir(root, 's'));
    processed = fullfile(root, 'RELAXProcessed');
    mkdir(fullfile(processed, 'Cleaned_Data'));
    copyfile(fullfile(fixtures, 'RELAXProcessed', 'PerParticipantStats'), fullfile(processed, 'PerParticipantStats'));
    writelines("not a set file", fullfile(processed, 'Cleaned_Data', 'boundary_fixture_RELAX.set'));
    warningState = warning('off', 'RELAX_screen:headerUnreadable'); restoreWarning = onCleanup(@() warning(warningState));
    opts = struct('myPath', [root filesep], 'outputFolderName', 'RELAXProcessed', 'outDir', fullfile(root, 'out'), ...
                  'cutoff', 3, 'sampleSize', 10, 'seed', 1);
    screen = RELAX_screen_run_score(opts);
    assert(screen.ScreenOutcome == "unscreened" && strlength(screen.SheetMessage) > 0);
    assert(isnan(screen.NInterpolated) && isnan(screen.BlinkRatio));
    assert(screen.MuscleProportion == healthyScreen.MuscleProportion && screen.MuscleStrength == healthyScreen.MuscleStrength);
end
