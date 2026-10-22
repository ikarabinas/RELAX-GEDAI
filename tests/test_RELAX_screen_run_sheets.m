function test_RELAX_screen_run_sheets()
% End-to-end sheets stage on the one-file smoke fixture plus a ghost row with no input file.
    outDir = tempname; cleanup = onCleanup(@() rmdir(outDir, 's'));
    fixtures = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
    opts = struct('myPath', [fixtures filesep], 'outputFolderName', 'RELAXProcessed', 'outDir', outDir, ...
                  'cutoff', 3, 'sampleSize', 10, 'seed', 20261022, 'windowSeconds', 3, 'nWindows', 3, ...
                  'muscleSmoothSeconds', 2, 'blinkHalfSeconds', 0.5, 'topoHalfSeconds', 0.05);
    RELAX_screen_run_score(opts);

    loaded = load(fullfile(outDir, 'residual_screen.mat'));
    ghost = loaded.screen(loaded.screen.File == "boundary_fixture", :);
    ghost.File = "ghost_file"; ghost.ScreenOutcome = "fail";
    screen = [loaded.screen; ghost]; meta = loaded.meta; %#ok<NASGU>
    save(fullfile(outDir, 'residual_screen.mat'), 'screen', 'meta');

    screen = RELAX_screen_run_sheets(opts, {});
    assert(screen.SheetStatus(screen.File == "ghost_file") == "failed" && strlength(screen.SheetMessage(screen.File == "ghost_file")) > 0);
    assert(screen.SheetStatus(screen.File == "boundary_fixture") == "ok", 'one failure does not stop the others');
    assert(isfile(fullfile(outDir, 'review', 'boundary_fixture.png')));

    saved = load(fullfile(outDir, 'residual_screen.mat'));
    assert(isequal(saved.screen.SheetStatus, screen.SheetStatus), 'the table on disk carries the sheet status');
    csv = readtable(fullfile(outDir, 'residual_screen.csv'));
    assert(isnumeric(csv.InReviewSample) && all(ismember(csv.InReviewSample, [0 1])), 'the rewritten CSV keeps 0/1 columns');
    assert(ismember("ghost_file", string(csv.File)) && ismember("failed", string(csv.SheetStatus)));

    onlyOne = RELAX_screen_run_sheets(opts, "boundary_fixture");
    assert(onlyOne.SheetStatus(onlyOne.File == "boundary_fixture") == "ok" && onlyOne.SheetStatus(onlyOne.File == "ghost_file") == "failed", ...
           'a string array selects exactly those files and leaves other rows alone');

    try
        RELAX_screen_run_sheets(opts, {'nope'}); error('should have thrown');
    catch err
        assert(strcmp(err.identifier, 'RELAX_screen:unknownFile'), err.message);
    end
    fprintf('test_RELAX_screen_run_sheets: PASS\n');
end
