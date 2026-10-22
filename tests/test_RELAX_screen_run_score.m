function test_RELAX_screen_run_score()
% End-to-end score stage on the one-file smoke fixture: every metric has zero spread.
    outDir = tempname; cleanup = onCleanup(@() rmdir(outDir, 's'));
    fixtures = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
    opts = struct('myPath', [fixtures filesep], 'outputFolderName', 'RELAXProcessed', 'outDir', outDir, ...
                  'cutoff', 3, 'sampleSize', 10, 'seed', 20261022);
    [screen, meta] = RELAX_screen_run_score(opts);
    assert(any(screen.File == "boundary_fixture"));
    for f = ["residual_screen.mat", "residual_screen.csv", "README.md"]
        assert(isfile(fullfile(outDir, f)), 'missing %s', f);
    end
    readme = fileread(fullfile(outDir, 'README.md'));
    assert(contains(readme, "cutoff 3") && contains(readme, "cannot flag any file in this run"));
    assert(~isempty(meta.zeroScaleMetrics));
    fprintf('test_RELAX_screen_run_score: PASS\n');
end
