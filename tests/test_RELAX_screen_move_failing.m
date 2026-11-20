function test_RELAX_screen_move_failing()
% Failing files move to Cleaned_Data/bad_QC_files; passing, unscreened and notCleaned files stay.
    processed = tempname; cleanup = onCleanup(@() rmdir(processed, 's'));
    cleanedDir = fullfile(processed, 'Cleaned_Data');
    mkdir(cleanedDir);
    for name = ["a_RELAX_GEDAI.set", "b_RELAX.set", "b_RELAX.fdt", "c_RELAX.set", "d_RELAX.set"]
        writelines("x", fullfile(cleanedDir, name));
    end
    outDir = fullfile(processed, 'Residual_screen');
    mkdir(outDir);
    screen = table(["a"; "b"; "c"; "d"; "e"], ["fail"; "fail"; "pass"; "unscreened"; "notCleaned"], ...
                   'VariableNames', {'File', 'ScreenOutcome'});
    meta = struct();
    save(fullfile(outDir, 'residual_screen.mat'), 'screen', 'meta');
    opts = struct('myPath', [fileparts(processed) filesep], 'outputFolderName', char(lastPart(processed)), 'outDir', outDir);

    moved = RELAX_screen_move_failing(opts);

    badDir = fullfile(cleanedDir, 'bad_QC_files');
    assert(isequal(sort(moved), ["a"; "b"]));
    for name = ["a_RELAX_GEDAI.set", "b_RELAX.set", "b_RELAX.fdt"]
        assert(isfile(fullfile(badDir, name)) && ~isfile(fullfile(cleanedDir, name)), '%s should have moved', name);
    end
    for name = ["c_RELAX.set", "d_RELAX.set"]
        assert(isfile(fullfile(cleanedDir, name)), '%s should have stayed', name);
    end

    assert(isempty(RELAX_screen_move_failing(opts)), 'a second move finds nothing left to move');

    threw = false;
    try
        RELAX_screen_assert_no_moved_files(opts);
    catch err
        threw = strcmp(err.identifier, 'RELAX_screen:movedFilesPresent');
    end
    assert(threw, 'scoring with files in bad_QC_files must refuse');
    fprintf('test_RELAX_screen_move_failing: PASS\n');
end

function name = lastPart(p)
    [~, name] = fileparts(p);
    name = string(name);
end
