function test_RELAX_QC_cleaned_suffix()
% A cleaned file is saved as _RELAX_GEDAI.set when GEDAI ran, else _RELAX.set.

    processed = tempname; cleanedDir = fullfile(processed, 'Cleaned_Data'); mkdir(cleanedDir);
    cleanup = onCleanup(@() rmdir(processed, 's'));
    fclose(fopen(fullfile(cleanedDir, 'a_RELAX.set'), 'w'));
    fclose(fopen(fullfile(cleanedDir, 'b_RELAX_GEDAI.set'), 'w'));

    assert(RELAX_QC_cleaned_suffix(cleanedDir, "a") == "_RELAX.set", 'plain RELAX file');
    assert(RELAX_QC_cleaned_suffix(cleanedDir, "b") == "_RELAX_GEDAI.set", 'GEDAI file');
    assert(RELAX_QC_cleaned_suffix(cleanedDir, "missing") == "_RELAX.set", 'absent file defaults to plain');

    inputs = struct('name', 'b.set', 'folder', processed);
    paths = RELAX_QC_file_paths(processed, inputs, "b");
    assert(strcmp(paths.cleaned, fullfile(cleanedDir, 'b_RELAX_GEDAI.set')), 'file_paths follows the suffix: %s', paths.cleaned);
    fprintf('test_RELAX_QC_cleaned_suffix: PASS\n');
end
