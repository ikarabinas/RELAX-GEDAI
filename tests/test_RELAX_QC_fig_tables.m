function test_RELAX_QC_fig_tables()
% Figures 3 and 4 save PNG+PDF and a caption, including when no file aligned; Figure 4
% shows cleaned and excluded/errored files but never files not yet processed.

    outDir = tempname; cleanup = onCleanup(@() rmdir(outDir, 's'));
    cache = make_qc_cache_fixture(4);

    caption = RELAX_QC_fig_paired_metrics(cache, outDir);
    assertSaved(outDir, 'fig3_paired_metrics');
    assert(startsWith(caption, "## Figure 3") && contains(caption, "4 files") ...
           && contains(caption, "floor") && contains(caption, "-0.31") && ~contains(caption, "exceed it"), 'caption: %s', caption);

    caption = RELAX_QC_fig_data_taken(cache, outDir);
    assertSaved(outDir, 'fig4_data_taken');
    assert(contains(caption, "6 files shown"), '4 ok + 1 skipped + 1 excluded, not the notProcessed one: %s', caption);
    assert(contains(caption, "1 excluded"), 'counts exclusions: %s', caption);
    assert(contains(caption, "rejection budget") && contains(caption, "25%"), 'states budget: %s', caption);

    empty = make_qc_cache_fixture(0);
    RELAX_QC_fig_paired_metrics(empty, fullfile(outDir, 'empty'));
    assertSaved(fullfile(outDir, 'empty'), 'fig3_paired_metrics');
    RELAX_QC_fig_data_taken(empty, fullfile(outDir, 'empty'));
    assertSaved(fullfile(outDir, 'empty'), 'fig4_data_taken');

    nothingShown = make_qc_cache_fixture(0);
    nothingShown.tables.Outcome(:) = "notProcessed";
    noneDir = fullfile(outDir, 'none');
    assert(startsWith(RELAX_QC_fig_paired_metrics(nothingShown, noneDir), "## Figure 3"));
    caption = RELAX_QC_fig_data_taken(nothingShown, noneDir);
    assertSaved(noneDir, 'fig3_paired_metrics');
    assertSaved(noneDir, 'fig4_data_taken');
    assert(startsWith(caption, "## Figure 4") && contains(caption, "0 files shown") && contains(caption, "25%"), 'no files: %s', caption);

    fprintf('test_RELAX_QC_fig_tables: PASS\n');
end

function assertSaved(outDir, name)
    for ext = [".png", ".pdf"]
        f = fullfile(outDir, name + ext);
        assert(isfile(f) && dir(f).bytes > 0, 'missing %s', f);
    end
end
