function test_RELAX_QC_fig_examples()
% One figure per example file; a window that could not be extracted shows its note; no
% examples at all still produces a figure saying so.

    outDir = tempname; cleanup = onCleanup(@() rmdir(outDir, 's'));
    cache = make_qc_cache_fixture(4);
    caption = RELAX_QC_fig_examples(cache, outDir);
    assertSaved(outDir, 'fig1_examples_1_f1');
    assertSaved(outDir, 'fig1_examples_2_f2');
    assert(startsWith(caption, "## Figure 1") && contains(caption, "reconstruction"), 'caption: %s', caption);

    empty = make_qc_cache_fixture(0);
    RELAX_QC_fig_examples(empty, fullfile(outDir, 'empty'));
    assertSaved(fullfile(outDir, 'empty'), 'fig1_examples_none');

    fprintf('test_RELAX_QC_fig_examples: PASS\n');
end

function assertSaved(outDir, name)
    for ext = [".png", ".pdf"]
        f = fullfile(outDir, name + ext);
        assert(isfile(f) && dir(f).bytes > 0, 'missing %s', f);
    end
end
