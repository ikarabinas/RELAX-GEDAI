function test_RELAX_QC_fig_spectra_topo()
% Figures 2 and 5 save and caption correctly, leave files with too few blinks out of the
% group topography, and survive a dataset with no aligned files.

    outDir = tempname; cleanup = onCleanup(@() rmdir(outDir, 's'));
    cache = make_qc_cache_fixture(4);

    caption = RELAX_QC_fig_spectra(cache, outDir);
    assertSaved(outDir, 'fig2_spectra');
    assert(startsWith(caption, "## Figure 2") && contains(caption, "4 files"), 'caption: %s', caption);

    caption = RELAX_QC_fig_blink_topo(cache, outDir);
    assertSaved(outDir, 'fig5_blink_topography');
    assert(contains(caption, "baseline-corrected") && contains(caption, "−500 to −300 ms"), 'states baseline: %s', caption);
    assert(contains(caption, "3 files") && contains(caption, "1 file with fewer than 5 blinks"), ...
           'f4 has 2 blinks and must be left out of the group row: %s', caption);

    noBlinkExample = cache;
    noBlinkExample.summaries(2).blinkN = 0;          % the "worst blink" pick has no blinks: fall back
    RELAX_QC_fig_blink_topo(noBlinkExample, fullfile(outDir, 'fallback'));
    assertSaved(fullfile(outDir, 'fallback'), 'fig5_blink_topography');

    missingInput = cache;
    for i = 1:4, missingInput.summaries(i).psdInput = nan(size(missingInput.summaries(i).freqs)); end
    RELAX_QC_fig_spectra(missingInput, fullfile(outDir, 'noinput'));
    assertSaved(fullfile(outDir, 'noinput'), 'fig2_spectra');

    empty = make_qc_cache_fixture(0);
    RELAX_QC_fig_spectra(empty, fullfile(outDir, 'empty'));
    assertSaved(fullfile(outDir, 'empty'), 'fig2_spectra');
    RELAX_QC_fig_blink_topo(empty, fullfile(outDir, 'empty'));
    assertSaved(fullfile(outDir, 'empty'), 'fig5_blink_topography');

    unprocessed = cache;
    unprocessed.tables.Outcome(:) = "notProcessed";
    for i = 1:numel(unprocessed.summaries), unprocessed.summaries(i) = RELAX_QC_empty_summary(); end
    caption = RELAX_QC_fig_spectra(unprocessed, fullfile(outDir, 'unprocessed'));
    assertSaved(fullfile(outDir, 'unprocessed'), 'fig2_spectra');
    assert(contains(caption, "0 files"), 'caption: %s', caption);
    caption = RELAX_QC_fig_blink_topo(unprocessed, fullfile(outDir, 'unprocessed'));
    assertSaved(fullfile(outDir, 'unprocessed'), 'fig5_blink_topography');
    assert(contains(caption, "0 files") && contains(caption, "no file with blinks"), 'caption: %s', caption);

    noneEligible = cache;
    for i = 1:4, noneEligible.summaries(i).blinkN = 1; end
    caption = RELAX_QC_fig_blink_topo(noneEligible, fullfile(outDir, 'noneligible'));
    assertSaved(fullfile(outDir, 'noneligible'), 'fig5_blink_topography');
    assert(contains(caption, "0 files") && contains(caption, "4 files with fewer than 5 blinks"), 'caption: %s', caption);

    fprintf('test_RELAX_QC_fig_spectra_topo: PASS\n');
end

function assertSaved(outDir, name)
    for ext = [".png", ".pdf"]
        f = fullfile(outDir, name + ext);
        assert(isfile(f) && dir(f).bytes > 0, 'missing %s', f);
    end
end
