function test_RELAX_QC_extract_robustness()
% One bad example must not lose a long run; checkpoint caches are flagged; the captions
% footer reads as one sentence; excised intervals are clipped to the example window.

    windows = RELAX_QC_try_example_windows(struct('file', "f1", 'centres', struct('blink', NaN, 'muscle', NaN, 'pop', NaN)), ...
        RELAX_QC_empty_summary(), struct('input', 'nope.set', 'extremes', 'nope_E.set', 'cleaned', 'nope_C.set'), struct(), struct());
    assert(numel(windows) == 4 && isequal([windows.kind], ["blink", "muscle", "pop", "popKept"]), 'four windows, in figure order');
    assert(~any([windows.available]), 'all unavailable after an error');
    assert(all(strlength([windows.note]) > 0) && all(contains([windows.note], "nope")), 'each note names the error: %s', windows(1).note);

    cache = make_qc_cache_fixture(2);
    cache.meta.partial = false;
    lastwarn('');
    RELAX_QC_warn_if_partial(cache, 'x.mat');
    assert(isempty(lastwarn()), 'complete cache does not warn');
    cache.meta.partial = true;
    RELAX_QC_warn_if_partial(cache, 'x.mat');
    [~, warnId] = lastwarn();
    assert(strcmp(warnId, 'RELAX_QC:partialCache'), 'checkpoint cache warns');

    cache = make_qc_cache_fixture(2);
    assert(contains(RELAX_QC_captions_footer(cache, 'c.mat'), "1 file skipped for alignment: f3 (pnts differ: baseline 10, cleaned 9)._"), 'one skipped');
    cache.summaries(1).status = "skipped"; cache.summaries(1).skipReason = "why";
    footer = RELAX_QC_captions_footer(cache, 'c.mat');
    assert(contains(footer, "2 files skipped for alignment: f1 (why); f3 (") && endsWith(strtrim(footer), "_") && startsWith(footer, "_"), 'many skipped: %s', footer);
    cache.summaries(:) = RELAX_QC_empty_summary(); [cache.summaries.status] = deal("ok");
    footer = RELAX_QC_captions_footer(cache, 'c.mat');
    assert(contains(footer, "No files were skipped for alignment._") && ~contains(footer, "0 file"), 'none skipped: %s', footer);

    periods = [1 500; 1001 1250; 5001 6000];
    inside = RELAX_QC_excised_intervals(periods, 500, [0.5 3]);
    assert(isequal(size(inside), [2 2]) && max(abs(inside(:) - [0.5; 2.0; 0.998; 2.498])) < 1e-9, 'clipped to window, overlapping only: %s', mat2str(inside));
    assert(isequal(size(RELAX_QC_excised_intervals(zeros(0, 2), 500, [0 10])), [0 2]), 'no periods');

    fprintf('test_RELAX_QC_extract_robustness: PASS\n');
end
