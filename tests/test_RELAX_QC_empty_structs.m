function test_RELAX_QC_empty_structs()
% Summaries and windows are built by filling templates, so differently-filled ones must
% still concatenate into struct arrays.
    a = RELAX_QC_empty_summary(); b = RELAX_QC_empty_summary();
    b.status = "ok"; b.pop = struct('label', "E1", 'centreSeconds', 1, 'shadeSeconds', [1 1.1]);
    b.metricsCleaned.blinkRatio = 1.2;
    both = [a; b];
    assert(numel(both) == 2 && both(2).status == "ok" && both(1).status == "notCleaned", 'summaries concatenate');
    w = [RELAX_QC_empty_window("blink"), RELAX_QC_empty_window("pop")];
    assert(~any([w.available]) && w(2).kind == "pop", 'windows concatenate');
    fprintf('test_RELAX_QC_empty_structs: PASS\n');
end
