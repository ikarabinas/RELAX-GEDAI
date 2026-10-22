function test_RELAX_QC_blink_average()
    srate = 500; pnts = 20000; chanlocs = make_qc_chanlocs(8);
    base = eeg_emptyset(); base.srate = srate; base.nbchan = 8; base.pnts = pnts; base.trials = 1; base.xmin = 0;
    base.chanlocs = chanlocs; base.data = zeros(8, pnts);
    peaks = [3000 9000 15000];
    for p = peaks, base.data(1:2, p-50:p+50) = 100; end
    base = eeg_checkset(base);
    clean = base; clean.data = 0.1 * base.data;
    avg = RELAX_QC_blink_average(base, clean, peaks, {'E1', 'E2'}, 250, 0.05);
    assert(abs(max(avg.erpBaseline) - 100) < 1e-9 && abs(max(avg.erpCleaned) - 10) < 1e-9);
    assert(abs(max(avg.erpRemoved) - 90) < 1e-9 && numel(avg.topoCleaned) == 8);
    avg = RELAX_QC_blink_average([], clean, peaks, {'E1', 'E2'}, 250, 0.05);
    assert(all(isnan(avg.erpBaseline)) && all(isnan(avg.topoRemoved)) && abs(max(avg.erpCleaned) - 10) < 1e-9);
    fprintf('test_RELAX_QC_blink_average: PASS\n');
end
