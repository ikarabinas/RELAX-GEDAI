function test_RELAX_QC_fill_window()
    EEG = eeg_emptyset(); EEG.srate = 100; EEG.nbchan = 2; EEG.pnts = 1000; EEG.trials = 1; EEG.xmin = 0;
    EEG.data = randn(2, 1000); EEG.chanlocs = struct('labels', {'E1', 'E2'}); EEG = eeg_checkset(EEG);
    w = RELAX_QC_fill_window(RELAX_QC_empty_window("muscle"), 5, EEG, EEG, {'E2', 'E1'}, 3);
    assert(numel(w.t) == 301 && abs(w.t(1) - 3.5) < 1e-12 && isequal(w.channels, {'E2', 'E1'}));
    w = RELAX_QC_fill_window(RELAX_QC_empty_window("muscle"), 0.2, EEG, EEG, {'E1'}, 3);
    assert(w.t(1) == 0, 'window clamps at the file start');
    w = RELAX_QC_fill_window(RELAX_QC_empty_window("muscle"), 5, [], EEG, {'E1', 'E9'}, 3);
    assert(isempty(w.before) && isempty(w.removed) && isequal(w.channels, {'E1'}) && w.available);
    fprintf('test_RELAX_QC_fill_window: PASS\n');
end
