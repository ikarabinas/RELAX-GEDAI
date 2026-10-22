function test_RELAX_screen_blink_residual()
    EEG = eeg_emptyset(); EEG.srate = 500; EEG.nbchan = 3; EEG.pnts = 10000; EEG.trials = 1; EEG.xmin = 0;
    EEG.chanlocs = struct('labels', {'E1', 'E2', 'E3'}); EEG.data = 50 * ones(3, 10000);
    EEG.data(1:2, 3000) = 60; EEG.data(1:2, 7000) = 70; EEG.data(3, [3000 7000]) = 1e4;
    s = RELAX_screen_blink_residual(eeg_checkset(EEG), [3000 7000], {'E1', 'E2'});
    assert(max(abs(s - [10 20])) < 1e-9, 'offset removed by pre-blink baseline; E3 ignored');
    assert(all(isnan(RELAX_screen_blink_residual(eeg_checkset(EEG), [3000 7000], {'E9'}))), 'no blink electrode present');
    fprintf('test_RELAX_screen_blink_residual: PASS\n');
end
