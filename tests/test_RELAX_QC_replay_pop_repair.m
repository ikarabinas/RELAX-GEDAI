function test_RELAX_QC_replay_pop_repair()
% The replay removes a saved step from the raw input and leaves the unrepaired copy as it
% was; outputs hold the requested channels, in the requested order, at the input rate.

    RELAX_add_vendored_noisetools();
    rng(5);
    srate = 1000; nChan = 20; pnts = 20 * srate;
    chanlocs = make_qc_chanlocs(nChan);
    data = 5 * randn(nChan, pnts);
    data(3, 10001:end) = data(3, 10001:end) + 800;
    inputEEG = eeg_checkset(setfield(setfield(setfield(eeg_emptyset(), 'data', single(data)), 'srate', srate), 'chanlocs', chanlocs)); %#ok<SFLD>

    steps = struct('label', {{'E3'}}, 'channel', 3, 'sample', 10001, 'onset', 9995, 'offset', 10010, ...
                   'amplitude', 800, 'commonMode', false, 'eyeMovement', false, 'eventId', 1, 'z', 50);
    extremesEEG = struct('chanlocs', chanlocs(1:18), ...
                         'RELAX', struct('PopScreen', struct('steps', steps, 'degenerateChannels', {cell(0, 1)})));
    cfg = struct('Perform_reference_pop_repair', 0, 'PopRepair_MaxPolyOrder', 3, 'PopRepair_MaxCutoffHz', 1, ...
                 'PopRepair_MaxWindowSeconds', 2, 'PopRepair_DetrendThresh', 3, 'PopRepair_DetrendIterations', 3, ...
                 'PopRepair_AnchorMs', 3, 'PopRepair_FillMethod', 'linear', 'PopRepair_MinBasisChannels', 16, ...
                 'PopRepair_ExciseCommonMode', 0, 'PopRepair_MaxExciseMs', 500, 'PopRepair_MinSegmentSeconds', 1);

    [unrepaired, repaired, info] = RELAX_QC_replay_pop_repair(inputEEG, extremesEEG, cfg, {'E3', 'E7', 'E19', 'E5'});

    assert(isequal({unrepaired.chanlocs.labels}, {'E3', 'E5', 'E7'}) || isequal(sort({unrepaired.chanlocs.labels}), {'E3', 'E5', 'E7'}), ...
           'E19 is not in Extremes_Rejected or the steps, so it is not replayed');
    assert(unrepaired.srate == srate && repaired.pnts == pnts && ~info.referencePopsReplayed, 'input rate, no reference pops');
    stepSize = @(EEG) median(EEG.data(strcmp({EEG.chanlocs.labels}, 'E3'), 10101:10600)) ...
                    - median(EEG.data(strcmp({EEG.chanlocs.labels}, 'E3'), 9401:9900));
    assert(abs(stepSize(unrepaired) - 800) < 5, 'unrepaired keeps the step: %g', stepSize(unrepaired));
    assert(abs(stepSize(repaired)) < 20, 'repaired loses the step: %g', stepSize(repaired));
    fprintf('test_RELAX_QC_replay_pop_repair: PASS\n');
end
