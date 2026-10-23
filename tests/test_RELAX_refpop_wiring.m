function test_RELAX_refpop_wiring()
% The wrapper's reference-pop sub-stage is a straight call sequence; this pins the
% contract it relies on: detector -> subtraction -> covariates, on a synthetic EEG, and
% that the auxiliary electrodes receive the same correction.
    srate = 1000; n = 20000; rng(2); N = 64;
    EEG = struct('srate', srate, 'pnts', n, 'trials', 1, 'xmin', 0, 'xmax', (n-1)/srate, ...
                 'event', struct('type', {}, 'latency', {}));
    EEG.chanlocs = struct('labels', arrayfun(@(k) sprintf('E%d', k), 1:N, 'UniformOutput', false));
    r = 80 ./ (1 + exp(-((1:n) - 8000) / 3.3));
    EEG.data = 10 * randn(N, n) - r; EEG.nbchan = N;
    aux = EEG; aux.data = EEG.data(1:2, :); aux.nbchan = 2; aux.chanlocs = EEG.chanlocs(1:2);

    [refSteps, refInfo] = RELAX_detect_reference_pops(EEG, 'Detect', {'MinAmplitudeMicrovolts', 30, 'SlopeNoiseFloorMicrovolts', 0.1});
    assert(refInfo.nSteps == 1);
    [EEG, repInfo] = RELAX_subtract_reference_pops(EEG, refSteps, refInfo.refChannel.data, 'WindowMs', refInfo.detectInfo.params.WindowMs);
    [aux, auxInfo] = RELAX_subtract_reference_pops(aux, refSteps, refInfo.refChannel.data, 'WindowMs', refInfo.detectInfo.params.WindowMs);
    assert(repInfo.available && auxInfo.available && isequal(aux.data, EEG.data(1:2, :)));

    cov = RELAX_refpop_covariates(refInfo, repInfo);
    for f = {'RefPopNCandidates','RefPopNSteps','RefPopNRejectedUniformity','RefPopNMerged','RefPopStepsPerMinute', ...
             'RefPopMedianAmplitude','RefPopMedianRiseMs','RefPopMedianSignAgreement', ...
             'RefPopMedianAmplitudeSpread','RefPopTransitionSeconds','RefPopNChannelsUsed','RefPopChannelRobustSd'}
        assert(isfield(cov, f{1}) && isnumeric(cov.(f{1})) && isscalar(cov.(f{1})), '%s', f{1});
    end
    assert(cov.RefPopNSteps == 1 && cov.RefPopNCandidates == 1 && cov.RefPopNRejectedUniformity == 0 && cov.RefPopNMerged == 0);
    assert(abs(cov.RefPopMedianAmplitude - 80) < 12 && cov.RefPopTransitionSeconds > 0);

    covEmpty = RELAX_refpop_covariates(struct('available', false, 'status', 'x'), struct('available', false));
    assert(isnan(covEmpty.RefPopMedianAmplitude) && covEmpty.RefPopNSteps == 0);
    fprintf('test_RELAX_refpop_wiring: PASS\n');
end