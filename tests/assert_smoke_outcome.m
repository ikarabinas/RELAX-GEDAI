function assert_smoke_outcome(cleanedFile, fixtureFile)
% Assert the boundary-aware smoke run's outcome. Separated from run_integration_smoke
% so it can be re-run against an existing cleaned file without repeating the pipeline.

    OUT = pop_loadset(cleanedFile);
    FIX = pop_loadset(fixtureFile);

    assert(isfield(OUT,'RELAXProcessingExtremeRejections'), 'missing RELAXProcessingExtremeRejections');
    assert(isempty(OUT.RELAXProcessingExtremeRejections.PREPBasedChannelToReject), ...
        'PREPBasedChannelToReject must be empty when PREP is skipped');

    assert(~isfield(OUT.RELAX,'NoiseMaskFullLengthR1'), ...
        'NoiseMaskFullLengthR1 must be absent when ComputeSER_ARR=0');

    % Boundaries must survive cleaning. eeg_eegrej can both add boundaries (rejecting
    % extreme periods) and merge them (when a rejected period spans a whole segment), so
    % exact preservation is not expected - but a near-total collapse means something went
    % badly wrong, and comparing against 1 would not catch that.
    fixSegments = RELAX_get_segments(FIX);
    outSegments = RELAX_get_segments(OUT);
    fprintf('Segments: fixture %d -> cleaned %d.\n', size(fixSegments,1), size(outSegments,1));
    assert(size(outSegments,1) >= ceil(size(fixSegments,1)/2), ...
        'cleaned file kept only %d of the fixture''s %d segments', ...
        size(outSegments,1), size(fixSegments,1));

    % Pop stage (when it ran): the info structs and every scalar covariate must be present,
    % and the repair must have stayed within its per-channel ceiling.
    if isfield(OUT.RELAX, 'PopScreen')
        assert(OUT.RELAX.PopScreen.available, 'pop screen reported unavailable: %s', OUT.RELAX.PopScreen.status);
        X = OUT.RELAXProcessingExtremeRejections;
        needed = {'PopNSteps','PopNChannelsAffected','PopNEvents','PopNChannelsRejected', ...
                  'PopNChannelsRejectedByRate','PopNChannelsRejectedByFraction', ...
                  'PopNChannelsRejectedByFlat','PopNEyeElectrodesProtected', ...
                  'PopNCommonModeSteps','PopNCommonModeEvents', ...
                  'PopNEyeMovementStepsSkipped','PopEyeMaskFraction','PopNIntervalsFilled', ...
                  'PopFractionSamplesInterpolated','PopMaxChannelFractionInterpolated', ...
                  'PopNFillsSpatial','PopNFillsLinear','PopNCommonModeExcised','PopExcisedSeconds'};
        for k = 1:numel(needed)
            assert(isfield(X, needed{k}), 'pop covariate %s missing', needed{k});
            assert(isscalar(X.(needed{k})) && isnumeric(X.(needed{k})), 'pop covariate %s is not a numeric scalar', needed{k});
        end
        refNeeded = {'RefPopNCandidates','RefPopNSteps','RefPopNRejectedUniformity','RefPopNMerged', ...
                     'RefPopStepsPerMinute','RefPopMedianAmplitude','RefPopMedianRiseMs', ...
                     'RefPopMedianSignAgreement','RefPopMedianAmplitudeSpread', ...
                     'RefPopTransitionSeconds','RefPopNChannelsUsed','RefPopChannelRobustSd'};
        for f = refNeeded
            assert(isfield(X, f{1}), 'missing reference-pop covariate %s', f{1});
            assert(isnumeric(X.(f{1})) && isscalar(X.(f{1})), '%s must be a numeric scalar', f{1});
        end
        assert(isfield(OUT.RELAX, 'RefPopScreen') && OUT.RELAX.RefPopScreen.available, ...
            'reference-pop screen must report available on the fixture');
        assert(~isfield(OUT.RELAX.RefPopScreen, 'refChannel'), 'the virtual channel trace must not be saved with the .set');
        assert(X.RefPopNSteps + X.RefPopNRejectedUniformity + X.RefPopNMerged == X.RefPopNCandidates, 'candidate bookkeeping');
        fprintf('Reference-pop stage: %d candidates, %d subtracted, %d rejected by uniformity, %d merged.\n', ...
            X.RefPopNCandidates, X.RefPopNSteps, X.RefPopNRejectedUniformity, X.RefPopNMerged);
        assert(OUT.RELAX.PopRepair.available, 'pop repair reported unavailable: %s', OUT.RELAX.PopRepair.status);
        % Blink-electrode fallback covariates (RELAX_blinks_IQR_method) are always written,
        % whether or not the fallback fired.
        for k = {'BlinkElectrodesUsed','BlinkElectrodeFallback'}
            assert(isfield(X, k{1}) && isscalar(X.(k{1})) && isnumeric(X.(k{1})), ...
                'blink covariate %s missing or not a numeric scalar', k{1});
        end
        assert(any(X.BlinkElectrodeFallback == [0 1 2]), 'BlinkElectrodeFallback must be 0, 1 or 2');
        assert(X.PopMaxChannelFractionInterpolated <= 0.05, ...
            'a channel had %.1f%% of its samples interpolated', 100*X.PopMaxChannelFractionInterpolated);
        fprintf('Pop stage: %d steps, %d filled (%d spatial), %d excised, %d eye-skipped, %d channels rejected.\n', ...
            X.PopNSteps, X.PopNIntervalsFilled, X.PopNFillsSpatial, X.PopNCommonModeExcised, ...
            X.PopNEyeMovementStepsSkipped, X.PopNChannelsRejected);
    end

    assert(OUT.RELAX.Data_has_been_cleaned == 1, 'Data_has_been_cleaned must be 1');
    assert(~any(isnan(OUT.data(:))), 'cleaned data must contain no NaNs');

    fprintf('Retained channels: %d\n', OUT.nbchan);
    fprintf('assert_smoke_outcome: PASS\n');
end
