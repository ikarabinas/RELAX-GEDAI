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

    assert(OUT.RELAX.Data_has_been_cleaned == 1, 'Data_has_been_cleaned must be 1');
    assert(~any(isnan(OUT.data(:))), 'cleaned data must contain no NaNs');

    fprintf('Retained channels: %d\n', OUT.nbchan);
    fprintf('assert_smoke_outcome: PASS\n');
end
