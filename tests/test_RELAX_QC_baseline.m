function test_RELAX_QC_baseline()
% The rebuilt pre-wICA baseline must line up sample for sample with the cleaned file,
% so baseline − cleaned is exactly what wICA removed; any misalignment must be reported.

    rng(1);
    extremesEEG = makeExtremesFixture([1001 1500; 4001 4200]);
    [cleanedEEG, expectedRemoved] = wrapperTail(extremesEEG);

    [baseline, info] = RELAX_QC_baseline(extremesEEG, cleanedEEG);
    assert(info.aligned, 'fixture should align: %s', info.reason);
    assert(info.reason == "", 'no reason when aligned');
    assert(isequal({baseline.chanlocs.labels}, {cleanedEEG.chanlocs.labels}), 'channel order follows the cleaned file');
    kept = ~isnan(expectedRemoved(:, 1));
    removed = double(baseline.data) - double(cleanedEEG.data);
    assert(max(abs(removed(kept, :) - expectedRemoved(kept, :)), [], 'all') < 1e-3, 'baseline − cleaned must equal what was removed');

    % Cleaned file with its channels in a different order still aligns, reordered to match
    perm = randperm(cleanedEEG.nbchan);
    shuffled = cleanedEEG; shuffled.data = cleanedEEG.data(perm, :); shuffled.chanlocs = cleanedEEG.chanlocs(perm);
    [shuffledBaseline, info] = RELAX_QC_baseline(extremesEEG, shuffled);
    assert(info.aligned, 'channel order alone is not a misalignment');
    assert(max(abs(double(shuffledBaseline.data) - double(baseline.data(perm, :))), [], 'all') < 1e-9, 'rows follow the reordering');

    % Misalignments are reported, not silently accepted
    truncated = pop_select(cleanedEEG, 'point', [1 cleanedEEG.pnts - 1]);
    [~, info] = RELAX_QC_baseline(extremesEEG, truncated);
    assert(~info.aligned && contains(info.reason, 'pnts'), 'sample-count mismatch: %s', info.reason);

    dropped = pop_select(cleanedEEG, 'nochannel', {'E1'});
    [~, info] = RELAX_QC_baseline(extremesEEG, dropped);
    assert(~info.aligned && contains(info.reason, 'channel'), 'channel-set mismatch: %s', info.reason);

    moved = cleanedEEG;
    b = find(strcmp({moved.event.type}, 'boundary'), 1);
    moved.event(b).latency = moved.event(b).latency + 50;
    [~, info] = RELAX_QC_baseline(extremesEEG, moved);
    assert(~info.aligned && contains(info.reason, 'boundary'), 'boundary mismatch: %s', info.reason);

    % A file with no extreme periods skips the excision
    extremesNone = makeExtremesFixture(zeros(0, 2));
    [cleanedNone, ~] = wrapperTail(extremesNone);
    [~, info] = RELAX_QC_baseline(extremesNone, cleanedNone);
    assert(info.aligned, 'no extreme periods: %s', info.reason);

    fprintf('test_RELAX_QC_baseline: PASS\n');
end

function EEG = makeExtremesFixture(periods)
    nChan = 32; srate = 500; pnts = 6000;
    EEG = eeg_emptyset();
    EEG.srate = srate; EEG.nbchan = nChan; EEG.pnts = pnts; EEG.trials = 1; EEG.xmin = 0;
    EEG.data = 10 * randn(nChan, pnts);
    EEG.chanlocs = make_qc_chanlocs(nChan);
    EEG.event = struct('type', {'boundary'}, 'latency', {3000.5}, 'duration', {0});
    EEG = eeg_checkset(EEG);
    fullMontage = EEG.chanlocs;
    EEG = pop_select(EEG, 'nochannel', {'E5', 'E17'});
    EEG.allchan = fullMontage;
    EEG.RELAX.Data_has_been_averagerereferenced = 0;
    EEG.RELAX.Data_has_been_cleaned = 0;
    EEG.RELAX.ListOfChannelsAfterRejections = {EEG.chanlocs.labels};
    EEG.RELAX.ExtremelyBadPeriodsForDeletion = periods;
end

function [cleanedEEG, expectedRemoved] = wrapperTail(extremesEEG)
% RELAX_Wrapper_slurm.m between the Extremes_Rejected save and the cleaned save, with
% random data standing in for what wICA removes.
    EEG = extremesEEG;
    if ~isempty(EEG.RELAX.ExtremelyBadPeriodsForDeletion)
        EEG = eeg_eegrej(EEG, EEG.RELAX.ExtremelyBadPeriodsForDeletion);
    end
    EEG = RELAX_average_rereference(EEG);
    wicaRemoved = 5 * randn(EEG.nbchan, EEG.pnts);
    EEG.data = double(EEG.data) - wicaRemoved;
    preInterpLabels = {EEG.chanlocs.labels};
    EEG = pop_interp(EEG, EEG.allchan, 'spherical');
    cleanedEEG = EEG;
    [~, where] = ismember(preInterpLabels, {cleanedEEG.chanlocs.labels});
    expectedRemoved = nan(cleanedEEG.nbchan, cleanedEEG.pnts);
    expectedRemoved(where, :) = wicaRemoved;
end
