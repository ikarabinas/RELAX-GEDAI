function verify_epoching_drops_boundary_epochs()
% Verify the assumption that pop_epoch discards every epoch containing a boundary
% event, which is what protects RELAX_epoching, RELAX_excluding_channels_and_epoching,
% RELAX_excluding_extreme_values, RELAX_muscle and RELAX_drift from discontinuities.
%
% Requires EEGLAB on the path and tests/fixtures/boundary_fixture.set to exist.

    fixturePath = fullfile(fileparts(mfilename('fullpath')), 'fixtures', 'boundary_fixture.set');
    assert(exist(fixturePath,'file') == 2, 'run make_boundary_fixture first');

    EEG = pop_loadset(fixturePath);
    segments = RELAX_get_segments(EEG);
    fprintf('Fixture has %d segments.\n', size(segments,1));
    assert(size(segments,1) > 1, 'fixture must contain boundaries for this check to mean anything');

    % Reproduce the epoching that RELAX_epoching performs
    EEG = eeg_regepochs(EEG,'recurrence',0.500,'eventtype','X','extractepochs','off');
    for e = 1:numel(EEG.event)
        EEG.event(e).originallatency = EEG.event(e).latency;
    end
    EEG = pop_selectevent(EEG, 'type', 'X', 'deleteevents', 'off');
    epoched = pop_epoch(EEG, {'X'}, [0 1.0], 'epochinfo', 'yes');

    % No surviving epoch may contain a boundary event
    boundaryInEpoch = 0;
    for ep = 1:numel(epoched.epoch)
        types = epoched.epoch(ep).eventtype;
        if ~iscell(types); types = {types}; end
        if any(strcmp(types, 'boundary'))
            boundaryInEpoch = boundaryInEpoch + 1;
        end
    end
    fprintf('Epochs retained: %d. Epochs containing a boundary: %d.\n', epoched.trials, boundaryInEpoch);
    assert(boundaryInEpoch == 0, ...
        'pop_epoch retained %d epoch(s) containing a boundary - the epoch-based stages are NOT protected and need their own fix', ...
        boundaryInEpoch);

    % Every retained epoch must lie wholly inside one segment
    epochLen = round(1.0 * epoched.srate);
    outside = 0;
    for ep = 1:numel(epoched.epoch)
        lat = epoched.epoch(ep).eventoriginallatency;
        if iscell(lat); lat = lat{1}; end
        inSegment = any(segments(:,1) <= lat & segments(:,2) >= lat + epochLen - 1);
        if ~inSegment
            outside = outside + 1;
        end
    end
    fprintf('Epochs not wholly inside a single segment: %d.\n', outside);
    assert(outside == 0, '%d epoch(s) straddle a segment boundary', outside);

    fprintf('verify_epoching_drops_boundary_epochs: PASS\n');
end
