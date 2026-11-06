function test_RELAX_interpolate_channel_epochs()
% A flagged channel is replaced over its whole flagged epoch by a spherical spline of the
% epoch's unflagged channels, blended in by a raised cosine over 0.25 s inside each edge.
% Nothing else changes.
%
% Fixture: the unflagged channels are all 0, so the spline is 0 and the output of a
% flagged channel holding 100 is exactly (1 - w) * 100.

    onsets = 1:50:451;
    epochSamples = 100;
    w = [0.5 * (1 - cos(pi * (1:25) / 26)), ones(1, 50), 0.5 * (1 + cos(pi * (1:25) / 26))];

    % spike
    EEG = zero_fixture();
    EEG.data(5, 201:300) = 100;
    flagged = false(32, numel(onsets));
    flagged(5, 5) = true;
    [out, changed] = RELAX_interpolate_channel_epochs(EEG, flagged, onsets, epochSamples, 0.25);
    assert(max(abs(out.data(5, 201:300) - (1 - w) * 100)) < 1e-9, 'flagged channel is not the raised-cosine blend of spline and signal');
    others = [1:4 6:32];
    assert(isequal(out.data(others, :), EEG.data(others, :)), 'an unflagged channel changed');
    assert(isequal(out.data(5, [1:200 301:end]), EEG.data(5, [1:200 301:end])), 'the flagged channel changed outside its epoch');
    assert(isequal(find(changed(5, :)), 201:300) && nnz(changed) == 100, 'changed mask is not exactly the flagged epoch');

    % overlap: two flagged epochs sharing 251-300 leave no dip back to the signal
    EEG = zero_fixture();
    EEG.data(5, 201:350) = 100;
    flagged = false(32, numel(onsets));
    flagged(5, [5 6]) = true;
    out = RELAX_interpolate_channel_epochs(EEG, flagged, onsets, epochSamples, 0.25);
    assert(max(abs(out.data(5, 226:325))) < 1e-9, 'overlapping flagged epochs dip back towards the signal');

    % two_channels: both replaced from the other 30
    EEG = zero_fixture();
    EEG.data([5 9], 201:300) = 100;
    flagged = false(32, numel(onsets));
    flagged([5 9], 5) = true;
    out = RELAX_interpolate_channel_epochs(EEG, flagged, onsets, epochSamples, 0.25);
    assert(max(abs(out.data([5 9], 226:275)), [], 'all') < 1e-9, 'two flagged channels were not both replaced from the rest');

    % keeps_class
    EEG = zero_fixture();
    EEG.data = single(EEG.data);
    EEG.data(5, 201:300) = 100;
    flagged = false(32, numel(onsets));
    flagged(5, 5) = true;
    out = RELAX_interpolate_channel_epochs(EEG, flagged, onsets, epochSamples, 0.25);
    assert(isa(out.data, 'single'), 'single data came back as %s', class(out.data));

    % nothing_flagged
    [out, changed] = RELAX_interpolate_channel_epochs(EEG, false(32, numel(onsets)), onsets, epochSamples, 0.25);
    assert(isequal(out.data, EEG.data) && ~any(changed, 'all'), 'nothing flagged, yet something changed');

    fprintf('test_RELAX_interpolate_channel_epochs: PASS\n');
end

function EEG = zero_fixture()
    EEG = eeg_emptyset();
    EEG.srate = 100;
    EEG.data = zeros(32, 600);
    EEG.nbchan = 32; EEG.pnts = 600; EEG.trials = 1;
    EEG.xmin = 0; EEG.xmax = 5.99;
    EEG.chanlocs = make_qc_chanlocs(32);
end
