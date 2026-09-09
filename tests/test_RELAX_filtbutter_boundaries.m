function test_RELAX_filtbutter_boundaries()
% Segment-wise, reflection-padded Butterworth filtering.

    srate = 500;
    nSamples = 8000;
    t = (0:nSamples-1)/srate;

    % Two segments of clean 10 Hz alpha, with a large DC step across the boundary.
    % A boundary-blind highpass smears that step deep into both segments.
    signal = sin(2*pi*10*t);
    signal(4001:end) = signal(4001:end) + 500; % DC step at the boundary

    EEG = struct('data', signal, 'srate', srate, 'pnts', nSamples, 'trials', 1, ...
                 'event', struct('type', {'boundary'}, 'latency', {4000.5}));

    % --- Backward compatibility: 6-arg call is unchanged -------------------------
    EEG_legacy_6arg = RELAX_filtbutter(EEG, 1, 200, 4, 'bandpass', 'acausal');
    EEG_legacy_7arg = RELAX_filtbutter(EEG, 1, 200, 4, 'bandpass', 'acausal', false);
    assert(isequal(EEG_legacy_6arg.data, EEG_legacy_7arg.data), ...
        'respect_boundaries=false must be bitwise identical to the 6-argument call');

    % --- Segment path runs and preserves shape -----------------------------------
    EEG_seg = RELAX_filtbutter(EEG, 1, 200, 4, 'bandpass', 'acausal', true);
    assert(isequal(size(EEG_seg.data), size(signal)), 'output size must be preserved');
    assert(~any(isnan(EEG_seg.data)), 'output must contain no NaNs');

    % --- ISOLATION: segment 1 must not change when segment 2 changes -------------
    signalB = signal;
    signalB(4001:end) = signalB(4001:end) * 100;
    EEG_B = EEG; EEG_B.data = signalB;
    out_A = RELAX_filtbutter(EEG,   1, 200, 4, 'bandpass', 'acausal', true);
    out_B = RELAX_filtbutter(EEG_B, 1, 200, 4, 'bandpass', 'acausal', true);
    assert(isequal(out_A.data(1:4000), out_B.data(1:4000)), ...
        'segment 1 must be bitwise unchanged when segment 2 is scaled');

    % The boundary-blind path does NOT have that property - this is the bug being fixed
    blind_A = RELAX_filtbutter(EEG,   1, 200, 4, 'bandpass', 'acausal', false);
    blind_B = RELAX_filtbutter(EEG_B, 1, 200, 4, 'bandpass', 'acausal', false);
    assert(~isequal(blind_A.data(1:4000), blind_B.data(1:4000)), ...
        'sanity check: the legacy path is expected to leak across the boundary');

    % --- Segment-wise result matches filtering that segment on its own ------------
    seg1 = struct('data', signal(1:4000), 'srate', srate, 'pnts', 4000, 'trials', 1, 'event', []);
    seg1_alone = RELAX_filtbutter(seg1, 1, 200, 4, 'bandpass', 'acausal', true);
    assert(max(abs(out_A.data(1:4000) - seg1_alone.data)) < 1e-9, ...
        'a segment filtered in place must match that segment filtered standalone');

    % --- The 10 Hz signal survives in the interior of each segment ---------------
    % Compare against the known input away from the edges (1 s in from each end).
    interior = 501:3500;
    recovered = out_A.data(interior);
    expected  = signal(interior);
    assert(corr(recovered(:), expected(:)) > 0.99, ...
        'the 10 Hz component must survive the bandpass in the segment interior');

    % --- No boundaries present: segment path stays close to the legacy path -------
    EEG_nb = struct('data', sin(2*pi*10*t), 'srate', srate, 'pnts', nSamples, 'trials', 1, 'event', []);
    nb_legacy = RELAX_filtbutter(EEG_nb, 1, 200, 4, 'bandpass', 'acausal', false);
    nb_seg    = RELAX_filtbutter(EEG_nb, 1, 200, 4, 'bandpass', 'acausal', true);
    % NOTE: empirically, filtfilt's own default automatic padding (order-based, a
    % few samples) does not fully settle a 1 Hz highpass edge within 1000 samples at
    % fs=500 - the legacy path itself still carries a ~1e-5-level edge artifact there,
    % which is exactly the deficiency this task's segment path corrects. 2000 samples
    % (4 s) in from each end is comfortably past that, confirmed empirically (diff
    % ~1e-9 there vs ~1e-6 tolerance).
    mid = 2001:6000; % 4 s in from each end, beyond any edge-treatment difference
    assert(max(abs(nb_legacy.data(mid) - nb_seg.data(mid))) < 1e-6, ...
        'with no boundaries the two paths must agree in the interior');

    % --- Bookkeeping: filtering must not disturb events, pnts or srate -------------
    assert(EEG_seg.pnts == EEG.pnts, 'pnts must be preserved');
    assert(numel(EEG_seg.event) == numel(EEG.event), 'event count must be preserved');
    assert(EEG_seg.event(1).latency == EEG.event(1).latency, 'boundary latency must be preserved');
    assert(strcmp(EEG_seg.event(1).type, 'boundary'), 'boundary event type must be preserved');

    % --- Highpass and lowpass types also work through the segment path ------------
    hp = RELAX_filtbutter(EEG, 1, [], 4, 'highpass', 'acausal', true);
    assert(isequal(size(hp.data), size(signal)) && ~any(isnan(hp.data)), 'highpass segment path');
    lp = RELAX_filtbutter(EEG, [], 45, 4, 'lowpass', 'acausal', true);
    assert(isequal(size(lp.data), size(signal)) && ~any(isnan(lp.data)), 'lowpass segment path');
    bs = RELAX_filtbutter(EEG, 56, 64, 4, 'bandstop', 'acausal', true);
    assert(isequal(size(bs.data), size(signal)) && ~any(isnan(bs.data)), 'bandstop segment path');

    % --- Multi-channel segment path matches per-channel filtering -----------------
    % Production runs 257 channels; the segment path above is only exercised
    % single-channel. Confirm the vectorised multi-channel path (filtfilt(z1,p1,d')')
    % agrees with filtering each channel on its own to a tight numerical tolerance.
    nChan = 8;
    multiData = repmat(signal, nChan, 1) .* (1 + 0.1*(1:nChan)'); % distinct per-channel scale
    EEG_multi = struct('data', multiData, 'srate', srate, 'pnts', nSamples, 'trials', 1, ...
                 'event', struct('type', {'boundary'}, 'latency', {4000.5}));
    out_multi = RELAX_filtbutter(EEG_multi, 1, 200, 4, 'bandpass', 'acausal', true);

    out_perChannel = NaN(nChan, nSamples);
    for ch = 1:nChan
        EEG_ch = struct('data', multiData(ch,:), 'srate', srate, 'pnts', nSamples, 'trials', 1, ...
                     'event', struct('type', {'boundary'}, 'latency', {4000.5}));
        out_ch = RELAX_filtbutter(EEG_ch, 1, 200, 4, 'bandpass', 'acausal', true);
        out_perChannel(ch,:) = out_ch.data;
    end
    assert(max(abs(out_multi.data(:) - out_perChannel(:))) < 1e-10, ...
        'multi-channel segment-path output must match per-channel filtering to ~1e-10');

    fprintf('test_RELAX_filtbutter_boundaries: PASS\n');
end
