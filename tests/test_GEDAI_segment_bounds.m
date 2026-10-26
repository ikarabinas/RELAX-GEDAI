function test_GEDAI_segment_bounds()
% GEDAI_segment_bounds (GEDAI fork, auxiliaries/) turns boundary events into the
% [start end] sample ranges GEDAI processes independently. A boundary only splits the
% data if it falls between two samples that exist. eeg_eegrej writes a boundary at
% latency 0.5 when it deletes from sample 1; GEDAI used to clamp that to sample 1 and
% process sample 1 as its own one-sample segment, returning it at ~36x its input.

    addpath('/home/nga4004/gitRepos/GEDAI-master/auxiliaries');

    EEG = struct('pnts', 100, 'event', struct('type', {}, 'latency', {}));
    assert(isequal(GEDAI_segment_bounds(EEG), [1 100]), 'no boundaries: one segment');

    EEG.event = struct('type', {'boundary'}, 'latency', {0.5});
    assert(isequal(GEDAI_segment_bounds(EEG), [1 100]), ...
        'boundary at 0.5 lies before the first sample and must not split off sample 1');

    EEG.event = struct('type', {'boundary'}, 'latency', {100.5});
    assert(isequal(GEDAI_segment_bounds(EEG), [1 100]), ...
        'boundary after the last sample must not split anything');

    EEG.event = struct('type', {'boundary', 'X', 'boundary', 'boundary'}, 'latency', {60.5, 10, 0.5, 30.5});
    assert(isequal(GEDAI_segment_bounds(EEG), [1 30; 31 60; 61 100]), ...
        'interior x.5 boundaries split between samples; non-boundary events and the 0.5 boundary are ignored');

    EEG.event = struct('type', {'boundary', 'boundary'}, 'latency', {30.5, 30.5});
    assert(isequal(GEDAI_segment_bounds(EEG), [1 30; 31 100]), 'duplicate boundaries split once');

    EEG.event = struct('type', {'boundary', 'boundary'}, 'latency', {30.5, 31.5});
    assert(isequal(GEDAI_segment_bounds(EEG), [1 30; 31 31; 32 100]), ...
        'a genuine one-sample segment between two adjacent boundaries is kept');

    EEG.event = struct('type', {'boundary'}, 'latency', {40});
    assert(isequal(GEDAI_segment_bounds(EEG), [1 40; 41 100]), ...
        'integer boundary latency ends the segment at that sample');

    fprintf('test_GEDAI_segment_bounds: PASS\n');
end
