function test_RELAX_QC_blink_peaks()
% Only blinks whose whole ±window sits inside one segment may enter a blink-locked
% average; an epoch spanning a boundary would mix data the pipeline joined after excision.

    EEG.pnts = 2000;
    EEG.event = struct('type', {'EyeBlinkMax', 'EyeBlinkMax', 'boundary', 'EyeBlinkMax', 'EyeBlinkMax', 7, 'EyeBlinkLeftBase'}, ...
                       'latency', {50, 500.4, 1000.5, 995, 1500, 1200, 1600});
    peaks = RELAX_QC_blink_peaks(EEG, 100);
    assert(isequal(peaks, [500 1500]), 'expected [500 1500], got %s', mat2str(peaks));

    EEG.event = struct('type', {}, 'latency', {});
    assert(isempty(RELAX_QC_blink_peaks(EEG, 100)), 'no events, no peaks');

    fprintf('test_RELAX_QC_blink_peaks: PASS\n');
end
