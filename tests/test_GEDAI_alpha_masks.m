function test_GEDAI_alpha_masks()
% GEDAI_alpha_masks decides which samples the GEDAI alpha figures may count: more than
% edgeSec from a boundary or file edge, in segments of at least minSegSec, with a rest
% block defined by the TMS_train_* events. GEDAI_alpha_mask_psd turns a mask into runs
% of true and passes them to RELAX_QC_segment_psd so no Welch window spans a gap.

    opts = GEDAI_alpha_figures('opts');
    srate = 250;
    nSamples = 60 * srate;

    boundaryTimes = [0.5 2500.5 2501.5];
    M = GEDAI_alpha_masks(fixtureWithEvents(srate, nSamples, boundaryTimes, {}, []), opts);
    assert(~M.valid(2501), '(a) the one-sample segment at 2501 must be invalid');
    assert(M.nTinySegments >= 1, '(a) the one-sample segment must be counted as tiny');
    assert(~any(M.valid(1:75)), '(a) the first 0.3 s (75 samples) must be invalid');
    assert(M.valid(7500), '(a) a sample far from any boundary must be valid');
    assert(M.nExcluded == sum(~M.valid), '(a) nExcluded must count the invalid samples');

    shortSegment = 3001:3375;
    assert(numel(shortSegment) == 1.5 * srate);
    M = GEDAI_alpha_masks(fixtureWithEvents(srate, nSamples, [3000.5 3375.5], {}, []), opts);
    assert(~any(M.valid(shortSegment)), '(b) a 1.5 s segment is shorter than minSegSec and must be entirely invalid');
    assert(M.valid(7500), '(b) long segments are unaffected');

    M = GEDAI_alpha_masks(fixtureWithEvents(srate, nSamples, [], {}, []), opts);
    assert(M.restIsWholeFile, '(c) no TMS_train events means the rest block is the whole file');
    assert(isequal(M.rest, M.valid), '(c) rest must equal valid when there are no trains');

    trainTypes = {'TMS_train_1', 'TMS_train_2'};
    trainLatencies = [10 20] * srate;
    M = GEDAI_alpha_masks(fixtureWithEvents(srate, nSamples, [], trainTypes, trainLatencies), opts);
    assert(~M.restIsWholeFile, '(d) trains present, so the rest block is not the whole file');
    assert(~any(M.rest(10*srate : 28*srate)), '(d) rest must be false from the first train to 8 s after the last');
    assert(M.rest(5*srate), '(d) rest must be true before the first train');
    assert(M.rest(40*srate), '(d) rest must be true after the last train plus 8 s');
    assert(all(M.rest <= M.valid), '(d) rest must be ANDed with valid');

    rng = RandStream('mt19937ar', 'Seed', 1);
    t = (0:nSamples-1) / srate;
    data = 0.1 * randn(rng, 4, nSamples);
    data(1, :) = data(1, :) + sin(2 * pi * 10 * t);
    gapStart = 29 * srate + 1;
    gapEnd = 32 * srate;
    mask = true(1, nSamples);
    mask(gapStart:gapEnd) = false;

    [f, P] = GEDAI_alpha_mask_psd(data, srate, mask, opts.winSec, opts.maxHz);
    assert(size(P, 1) == 4 && size(P, 2) == numel(f), '(e) P must be nch x nf');
    [~, peak] = max(P(1, :));
    assert(abs(f(peak) - 10) <= 0.5, '(e) peak of the 10 Hz sine is at %g Hz', f(peak));

    stepped = data;
    stepMiddle = round((gapStart + gapEnd) / 2);
    stepped(:, stepMiddle:end) = stepped(:, stepMiddle:end) + 20000;
    [~, Pstep] = GEDAI_alpha_mask_psd(stepped, srate, mask, opts.winSec, opts.maxHz);
    lowBand = f >= 1 & f <= 4;
    withoutStep = sum(P(1, lowBand));
    withStep = sum(Pstep(1, lowBand));
    assert(abs(withStep - withoutStep) <= 0.05 * withoutStep, ...
        '(e) a 20 mV step inside the gap changed 1-4 Hz power from %g to %g', withoutStep, withStep);

    fprintf('test_GEDAI_alpha_masks: PASS\n');
end

function EEG = fixtureWithEvents(srate, nSamples, boundaryLatencies, trainTypes, trainLatencies)
    EEG = eeg_emptyset();
    EEG.srate = srate;
    EEG.nbchan = 4;
    EEG.pnts = nSamples;
    EEG.trials = 1;
    EEG.xmin = 0;
    EEG.xmax = (nSamples - 1) / srate;
    EEG.data = zeros(4, nSamples);
    EEG.times = (0:nSamples-1) * 1000 / srate;
    types = [repmat({'boundary'}, 1, numel(boundaryLatencies)), trainTypes];
    latencies = num2cell([boundaryLatencies(:)', trainLatencies(:)']);
    EEG.event = struct('type', types, 'latency', latencies);
end
