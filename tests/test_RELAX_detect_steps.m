function test_RELAX_detect_steps()
% Two-window step detector: recovers slow-rising pops with their transition interval,
% rejects drift/noise/oscillation, groups co-occurring steps into events, and ranks bad
% contacts. Every section uses function defaults unless it says otherwise.

    srate = 1000;
    n     = 10000;          % 10 s
    rng(7);

    % ---- helpers ----
    function EEG = mkEEG(data, boundaries)
        if nargin < 2; boundaries = []; end
        EEG = struct();
        EEG.data     = data;
        EEG.nbchan   = size(data,1);
        EEG.pnts     = size(data,2);
        EEG.srate    = srate;
        EEG.trials   = 1;
        EEG.event    = struct('type', {}, 'latency', {});
        for b = boundaries(:)'
            EEG.event(end+1) = struct('type', 'boundary', 'latency', b);
        end
        EEG.chanlocs = struct('labels', ...
            arrayfun(@(k) sprintf('E%d', k), 1:size(data,1), 'UniformOutput', false));
    end
    numericCols = {'channel','sample','onset','offset','amplitude','z','channelZ', ...
                   'riseMs','eventId'};

    % The operating point section 3 was written against: its bracketing shapes need the long
    % TransitionMs of the superseded defaults (2026-09-23).
    legacy = {'WindowMs', 150, 'TransitionMs', 100, 'ZThreshold', 12, 'ChannelZThreshold', 10, ...
              'MinAmplitudeMicrovolts', 50, 'RefractoryMs', 20, 'SlopeBaselineMs', 2000, ...
              'CommonModeChannelFraction', 0.08};

    %% ---- 1a. empty result has 0x1 columns of the right class ----
    EEG = mkEEG(randn(2, n) * 1e-9);
    [steps, info] = RELAX_detect_steps(EEG);
    assert(info.available, 'quiet data is still available: %s', info.status);
    for c = numericCols
        assert(isnumeric(steps.(c{1})) && isequal(size(steps.(c{1})), [0 1]), ...
            'steps.%s must be a 0x1 numeric column on an empty result', c{1});
    end
    assert(iscell(steps.label) && isequal(size(steps.label), [0 1]), 'label must be 0x1 cell');
    assert(islogical(steps.commonMode) && isequal(size(steps.commonMode), [0 1]), ...
        'commonMode must be a 0x1 logical column');
    assert(info.nSteps == 0 && info.nChannelsAffected == 0 && info.nCommonModeSteps == 0 ...
        && info.nCommonModeEvents == 0);
    for f = {'badChannels','badChannelReason','degenerateChannels'}
        assert(iscell(info.(f{1})) && isequal(size(info.(f{1})), [0 1]), ...
            'info.%s must be a 0x1 cell', f{1});
    end
    for f = {'stepsPerChannel','stepRatePerMinute','transitionFraction','channelScaleSlope'}
        assert(isequal(size(info.(f{1})), [EEG.nbchan 1]), 'info.%s must be nbchan x 1', f{1});
    end
    assert(~isfield(info, 'unstableFraction'), ...
        'the sample-level instability criterion was retired on 2026-09-22');
    assert(isnan(info.medianRiseMs), 'medianRiseMs is NaN when nothing was found');
    assert(isstruct(info.events) || istable(info.events));

    %% ---- 1b. an exactly flat channel is degenerate, not an error ----
    EEG = mkEEG([zeros(1, n); randn(1, n) * 5]);
    [steps, info] = RELAX_detect_steps(EEG);
    assert(any(strcmp(info.degenerateChannels, 'E1')), 'flat channel must be listed as degenerate');
    assert(isequal(size(info.degenerateChannels), [1 1]));
    assert(~any(isnan(steps.z)), 'no NaN z may leak out of a degenerate channel');
    assert(isempty(steps.sample(steps.channel == 1)));

    % a railed channel: constant within EVERY inter-boundary segment, at different levels
    % (m544 day 2, E154: 0 uV in some segments, +333333 uV in the others, 0 steps). The
    % global range is not zero, but the channel carries nothing and the repair's per-segment
    % demean would turn it into exact zeros, which crashes pop_jointprob downstream.
    rs = RandStream('twister', 'Seed', 11);             % own stream: leave rng(7) untouched for later sections
    x = [zeros(1, 5000), 333333 * ones(1, 5000)];
    [steps, info] = RELAX_detect_steps(mkEEG([x; randn(rs, 1, n) * 5], 5000.5));
    assert(isequal(info.degenerateChannels, {'E1'}), 'per-segment constant channel must be degenerate');
    assert(isempty(steps.sample(steps.channel == 1)), 'no steps are reported on a degenerate channel');
    % constant in one segment only is NOT degenerate (the other segment is real data)
    x = [zeros(1, 5000), randn(rs, 1, 5000) * 5];
    [~, info] = RELAX_detect_steps(mkEEG([x; randn(rs, 1, n) * 5], 5000.5));
    assert(isempty(info.degenerateChannels), 'a channel flat in only one segment is not degenerate');

    %% ---- 1c. the gate scores how fast the median level moves, not how far ----
    % Fixtures measured against the pinned level gate on 2026-09-22 (jobs 19866817, 19866836).
    % The drift envelope is tapered because a rectangular onset is itself a fast level change,
    % which the gate is right to find - the claim under test is about sustained drift.
    rng(7);
    n60c = 60 * srate; t60c = 0:n60c-1;
    quietc = randn(1, n60c) * 5;
    atc = 25000; Lc = 10000;
    envc   = sin(pi*(t60c - atc)/Lc).^2 .* (t60c >= atc & t60c < atc + Lc);
    phasec = 2*pi*(t60c - atc)/srate;

    % a 2 Hz slow wave is delta-band EEG, not 25 electrode pops
    x = quietc + 100*sin(2*phasec).*envc;
    steps = RELAX_detect_steps(mkEEG(x));
    assert(isempty(steps.sample), 'a 2 Hz slow wave is not a pop, %d step(s) found', numel(steps.sample));

    % a slow random walk likewise (the level gate found 2)
    walkc = cumsum(randn(1, n60c)); walkc = walkc - movmean(walkc, 4001);
    walkc = 150 * walkc / max(abs(walkc));
    x = quietc + walkc.*envc;
    steps = RELAX_detect_steps(mkEEG(x));
    assert(isempty(steps.sample), 'a slow random walk is not a pop, %d step(s) found', numel(steps.sample));

    % a real pop is still found, with its amplitude and its position
    x = quietc + 200 * (t60c >= 30000);
    steps = RELAX_detect_steps(mkEEG(x));
    assert(numel(steps.sample) == 1, 'one step expected, got %d', numel(steps.sample));
    assert(abs(steps.amplitude - 200) < 20, 'amplitude %g should be near 200', steps.amplitude);
    assert(abs(steps.sample - 30000) < 100, 'step located at %d, expected near 30000', steps.sample);

    % a 50 ms sigmoid is a pop
    x = quietc + 100 ./ (1 + exp(-(t60c - 30000) / 10));
    steps = RELAX_detect_steps(mkEEG(x));
    assert(numel(steps.sample) == 1, 'a 50 ms sigmoid is a pop, got %d step(s)', numel(steps.sample));

    % two steps 300 ms apart do not mask each other through the shared baseline span
    x = quietc + 150 * (t60c >= 30000) + 150 * (t60c >= 30300);
    steps = RELAX_detect_steps(mkEEG(x));
    assert(numel(steps.sample) == 2, 'two steps 300 ms apart, got %d', numel(steps.sample));

    % the red line: a genuine step inside a drifting stretch is found, and found alone - the
    % level gate reported it with two spurious neighbours. The tolerated drift is bounded: at
    % the provisional ZThreshold of 4 the step's slope must clear the drift's own by roughly a
    % factor of two, which is what the 2026-09-22 threshold sweep exists to price.
    x = quietc + 75*sin(phasec).*envc + 200 * (t60c >= 30000);
    steps = RELAX_detect_steps(mkEEG(x));
    assert(numel(steps.sample) == 1, 'a step on a drifting baseline must be found alone, got %d', ...
        numel(steps.sample));
    assert(abs(steps.sample - 30000) < 150, 'step on drift located at %d', steps.sample);

    % new options are echoed, and the channel scale is the slope statistic's
    [~, info] = RELAX_detect_steps(mkEEG(quietc), 'SlopeBaselineMs', 1500);
    assert(info.params.SlopeBaselineMs == 1500 && info.params.SlopeNoiseFloorMicrovolts == 0.1);
    assert(mod(info.params.derived.baseline, 2) == 1, 'the baseline span must be an odd sample count');
    assert(~isfield(info, 'channelScaleD') && isfield(info, 'channelScaleSlope'));

    %% ---- 1d. derived sample counts are odd and echoed ----
    [~, info] = RELAX_detect_steps(mkEEG(randn(1, n) * 5), 'WindowMs', 150, 'SmoothMs', 10);
    assert(info.params.WindowMs == 150);
    assert(mod(info.params.derived.W, 2) == 1 && info.params.derived.W >= 150, ...
        'W must be forced odd (150 samples -> 151)');
    assert(mod(info.params.derived.smooth, 2) == 1 && info.params.derived.smooth >= 10);
    [~, info] = RELAX_detect_steps(mkEEG(randn(1, n) * 5), 'ZThreshold', 33);
    assert(info.params.ZThreshold == 33, 'parameters must be echoed in info.params');

    %% ---- 1e. no usable segment -> unavailable with a status, not an error ----
    EEG = mkEEG(randn(1, 1000) * 5, 200.5:200:900);   % five 200-sample segments
    [steps, info] = RELAX_detect_steps(EEG);
    assert(~info.available && ~isempty(info.status), 'short segments must yield available=false');
    assert(isequal(size(steps.sample), [0 1]));

    %% ---- 1f. TransitionMs is Dslope's exact lag, down to an uncentred one-sample difference ----
    rng(7);
    q = randn(1, 20000) * 5;
    for tr = [1 2 4 6 10 100]
        [~, info] = RELAX_detect_steps(mkEEG(q), 'TransitionMs', tr);
        d = info.params.derived;
        assert(d.lag == tr, 'TransitionMs %d at 1 kHz must be a %d-sample lag, got %d', tr, tr, d.lag);
        assert(d.gPre == floor(tr/2) && d.gPost == tr - floor(tr/2), ...
            'lag %d must split floor/ceil, got %d + %d', tr, d.gPre, d.gPost);
    end
    assert(~isfield(info.params.derived, 'g'), 'the symmetric half-gap g was replaced by gPre/gPost');

    % an ideal step is found once at every lag in the grid, at its edge, at full amplitude, and
    % its interval contains the edge and stays local
    x = q + 200 * ((0:19999) >= 10000);
    for tr = [1 2 4 6 10]
        steps = RELAX_detect_steps(mkEEG(x), 'TransitionMs', tr);
        assert(numel(steps.sample) == 1, 'lag %d: one step expected, got %d', tr, numel(steps.sample));
        assert(abs(steps.sample - 10000) <= tr, 'lag %d: step at %d, expected within %d of 10000', ...
            tr, steps.sample, tr);
        assert(abs(steps.amplitude - 200) < 10, 'lag %d: amplitude %g, expected ~200', tr, steps.amplitude);
        assert(steps.onset <= 10000 && steps.offset >= 9999, ...
            'lag %d: interval [%d %d] must contain the edge', tr, steps.onset, steps.offset);
        assert(steps.offset - steps.onset < 50, ...
            'lag %d: the gate run must stay local, got %d samples', tr, steps.offset - steps.onset);
    end

    % the lag-1 windows are indexed from the segment start: a step 600 samples in (two windows)
    % is found at its edge, not shifted and not an index error. Z 20, not the provisional 12: at
    % W 301 / lag 1 a median-bend shoulder 12 samples past this edge reaches z 14 for this noise
    % draw (1 of 40 seeds; none at Z 20); how often shoulders occur is measured by the survey.
    x = q; x(1:600) = x(1:600) + 100;              % a step 600 samples in: inside the valid span
    steps = RELAX_detect_steps(mkEEG(x), 'TransitionMs', 1, 'ZThreshold', 20);
    assert(numel(steps.sample) == 1 && abs(steps.sample - 600) <= 1, ...
        'lag 1: a step 600 samples in must be found at its edge');

    %% ---- 1g. diagnostics for the noise floor and the channel's raw noise ----
    [~, iq]  = RELAX_detect_steps(mkEEG(q), 'TransitionMs', 1);
    [~, iq5] = RELAX_detect_steps(mkEEG(q * 100), 'TransitionMs', 1);
    for f = {'floorBoundFraction', 'rawScale', 'diffScale'}
        assert(isequal(size(iq.(f{1})), [1 1]), 'info.%s must be nbchan x 1', f{1});
    end
    assert(iq.floorBoundFraction >= 0 && iq.floorBoundFraction <= 1);
    assert(iq5.floorBoundFraction < iq.floorBoundFraction, ...
        'a fixed floor binds less on a noisier channel: %g vs %g', iq5.floorBoundFraction, iq.floorBoundFraction);
    assert(abs(iq.rawScale - 5) < 0.5, 'rawScale %g of white noise sd 5', iq.rawScale);
    assert(abs(iq.diffScale - 5 * sqrt(2)) < 0.7, 'diffScale %g of white noise sd 5', iq.diffScale);

    %% ---- 1h. the defaults are the detection-only operating point of the entry script ----
    [~, info] = RELAX_detect_steps(mkEEG(q));
    P = info.params;
    assert(isequal(P.RefineTransition, true), ...
        'the shipped refinement is a single bracketing pass (user decision, 2026-09-25)');
    threw = false;
    try
        RELAX_detect_steps(mkEEG(q), 'RefineTransition', 'once');
    catch
        threw = true;
    end
    assert(threw, 'RefineTransition is a logical: the retired ''once'' spelling must error');
    assert(P.WindowMs == 151 && P.TransitionMs == 4 && P.ZThreshold == 20 && ...
           P.ChannelZThreshold == 35 && P.MinAmplitudeMicrovolts == 15 && P.RefractoryMs == 10 && ...
           P.CommonModeChannelFraction == 0.15 && ...
           P.SlopeBaselineMs == 1000 && P.SlopeNoiseFloorMicrovolts == 0.1 && ...
           P.SlopeSpreadQuantile == 0.85 && P.CommonModeGapMs == 10 && ...
           ~isfield(P, 'CommonModeWindowMs'), ...
        'the shipped defaults are the ones the 2026-09-23 survey chose');
    for f = {'MaxTransitionMs', 'RefineIterations', 'MaxPlateauDriftFraction', 'HoldMs', 'MaxHoldDriftFraction'}
        assert(~isfield(P, f{1}), 'the %s screen was removed (Task 4D, 2026-09-25)', f{1});
    end

    % At the defaults an exponential approach (tau 20 ms, 200 uV) is bracketed, not reported as
    % the gate's raw run: the step sits at the mid-level crossing (20 ln 2 = 14 ms in) rather than
    % at the steepest sample, and the interval spans the rise to within 10 % of the new level.
    % (the draw is taken off the file's stream, so the seed-dependent drift section in 3 sees the
    % same noise it was written against)
    streamBefore = rng;
    tt = 0:n-1;
    x = randn(1, n) * 5; x(4001:end) = x(4001:end) + 200 * (1 - exp(-(tt(4001:end) - 4000) / 20));
    steps = RELAX_detect_steps(mkEEG(x));
    rng(streamBefore);
    assert(numel(steps.sample) == 1, 'default exp approach: expected 1 step, got %d', numel(steps.sample));
    assert(steps.sample >= 4011 && steps.sample <= 4020, ...
        'the default step sits at the mid-level crossing, got %d', steps.sample);
    assert(steps.onset >= 3990 && steps.onset <= 4008 && steps.offset >= 4030 && steps.offset <= 4070, ...
        'the default interval brackets the rise, got [%d %d]', steps.onset, steps.offset);
    assert(abs(steps.amplitude - 200) < 15, 'default exp amplitude %g', steps.amplitude);

    %% ---- 1i. the spread is a high quantile, so it survives a mostly-zero short-lag Dslope ----
    % Neighbouring centred medians are equal whenever the entering and leaving samples fall on the
    % same side of the median, so at lag 1 most of Dslope is exactly 0 and its MAD collapses.
    [~, iq9]   = RELAX_detect_steps(mkEEG(q), 'TransitionMs', 1);
    [~, iqMad] = RELAX_detect_steps(mkEEG(q), 'TransitionMs', 1, 'SlopeSpreadQuantile', 0.5);
    assert(iq9.channelScaleSlope > 0.01, ...
        'the default-q channel scale must not collapse at lag 1, got %g', iq9.channelScaleSlope);
    assert(iqMad.channelScaleSlope < iq9.channelScaleSlope / 10, ...
        'q 0.5 is the MAD, which does collapse at lag 1: %g vs %g', ...
        iqMad.channelScaleSlope, iq9.channelScaleSlope);

    % a DC offset between inter-boundary segments is not channel noise
    x = q; x(10001:end) = x(10001:end) + 1000;
    [~, iseg] = RELAX_detect_steps(mkEEG(x, 10000.5), 'TransitionMs', 1);
    assert(abs(iseg.rawScale - 5) < 0.5, ...
        'rawScale %g must measure noise within segments, not the 1000 uV offset between them', ...
        iseg.rawScale);

    %% ---- 2a. a clean step is found once, located, signed and scored ----
    x = randn(1, n) * 5;
    x(4001:end) = x(4001:end) + 500;
    EEG = mkEEG([x; randn(1, n) * 5]);
    [steps, info] = RELAX_detect_steps(EEG);
    assert(info.available, '%s', info.status);
    k = find(steps.channel == 1);
    assert(numel(k) == 1, 'expected exactly 1 step on channel 1, got %d', numel(k));
    assert(abs(steps.sample(k) - 4001) <= 3, 'step sample should be ~4001, got %d', steps.sample(k));
    assert(abs(steps.amplitude(k) - 500) < 25, 'amplitude should be ~+500, got %g', steps.amplitude(k));
    assert(steps.z(k) > 20, 'a 100-sigma step must score far above threshold');
    assert(steps.channelZ(k) > 4, 'channel-wide z must also clear its threshold');
    assert(isempty(steps.sample(steps.channel == 2)), 'clean channel 2 must yield no steps');
    assert(strcmp(steps.label{k}, 'E1'));
    assert(info.nSteps == 1 && info.nChannelsAffected == 1);

    x2 = randn(1, n) * 5;
    x2(4001:end) = x2(4001:end) - 500;
    steps = RELAX_detect_steps(mkEEG(x2));
    assert(numel(steps.sample) == 1 && steps.amplitude(1) < 0, 'downward step keeps its sign');

    %% ---- 2b. two same-sign edges closer than the window are one candidate ----
    x = randn(1, n) * 5;
    x(4001:end) = x(4001:end) + 400;
    x(4006:end) = x(4006:end) + 400;
    steps = RELAX_detect_steps(mkEEG(x));
    assert(numel(steps.sample) == 1, 'two edges 5 ms apart are one step, got %d', numel(steps.sample));
    assert(abs(steps.amplitude(1) - 800) < 40, 'merged step carries the full amplitude, got %g', steps.amplitude(1));

    x = randn(1, n) * 5;
    x(3001:end) = x(3001:end) + 400;
    x(3501:end) = x(3501:end) - 400;
    steps = RELAX_detect_steps(mkEEG(x));
    assert(numel(steps.sample) == 2, 'opposite edges 500 ms apart are two steps, got %d', numel(steps.sample));
    assert(steps.amplitude(1) > 0 && steps.amplitude(2) < 0);

    %% ---- 2c. the absolute floor stops a tiny step on a silent channel ----
    % On a silent channel z is capped near amplitude / SlopeNoiseFloorMicrovolts (about 20 for
    % this 2 uV step), so the local gate is pinned at 12 to leave the floor as the only test.
    x = randn(1, n) * 0.01;
    x(4001:end) = x(4001:end) + 2;
    steps = RELAX_detect_steps(mkEEG(x), 'ZThreshold', 12, 'MinAmplitudeMicrovolts', 10);
    assert(isempty(steps.sample), 'a 2 uV step must be below a 10 uV floor');
    steps = RELAX_detect_steps(mkEEG(x), 'ZThreshold', 12, 'MinAmplitudeMicrovolts', 1);
    assert(numel(steps.sample) == 1, 'the same step is found once the floor is lowered');

    %% ---- 2d. one transition is reported once, however many gate runs it produced ----
    % An exponential approach's edge and its median-filter shoulder can be separate gate runs
    % that bracket to the same transition (final review, 2026-09-25: 12 overlapping pairs in
    % 100 such pops; 99-572 same-sample pairs per file on the four real files).
    nOverlap = 0;
    for seed = 1:100
        rs = RandStream('twister', 'Seed', seed);
        t = 0:7999;
        x = randn(rs, 1, 8000) * (3 + 10*rand(rs));
        tau = 5 + 60*rand(rs); A = (30 + 300*rand(rs)) * sign(randn(rs));
        x(4001:end) = x(4001:end) + A * (1 - exp(-(t(4001:end) - 4000) / tau));
        steps = RELAX_detect_steps(mkEEG(x));
        [on, ord] = sort(steps.onset); off = steps.offset(ord);
        nOverlap = nOverlap + sum(on(2:end) <= cummax(off(1:end-1)));
    end
    assert(nOverlap == 0, '%d same-channel steps overlap an earlier one', nOverlap);

    % Section 3 tests the transition bracketing at the legacy operating point, whose long
    % TransitionMs the slow shapes need. The plateau- and hold-drift gates it once also tested
    % were removed (Task 4D, 2026-09-25).

    %% ---- 3. transition bracketing on the shapes seen in the survey (10 uV noise) ----
    function [steps, k] = one(x, varargin)
        steps = RELAX_detect_steps(mkEEG(x), legacy{:}, varargin{:});
        k = find(steps.channel == 1);
    end
    t = 0:n-1;

    % clean step: bracket is [4000, 4001], zero rise
    x = randn(1, n) * 10; x(4001:end) = x(4001:end) + 500;
    [steps, k] = one(x);
    assert(numel(k) == 1);
    assert(abs(steps.onset(k) - 4000) <= 1 && abs(steps.offset(k) - 4001) <= 1, ...
        'clean step bracket should be [4000 4001], got [%d %d]', steps.onset(k), steps.offset(k));
    assert(steps.riseMs(k) <= 1, 'clean step reports ~0 ms rise, got %g', steps.riseMs(k));
    assert(abs(steps.amplitude(k) - 500) < 15);

    % exponential approach, tau 40 ms, 200 uV
    x = randn(1, n) * 10; x(4001:end) = x(4001:end) + 200 * (1 - exp(-(t(4001:end) - 4000) / 40));
    [steps, k] = one(x);
    assert(numel(k) == 1, 'exp approach: expected 1 step, got %d', numel(k));
    assert(steps.onset(k) >= 3998 && steps.onset(k) <= 4008, 'exp onset %d', steps.onset(k));
    assert(steps.offset(k) >= 4065 && steps.offset(k) <= 4100, 'exp offset %d', steps.offset(k));
    assert(abs(steps.amplitude(k) - 200) < 15, 'exp amplitude %g', steps.amplitude(k));
    assert(steps.riseMs(k) >= 55 && steps.riseMs(k) <= 100);

    % staircase: two 100 uV rungs 80 ms apart is ONE step spanning both
    x = randn(1, n) * 10; x(4001:end) = x(4001:end) + 100; x(4081:end) = x(4081:end) + 100;
    [steps, k] = one(x);
    assert(numel(k) == 1, 'staircase: expected 1 step, got %d', numel(k));
    assert(steps.onset(k) <= 4001 && steps.offset(k) >= 4080, 'staircase bracket [%d %d]', steps.onset(k), steps.offset(k));
    assert(abs(steps.amplitude(k) - 200) < 15, 'staircase amplitude %g', steps.amplitude(k));

    % 200 ms linear ramp of 200 uV. A single bracketing pass settles against levels measured at
    % the gate's strongest sample, so on a ramp it reports the onset up to ~45 ms late (4018-4045
    % over 20 seeds, 2026-09-25); the iterative bracketing this bound was written for is removed.
    x = randn(1, n) * 10; ramp = min(1, max(0, (t - 4000) / 200)); x = x + 200 * ramp;
    [steps, k] = one(x);
    assert(numel(k) == 1, 'ramp: expected 1 step, got %d', numel(k));
    assert(steps.onset(k) >= 4001 && steps.onset(k) <= 4050, 'ramp onset %d', steps.onset(k));
    assert(steps.offset(k) >= 4165 && steps.offset(k) <= 4200, 'ramp offset %d', steps.offset(k));
    assert(abs(steps.amplitude(k) - 200) < 15, 'ramp amplitude %g', steps.amplitude(k));

    % 40 uV / 300 ms pre-ramp followed by a +200 jump
    x = randn(1, n) * 10; pre = 40 * min(1, max(0, (t - 3700) / 300)); x = x + pre; x(4001:end) = x(4001:end) + 200;
    [steps, k] = one(x);
    assert(numel(k) == 1, 'pre-ramp: expected 1 step, got %d', numel(k));
    assert(abs(steps.offset(k) - 4001) <= 2, 'pre-ramp offset %d', steps.offset(k));   % single pass: 4003 on 1 of 20 seeds
    assert(steps.amplitude(k) >= 195 && steps.amplitude(k) <= 225, 'pre-ramp amplitude %g', steps.amplitude(k));

    % 150 ms / 200 uV pulse -> two opposite steps; 50 ms pulse -> nothing
    x = randn(1, n) * 10; x(4001:4150) = x(4001:4150) + 200;
    [steps, k] = one(x);
    assert(numel(k) == 2, '150 ms pulse: expected 2 steps, got %d', numel(k));
    assert(abs(steps.sample(k(1)) - 4001) <= 2 && abs(steps.sample(k(2)) - 4151) <= 2);
    assert(abs(steps.amplitude(k(1)) - 200) < 20 && abs(steps.amplitude(k(2)) + 200) < 20);
    x = randn(1, n) * 10; x(4001:4050) = x(4001:4050) + 200;
    [~, k] = one(x);
    assert(isempty(k), '50 ms pulse must not be a step (got %d)', numel(k));

    % continuous oscillation is not a step
    x = randn(1, n) * 10 + 200 * sin(2*pi*1.0*t/srate);
    [~, k] = one(x); assert(isempty(k), '1 Hz 200 uV sine produced %d steps', numel(k));
    x = randn(1, n) * 10 + 300 * sin(2*pi*0.5*t/srate);
    [~, k] = one(x); assert(isempty(k), '0.5 Hz 300 uV sine produced %d steps', numel(k));

    % random-walk drift (1/f^2, 150 uV rms, 60 s) with five pops of mixed shape (a 150 ms ramp
    % was removed 2026-09-24: transitions that slow are not what this detector is for). At 200 uV
    % rms the channel gate sits at z~4 for a 200 uV pop and recovery becomes seed-dependent.
    n60 = 60 * srate; f = (0:n60-1)'; f(1) = 1; f = min(f, n60 - f + 1);
    w = real(ifft(fft(randn(n60,1)) ./ f)); w = 150 * w / std(w); w = w';
    t60 = 0:n60-1;
    at = [8 18 28 38 55] * srate + 1;
    w(at(1):end) = w(at(1):end) + 300;
    w(at(2):end) = w(at(2):end) - 250 * (1 - exp(-(t60(at(2):end) - at(2)) / 40));
    w(at(3):end) = w(at(3):end) + 120; w(at(3)+60:end) = w(at(3)+60:end) + 120;
    w(at(4):end) = w(at(4):end) + 300;
    w(at(5):end) = w(at(5):end) - 250 * (1 - exp(-(t60(at(5):end) - at(5)) / 40));
    % Z 20, not the provisional 12: at q 0.85 a large edge can carry a median-bend shoulder at z 12-14,
    % and 1/f drift gives z 13-18 false steps (job 19909316); this section is not testing Z.
    steps = RELAX_detect_steps(mkEEG(w + randn(1, n60) * 10), 'ZThreshold', 20);
    hit = arrayfun(@(a) any(abs(steps.sample - a) <= 200), at);
    false_ = sum(~any(abs(steps.sample - at) <= 200, 2));
    assert(sum(hit) >= 4, 'drift channel: recovered only %d/5 pops', sum(hit));
    assert(false_ <= 2, 'drift channel: %d false steps', false_);

    % a NaN block elsewhere on the channel does not break detection
    x = randn(1, n) * 10; x(8001:end) = x(8001:end) + 300; x(5000:5010) = NaN;
    [steps, k] = one(x);
    assert(numel(k) == 1 && abs(steps.sample(k) - 8001) <= 3, 'NaN block: expected the 8001 step');

    %% ---- 4. boundaries ----
    x = randn(1, n) * 5; x(5001:end) = x(5001:end) + 800;
    steps = RELAX_detect_steps(mkEEG(x, 5000.5));
    assert(isempty(steps.sample), 'a jump exactly at a boundary is a splice, not a step');
    x(7001:end) = x(7001:end) + 400;
    steps = RELAX_detect_steps(mkEEG(x, 5000.5));
    assert(numel(steps.sample) == 1 && abs(steps.sample(1) - 7001) <= 3, 'step in segment 2 must be found');
    % The median-slope gate reads the pre-level at mc(t-g), centred rather than trailing, so a
    % step must sit about W samples into the segment to be seen, where the superseded level gate
    % reached back about W/2 (~140 samples at W 150, measured 2026-09-22, job 19867003). At the
    % W 151 default, 350 samples in leaves a margin.
    x = randn(1, n) * 5; x(5350:end) = x(5350:end) + 400;           % 350 samples into segment 2
    steps = RELAX_detect_steps(mkEEG(x, 5000.5));
    assert(numel(steps.sample) == 1 && abs(steps.sample(1) - 5350) <= 3, ...
        'a step well inside a segment is still bracketed from the first valid centre');
    x = randn(1, n) * 5; x(5031:end) = x(5031:end) + 400;           % 30 samples in: no pre-level to measure
    steps = RELAX_detect_steps(mkEEG(x, 5000.5));
    assert(isempty(steps.sample), 'with no pre-level inside the segment the step cannot be reported');
    x = randn(1, n) * 5; x(7001:end) = x(7001:end) + 400;
    [steps, info] = RELAX_detect_steps(mkEEG(x, [100.5 300.5]));   % two short segments + one long
    assert(info.available && numel(steps.sample) == 1, 'short segments are skipped, the long one screened');

    %% ---- 5. co-occurrence events and the common-mode flag ----
    nc = 40; X = randn(nc, n) * 10;
    for c = 1:30                                                   % rim-wide, mixed sign, jittered
        at = 5000 + round(40 * (rand - 0.5));
        X(c, at:end) = X(c, at:end) + 300 * sign(randn);
    end
    X(35, 2000:end) = X(35, 2000:end) + 300;                        % private pop
    for c = 36:38; X(c, 8000:end) = X(c, 8000:end) + 250; end       % 3-channel same-sign cluster
    % Z 20, not the provisional 12: at q 0.85 a large edge can carry a median-bend shoulder at z 12-14,
    % and 1/f drift gives z 13-18 false steps (job 19909316); this section is not testing Z.
    [steps, info] = RELAX_detect_steps(mkEEG(X), 'ZThreshold', 20);
    assert(numel(steps.sample) == 34, 'expected 34 steps, got %d', numel(steps.sample));
    cm = steps.commonMode;
    assert(sum(cm) == 30 && all(steps.channel(cm) <= 30), 'the 30 rim-wide steps are common mode');
    assert(~any(cm(steps.channel >= 35)), 'private / cluster steps are not common mode');
    assert(info.nCommonModeSteps == 30 && info.nCommonModeEvents == 1);
    clusterIds = unique(steps.eventId(steps.channel >= 36));
    assert(numel(clusterIds) == 1, 'the cluster shares one eventId');
    ev = info.events;
    assert(ev.nChannels(clusterIds) == 3 && ~ev.commonMode(clusterIds) && ev.sameSignFraction(clusterIds) == 1);
    assert(ev.lo(clusterIds) <= 8000 && ev.hi(clusterIds) >= 8000);
    assert(isempty(info.badChannels), 'a single pop per channel is not a bad contact');
    assert(all(info.stepRatePerMinute(1:30) == 0), 'common-mode steps do not count against the channel');

    %% ---- 5b. events chain steps no more than CommonModeGapMs apart, with no fixed bins ----
    % At lag 4 a step is placed within [edge-2, edge+1], so two placements differ from their
    % edges' separation by at most 3 samples; every separation below keeps that margin from 10.
    X = randn(12, n) * 5;
    edges = [3000 3005, 5000 5016, 7000 7006 7012 7018 7024, 8960 9020];
    for c = 1:numel(edges); X(c, edges(c):end) = X(c, edges(c):end) + 300; end
    % Z 20, not the provisional 12: at q 0.85 a large edge can carry a median-bend shoulder at z 12-14,
    % and 1/f drift gives z 13-18 false steps (job 19909316); this section is not testing Z.
    [steps, info] = RELAX_detect_steps(mkEEG(X), 'ZThreshold', 20);
    assert(isequal(sort(steps.channel)', 1:numel(edges)), 'fixture: one step on each of channels 1-11');
    ev = @(c) steps.eventId(steps.channel == c);
    assert(ev(1) == ev(2), 'steps 5 ms apart are one event');
    assert(ev(3) ~= ev(4), 'steps 16 ms apart are two events');
    assert(all(arrayfun(ev, 5:9) == ev(5)), 'a chain of 6 ms gaps spanning 24 ms is one event');
    assert(ev(10) ~= ev(11), ['steps 60 ms apart are two events, though 50 ms bins ' ...
        '(8951-9000, 9001-9050) would have joined them']);
    assert(numel(info.events.id) == 6, 'expected 6 events, got %d', numel(info.events.id));
    steps = RELAX_detect_steps(mkEEG(X), 'ZThreshold', 20, 'CommonModeGapMs', 20);
    ev = @(c) steps.eventId(steps.channel == c);
    assert(ev(3) == ev(4) && ev(10) ~= ev(11), 'CommonModeGapMs sets the gap');

    %% ---- 6. bad contacts, worst-first, with the binding criterion named ----
    n60 = 60 * srate;
    x = randn(1, n60) * 5; lvl = 0;
    for a = 2001:2000:n60; lvl = 300 - lvl; x(a:end) = x(a:end) + (lvl - (300 - lvl)); end  % telegraph every 2 s
    [steps, info] = RELAX_detect_steps(mkEEG([x; randn(1, n60) * 5]));
    assert(numel(steps.sample) >= 28, 'telegraph: expected ~29 toggles, got %d', numel(steps.sample));
    assert(isequal(info.badChannels, {'E1'}), 'the telegraph channel is a bad contact');
    assert(contains(info.badChannelReason{1}, 'steps/min'));

    % ordering: 108 / 72 / 36 steps per minute on 60 s (rate default 6/min) -> worst first
    rates = [108 72 36]; X = randn(3, n60) * 5;
    for c = 1:3
        at = round(linspace(2000, n60 - 2000, rates(c)));
        for j = 1:numel(at); X(c, at(j):end) = X(c, at(j):end) + 300 * (-1)^j; end
    end
    [~, info] = RELAX_detect_steps(mkEEG(X([2 1 3], :)));   % E1=72/min, E2=108/min, E3=36/min
    assert(isequal(info.badChannels, {'E2'; 'E1'; 'E3'}), 'bad channels must be worst-first');
    assert(isequal(size(info.badChannelReason), [3 1]));

    % transition fraction: 90 alternating 15 ms ramps in a minute (transitions no longer than
    % pops' ~15 ms; spaced ~660 ms, so each has a full W-sample window on both sides). A rate cap
    % of 1000/min cannot bind, so the fraction criterion must be the one that names itself.
    x = randn(1, n60) * 5; tt = 0:n60-1;
    at = round(linspace(2000, n60 - 2000, 90));
    for j = 1:numel(at)
        x = x + 300 * min(1, max(0, (tt - at(j)) / 15)) * (-1)^j;
    end
    % Z 20 as in sections 3 and 5: at Z 12 a few 300 uV edges carry a shoulder (93 steps, job 19909403)
    [steps, info] = RELAX_detect_steps(mkEEG(x), 'ZThreshold', 20, 'StepsPerMinuteToRejectChannel', 1000);
    assert(abs(numel(steps.sample) - 90) <= 2, 'fixture: expected ~90 steps, got %d', numel(steps.sample));
    assert(info.transitionFraction(1) > 0.02, 'transition fraction should exceed 2%%, got %g', info.transitionFraction(1));
    assert(isequal(info.badChannels, {'E1'}) && contains(info.badChannelReason{1}, 'transition'));

    % telegraph every 250 ms: every transition is abrupt and large, so the channel rejects on
    % step rate; the sample-level instability criterion it used to need was retired 2026-09-22
    x = randn(1, n60) * 5; sq = 300 * (mod(floor(tt / 250), 2) == 0); x = x + sq;
    [~, info] = RELAX_detect_steps(mkEEG(x));
    assert(isequal(info.badChannels, {'E1'}), 'telegraph channel must still be rejected');
    assert(contains(info.badChannelReason{1}, 'steps/min'), ...
        'telegraph channel must now reject by rate, got: %s', info.badChannelReason{1});
    assert(isequal(size(info.badChannels), size(info.badChannelReason)));

    fprintf('test_RELAX_detect_steps: PASS\n');
end
