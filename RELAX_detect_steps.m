% RELAX_detect_steps() - Detect electrode-pop level shifts, including slow-rising ones, in
%                        continuous EEG, and estimate each shift's transition interval.
%
% Usage:
%   >> [steps, info] = RELAX_detect_steps(EEG);
%   >> [steps, info] = RELAX_detect_steps(EEG, 'ZThreshold', 5, 'TransitionMs', 150);
%
% On a saline-sponge net (EGI HydroCel GSN) an electrode pop is a change in a channel's DC
% level that takes anywhere from one sample to ~100 ms to complete - staircases, sigmoids
% and exponential approaches all occur - and then holds. A derivative detector spreads such
% a transition over many samples and never fires, and a plain before/after level comparison
% cannot tell a pop from sustained drift, which reaches the same level change given enough
% time. This detector instead scores how FAST the level moved: a centred median filter (the
% mc(t) = mPre(t+h) construction below, half-width h = floor(W/2)) is differenced across the
% TransitionMs lag alone (Dslope(t) = mc(t+gPost) - mc(t-gPre), gPre = floor(lag/2),
% gPost = lag - gPre), and that slope is compared to
% a SlopeBaselineMs-long local TREND in the channel's own recent Dslope (a signed median b(t),
% plus a robust spread s(t) around it): monotonic drift sets its own Dslope trend and is
% absorbed by b, oscillating drift is absorbed by the widened s, while a pop's slope stands out
% as a departure from both. Dslope's sign is the candidate's sign. The signed, longer-lag
% level difference D(t) (post-window median minus pre-window median, W samples either side of
% a TransitionMs gap) is checked against the absolute floor.
%
% Two gates decide a candidate:
%   local   |Dslope - local trend b| / max(local robust spread s, SlopeNoiseFloorMicrovolts) >
%           ZThreshold, and |D| >= MinAmplitudeMicrovolts;
%   channel |Dslope| / (SlopeSpreadQuantile-quantile spread of Dslope over this channel's
%           locally-quiet samples, scaled to Gaussian sd units) > ChannelZThreshold. Both the
%           local spread s(t) and the channel spread are the q-quantile of |deviation| divided
%           by the Gaussian quantile at q (q = 0.5 reduces to the familiar 1.4826 * MAD); a high
%           q (default 0.85) is used because at short TransitionMs lags two neighbouring centred
%           medians are equal whenever the entering and leaving samples fall on the same side of
%           the median, so most of Dslope is exactly 0 and its MAD collapses to (near) zero.
%           A random-walk drift (sweat, an electrode drying out) makes a
%           plain level-difference test fire continuously because nothing inside a window is
%           stationary; scoring the slope against its own local and channel-wide baseline
%           asks whether THIS level change happened unusually FAST, not merely how far it went.
% Same-sign runs of candidate samples are one step each. Around every run the transition
% is then bracketed once: the smoothed trace (movmedian SmoothMs) is walked outward from its
% mid-level crossing until it settles within a band of the pre-level (onset) and of the
% post-level (offset), and the levels are re-measured from windows adjacent to those samples.
% Intermediate staircase rungs lie outside both bands, so they stay inside the interval. An
% edge run and a separate gate run more than RefractoryMs away (a median-filter "shoulder" of
% the same edge, say) can both bracket to the same or an overlapping transition; when
% RefineTransition is true those same-sign, same-channel brackets are merged into one step
% with RELAX_merge_overlapping_steps (opposite-sign overlaps are left alone - they are a brief
% excursion, not a duplicate).
%
% Run this BEFORE any filtering, at the native rate: a zero-phase highpass converts a step
% into a bipolar excursion spread over seconds in BOTH time directions.
%
% Inputs:
%   EEG - EEGLAB EEG structure (uses .data, .nbchan, .pnts, .srate, .chanlocs, .event)
%   Name-value pairs (all optional):
%     WindowMs             - W, level window either side (default 151; forced to an odd
%                             sample count: an even movmedian window averages its two
%                             middle samples and splits an edge into two half-steps)
%     TransitionMs         - the exact lag of Dslope in samples (min 1), split floor/ceil
%                             either side of t (default 4). 1 is the uncentred difference
%                             mc(t+1) - mc(t) with adjacent level windows.
%     ZThreshold           - local z threshold (default 20, the 2026-09-23 slope-gate survey;
%                             see docs/superpowers/reports/2026-09-23-slope-gate-survey.md, section 6)
%     ChannelZThreshold    - channel-wide z threshold (default 35, same survey)
%     MinAmplitudeMicrovolts - absolute floor on the level change (default 15)
%     RefractoryMs         - same-sign candidate runs separated by less than this are one
%                             step (default 10). Gap closing only; it is NOT a minimum
%                             inter-step spacing (runs already span ~W+G).
%     SmoothMs             - movmedian length for the transition bracket (default 11, odd)
%     BoundaryFraction     - settle band as a fraction of |amplitude| (default 0.10)
%     BoundaryNoiseMultiple - settle band floor, in robust sd of the smoothed window
%                             (default 2)
%     RefineTransition     - true (default) brackets the transition as described above: the
%                             run's strongest sample is the starting midpoint (the run is a
%                             plateau of the gate statistic, so its geometric centre does not
%                             locate the edge), and a walk that never settles falls back to the
%                             run's own end rather than rejecting. The only rejections left are
%                             a re-measured level change under MinAmplitudeMicrovolts and a
%                             bracket too close to a segment edge to re-measure (info.nRejected;
%                             neither fired on the four surveyed files). false skips bracketing entirely: onset and
%                             offset become the candidate run's own ends, the step is placed at
%                             the run's strongest sample, and amplitude is the signed D there.
%                             The interval is then about TransitionMs wide and centred on the
%                             edge rather than an estimate of the true transition. Diagnostic:
%                             it shows the gate's raw output with no interpretation.
%     StepsPerMinuteToRejectChannel - per-channel step rate (excluding common-mode steps)
%                             above which the channel is a bad contact (default 6)
%     MaxTransitionFractionToRejectChannel - fraction of the recording a channel spends in
%                             transition above which it is a bad contact (default 0.02)
%     CommonModeGapMs      - steps whose .sample are within this many ms of the previous
%                             step (across all channels) chain into one event; a larger gap
%                             starts a new one (default 10)
%     CommonModeChannelFraction - an event with more than this fraction of channels (and
%                             at least two) is common-mode - a whole-head shift (movement / reference), not
%                             a contact pop (default 0.15, ~38 channels of 256: above the
%                             local clusters of 2-12 seen in the survey, below the rim-wide
%                             events of 30+).
%     SlopeBaselineMs      - length of the local window used to estimate Dslope's own trend
%                             (signed median) and spread (robust sd), which a candidate is
%                             scored as a departure from (default 1000; forced odd)
%     SlopeNoiseFloorMicrovolts - floor on the local trend's robust spread, so z does not blow
%                             up where the channel is dead quiet (default 0.1)
%     SlopeSpreadQuantile  - the local spread s(t) and the channel scale are the q-quantile of
%                             |deviation| divided by the Gaussian quantile at q, so they stay in
%                             sd units; q = 0.5 is exactly 1.4826 * MAD (default 0.85). A high q
%                             is needed because at short TransitionMs lags most of Dslope is
%                             exactly 0 (two neighbouring centred medians agree whenever the
%                             entering and leaving samples fall on the same side of the median),
%                             which collapses the MAD (q = 0.5) to (near) zero
%     ThresholdPairs       - N x 2 matrix of [ZThreshold ChannelZThreshold] rows (default [],
%                             meaning use the ZThreshold/ChannelZThreshold options above).
%                             When given, the per-channel statistics that do not depend on the
%                             thresholds (levelStatistic, floorBoundFraction, rawScale,
%                             diffScale) are computed once per channel and reused for every row,
%                             steps becomes an N x 1 struct array and info an N x 1 struct
%                             array, and element k is identical to a separate call with
%                             'ZThreshold', ThresholdPairs(k,1), 'ChannelZThreshold',
%                             ThresholdPairs(k,2) (info(k).params carries that pair's own
%                             ZThreshold/ChannelZThreshold, alongside ThresholdPairs itself).
%                             Malformed values (not N x 2 with N >= 1, or non-finite) error.
%
% Outputs:
%   steps - struct of parallel [N x 1] columns, sorted by channel then sample. An empty
%           result has 0x1 / {} / false(0,1) fields, never a bare [].
%           .channel   row index into EEG.data     .label    channel label
%           .sample    mid-level crossing (best single-sample location)
%           .onset     last sample still at the pre-level
%           .offset    first sample at the post-level (the transition proper is the open
%                      interval (onset, offset); repair replaces [onset, offset] inclusive)
%           .amplitude signed uV, post-level minus pre-level
%           .z         max local z in the run (Dslope's departure from its local trend, in
%                      units of its local robust spread)
%           .channelZ  max channel-wide slope z in the run (|Dslope| over this channel's own
%                      quiet-sample SlopeSpreadQuantile scale)
%           .riseMs    (offset - onset - 1) / srate * 1000, so a one-sample step reports 0
%                      (with RefineTransition false a one-sample gate run has onset == offset
%                      and reports -1 / srate * 1000)
%           .eventId   id of the co-occurrence event this step belongs to (info.events)
%           .commonMode true when that event spans more than CommonModeChannelFraction of
%                      the channels
%   info  - .available, .status, .params (with .derived sample counts), .nSteps,
%           .nChannelsAffected, .stepsPerChannel, .stepRatePerMinute, .transitionFraction,
%           .channelScaleSlope, .floorBoundFraction (fraction of valid samples where the local
%           robust spread s was below SlopeNoiseFloorMicrovolts), .rawScale (1.4826 * MAD of
%           the raw channel over valid samples, after subtracting each inter-boundary segment's
%           own median - a stitched recording's segments can carry different DC offsets, which
%           are not channel noise), .diffScale (the same of diff(x) restricted to
%           pairs whose samples are both valid), .zeroSlopeFraction (fraction of valid samples
%           where Dslope is exactly 0 - at short TransitionMs lags two neighbouring centred
%           medians are often equal; NaN for degenerate channels) (all nbchan x 1), .badChannels /
%           .badChannelReason (worst-first, as RELAX_channel_budget 'spend' requires;
%           re-rank with RELAX_rank_pop_channels once eye-movement steps are flagged),
%           .labels / .minutes / .pnts (what that re-ranking needs),
%           .degenerateChannels, .nRejected (counts per rejection reason), .events
%           (columns id, nChannels, lo, hi, commonMode, sameSignFraction),
%           .nCommonModeSteps, .nCommonModeEvents, .medianRiseMs.
%
% Detection runs inside each inter-boundary segment (RELAX_get_segments); a splice at a
% boundary is never scored. Each segment edge has a placement margin of W + gPre samples on the
% pre side and W + gPost samples on the post side where no candidate can even be centred, and
% segments shorter than 2W + lag - 1 are skipped. A step's own visibility reach is wider still:
% the candidate statistic Dslope compares mc(t-gPre) to mc(t+gPost), both centred median
% filters, so a step needs to sit roughly W samples into the segment - not just past the raw
% placement margin - before some valid centre's window pair still straddles it; that is about
% double the reach the superseded level gate had, which compared trailing windows instead
% (measured 2026-09-22, job 19867003: reliably found from ~140-160 samples in, vs ~80-100
% before).
%
% Known limits: a pulse shorter than ~W/2 is not a step and is not found (left to the
% extreme-period marking); a pulse of ~W/2 or longer is reported as two opposite steps; an
% exactly constant channel has
% no scale and is reported in info.degenerateChannels rather than scored. "Degenerate" means
% constant (or non-finite) within EVERY usable segment - a globally flat channel, but also a
% railed one that sits at a different constant level in each inter-boundary segment (m544
% day 2, E154: 0 / +333333 uV, 0 steps). Such a channel carries nothing; the wrapper rejects
% it through RELAX_rank_pop_channels (reason 'flat'). A rise longer than TransitionMs is
% attenuated by L/T (L = the TransitionMs lag in samples, T = the rise's own duration in
% samples): Dslope, taken across only the TransitionMs lag, sees a fraction L/T of the rise's
% full amplitude. A slow rise also masks itself: once it fills more than (1 - q) of the
% SlopeBaselineMs window, its own slope inflates the local spread s it is scored against. At
% the defaults (1000 ms, q 0.85) a linear rise of ~150 ms or longer is not found at any
% amplitude (final review, 2026-09-25: 100-125 ms ramps found, 150-500 ms never, at 200 and
% 1000 uV). That limit is set by SlopeBaselineMs and SlopeSpreadQuantile, not by TransitionMs.
%
% Per-channel work splits into a statistics stage (independent of ZThreshold/
% ChannelZThreshold: the degenerate check, levelStatistic, floorBoundFraction, rawScale,
% diffScale, zeroSlopeFraction, the smoothed trace) computed once, and a gating stage (gateCandidates onward)
% repeated for every ThresholdPairs row. Channels are distributed over any already-open
% parallel pool (gcp('nocreate')); this function never opens one itself.

function [steps, info] = RELAX_detect_steps(EEG, varargin)

    p = inputParser;
    p.addParameter('WindowMs', 151);
    p.addParameter('TransitionMs', 4);
    p.addParameter('ZThreshold', 20);
    p.addParameter('ChannelZThreshold', 35);
    p.addParameter('MinAmplitudeMicrovolts', 15);
    p.addParameter('RefractoryMs', 10);
    p.addParameter('SmoothMs', 11);
    p.addParameter('BoundaryFraction', 0.10);
    p.addParameter('BoundaryNoiseMultiple', 2);
    p.addParameter('RefineTransition', true, @isRefineTransitionValue);
    p.addParameter('StepsPerMinuteToRejectChannel', 6);
    p.addParameter('MaxTransitionFractionToRejectChannel', 0.02);
    p.addParameter('CommonModeGapMs', 10);
    p.addParameter('CommonModeChannelFraction', 0.15);
    p.addParameter('SlopeBaselineMs', 1000);
    p.addParameter('SlopeNoiseFloorMicrovolts', 0.1);
    p.addParameter('SlopeSpreadQuantile', 0.85);
    p.addParameter('ThresholdPairs', []);
    p.parse(varargin{:});
    opt = p.Results;
    opt.RefineTransition = logical(opt.RefineTransition);
    opt.derived = derivedSamples(opt, EEG.srate);
    d = opt.derived;

    pairs = validateThresholdPairs(opt.ThresholdPairs, opt.ZThreshold, opt.ChannelZThreshold);
    nPairs = size(pairs, 1);

    nb = EEG.nbchan;
    labels = {EEG.chanlocs.labels};

    steps = repmat(emptySteps(), nPairs, 1);
    infoCell = cell(nPairs, 1);
    for k = 1:nPairs
        infoCell{k} = makeInfo(opt, nb, EEG, labels, pairs(k,:));
    end
    info = vertcat(infoCell{:});

    segments = RELAX_get_segments(EEG);
    minSegLen = 2*d.W + d.lag - 1;      % one valid centre needs both windows inside
    usable = (segments(:,2) - segments(:,1) + 1) >= minSegLen;
    if ~any(usable)
        status = sprintf(['no inter-boundary segment reaches the %d samples needed for ' ...
            'two %g ms windows around a %g ms transition'], minSegLen, opt.WindowMs, opt.TransitionMs);
        for k = 1:nPairs
            info(k).status = status;
        end
        return
    end
    segments = segments(usable, :);
    valid = validMask(segments, EEG.pnts, d.W, d.gPre, d.gPost);

    data = double(EEG.data);
    srate = EEG.srate;
    pool = gcp('nocreate');
    if isempty(pool); M = 0; else; M = pool.NumWorkers; end

    chanResults = cell(nb, 1);
    parfor (ch = 1:nb, M)
        chanResults{ch} = detectStepsOneChannel(data(ch,:), segments, valid, d, opt, pairs, srate); %#ok<PFBNS>
    end

    rowsAll = cell(nPairs, 1);
    rej = cell(nPairs, 1);
    for k = 1:nPairs
        rowsAll{k} = {};
        rej{k} = emptyRejected();
    end
    degenerate = cell(0,1);

    for ch = 1:nb
        r = chanResults{ch};
        if r.degenerate
            degenerate{end+1,1} = labels{ch}; %#ok<AGROW>
            continue
        end
        for k = 1:nPairs
            info(k).channelScaleSlope(ch) = r.channelScaleSlope(k);
            info(k).floorBoundFraction(ch) = r.floorBoundFraction;
            info(k).rawScale(ch) = r.rawScale;
            info(k).diffScale(ch) = r.diffScale;
            info(k).zeroSlopeFraction(ch) = r.zeroSlopeFraction;
            for f = fieldnames(r.rej{k})'
                rej{k}.(f{1}) = rej{k}.(f{1}) + r.rej{k}.(f{1});
            end
            for m = 1:numel(r.rows{k})
                row = r.rows{k}{m};
                row.channel = ch;
                rowsAll{k}{end+1} = row; %#ok<AGROW>
            end
        end
    end

    for k = 1:nPairs
        info(k).degenerateChannels = degenerate;
        info(k).nRejected = rej{k};
        info(k).available = true;
        info(k).status    = 'ok';
        if isempty(rowsAll{k}); continue; end

        s = emptySteps();
        tbl = [rowsAll{k}{:}];
        s.channel   = [tbl.channel]';
        s.sample    = [tbl.sample]';
        s.onset     = [tbl.onset]';
        s.offset    = [tbl.offset]';
        s.amplitude = [tbl.amplitude]';
        s.z         = [tbl.z]';
        s.channelZ  = [tbl.channelZ]';
        s.riseMs    = [tbl.riseMs]';
        [~, ord] = sortrows([s.channel s.sample]);
        for f = fieldnames(s)'
            if ~isempty(s.(f{1})); s.(f{1}) = s.(f{1})(ord); end
        end
        s.label = reshape(labels(s.channel), [], 1);

        % co-occurrence events and common-mode flag
        [s.eventId, s.commonMode, info(k).events] = groupEvents(s, nb, d.cmGap, ...
            opt.CommonModeChannelFraction);
        info(k).nCommonModeSteps  = sum(s.commonMode);
        info(k).nCommonModeEvents = sum(info(k).events.commonMode);

        info(k).nSteps            = numel(s.sample);
        info(k).stepsPerChannel   = accumarray(s.channel, 1, [nb 1]);
        info(k).nChannelsAffected = sum(info(k).stepsPerChannel > 0);
        info(k).medianRiseMs      = median(s.riseMs);

        % bad contacts, from the channel's own (non-common-mode) steps; the wrapper re-ranks
        % with RELAX_rank_pop_channels once eye-movement steps are known
        info(k) = RELAX_rank_pop_channels(s, info(k), ~s.commonMode);
        steps(k) = s;
    end
end

% ------------------------------------------------------------------------------------------
function pairs = validateThresholdPairs(tp, zDefault, czDefault)
    if isequal(size(tp), [0 0])
        pairs = [zDefault czDefault];
        return
    end
    if ~isnumeric(tp) || ~ismatrix(tp) || size(tp,2) ~= 2 || size(tp,1) < 1
        error('RELAX_detect_steps:ThresholdPairs', ...
            'ThresholdPairs must be an N x 2 numeric matrix with N >= 1');
    end
    if ~all(isfinite(tp(:)))
        error('RELAX_detect_steps:ThresholdPairs', 'ThresholdPairs must contain only finite values');
    end
    pairs = tp;
end

function tf = isRefineTransitionValue(v)
% RefineTransition is a logical scalar: true (bracket once) or false (raw candidate). A 0/1
% numeric scalar is accepted; anything else, including the retired string 'once', is not.
    tf = isscalar(v) && (islogical(v) || (isnumeric(v) && (v == 0 || v == 1)));
end

function info = makeInfo(opt, nb, EEG, labels, pair)
    optK = opt;
    optK.ZThreshold = pair(1);
    optK.ChannelZThreshold = pair(2);
    info = struct('available', false, 'status', '', 'params', optK, ...
                   'nSteps', 0, 'nChannelsAffected', 0, ...
                   'stepsPerChannel', zeros(nb,1), 'stepRatePerMinute', zeros(nb,1), ...
                   'transitionFraction', zeros(nb,1), ...
                   'channelScaleSlope', nan(nb,1), ...
                   'floorBoundFraction', nan(nb,1), 'rawScale', nan(nb,1), 'diffScale', nan(nb,1), ...
                   'zeroSlopeFraction', nan(nb,1), ...
                   'badChannels', {cell(0,1)}, 'badChannelReason', {cell(0,1)}, ...
                   'degenerateChannels', {cell(0,1)}, ...
                   'labels', {reshape(labels, [], 1)}, ...
                   'minutes', (EEG.pnts / EEG.srate) / 60, 'pnts', EEG.pnts, ...
                   'nRejected', emptyRejected(), ...
                   'events', emptyEvents(), ...
                   'nCommonModeSteps', 0, 'nCommonModeEvents', 0, 'medianRiseMs', NaN);
end

function r = detectStepsOneChannel(x, segments, valid, d, opt, pairs, srate)
% One channel's statistics stage (once) plus gating stage (per ThresholdPairs row).
    nPairs = size(pairs, 1);
    r.degenerate = isDegenerate(x, segments);
    r.channelScaleSlope = nan(1, nPairs);
    r.floorBoundFraction = NaN;
    r.rawScale = NaN;
    r.diffScale = NaN;
    r.zeroSlopeFraction = NaN;
    r.rows = cell(nPairs, 1);
    r.rej  = cell(nPairs, 1);
    for k = 1:nPairs
        r.rows{k} = {};
        r.rej{k}  = emptyRejected();
    end
    if r.degenerate
        return
    end

    [D, Dslope, z, s] = levelStatistic(x, d.W, d.gPre, d.gPost, valid, d.baseline, ...
                                    opt.SlopeNoiseFloorMicrovolts, opt.SlopeSpreadQuantile);
    r.floorBoundFraction = mean(s(valid) < opt.SlopeNoiseFloorMicrovolts);
    r.rawScale = robustSd(x(valid) - segmentMedians(x, segments, valid));
    dx = diff(x);
    r.diffScale = robustSd(dx(valid(1:end-1) & valid(2:end)));
    r.zeroSlopeFraction = mean(Dslope(valid) == 0);

    xs = [];
    for k = 1:nPairs
        optK = opt;
        optK.ZThreshold = pairs(k,1);
        optK.ChannelZThreshold = pairs(k,2);
        [cand, chanZ, scaleSlope] = gateCandidates(D, Dslope, z, valid, optK);
        r.channelScaleSlope(k) = scaleSlope;
        if ~any(cand); continue; end

        if isempty(xs)
            xs = movmedian(x, d.smooth, 'omitnan');
        end
        for sgn = [1 -1]
            [rs, re] = signedRuns(cand, Dslope, sgn, d.refr);
            for m = 1:numel(rs)
                segIdx = find(segments(:,1) <= rs(m) & segments(:,2) >= re(m), 1);
                if isempty(segIdx); continue; end
                if optK.RefineTransition
                    [row, why] = refineCandidate(x, xs, z, chanZ, rs(m), re(m), ...
                                                 segments(segIdx,1), segments(segIdx,2), optK, srate);
                else
                    [row, why] = rawCandidate(x, D, z, chanZ, rs(m), re(m), ...
                                              segments(segIdx,1), segments(segIdx,2), optK, srate);
                end
                if isempty(row)
                    r.rej{k}.(why) = r.rej{k}.(why) + 1;
                else
                    r.rows{k}{end+1} = row; %#ok<AGROW>
                end
            end
        end
        if optK.RefineTransition
            [r.rows{k}, r.rej{k}] = mergeSameSignRows(r.rows{k}, x, segments, d.W, srate, ...
                                                       optK.MinAmplitudeMicrovolts, r.rej{k});
        end
    end
end

% ------------------------------------------------------------------------------------------
function [rows, rej] = mergeSameSignRows(rows, x, segments, W, srate, minAmp, rej)
% Merge this channel's overlapping brackets (RELAX_merge_overlapping_steps), positive- and
% negative-amplitude rows separately so a genuine up/down excursion is never cancelled.
    if isempty(rows); return; end
    tbl = [rows{:}];
    amp = [tbl.amplitude];
    rows = {};
    for sgn = [1 -1]
        idx = find(sign(amp) == sgn);
        if isempty(idx); continue; end
        sub = tbl(idx);
        s = struct('channel', [sub.channel]', 'sample', [sub.sample]', 'onset', [sub.onset]', ...
                    'offset', [sub.offset]', 'amplitude', [sub.amplitude]', 'z', [sub.z]', ...
                    'channelZ', [sub.channelZ]', 'riseMs', [sub.riseMs]');
        merged = RELAX_merge_overlapping_steps(s, x, segments, W, srate);
        for i = 1:numel(merged.sample)
            if abs(merged.amplitude(i)) < minAmp
                rej.amplitude = rej.amplitude + 1;
                continue
            end
            rows{end+1} = struct('channel', merged.channel(i), 'sample', merged.sample(i), ...
                'onset', merged.onset(i), 'offset', merged.offset(i), ...
                'amplitude', merged.amplitude(i), 'z', merged.z(i), ...
                'channelZ', merged.channelZ(i), 'riseMs', merged.riseMs(i)); %#ok<AGROW>
        end
    end
end

% ------------------------------------------------------------------------------------------
function tf = isDegenerate(x, segments)
% True when the channel has no variation inside any usable segment.
    tf = RELAX_flat_channels(x, segments);
end

function s = emptySteps()
    z = zeros(0,1);
    s = struct('channel', z, 'label', {cell(0,1)}, 'sample', z, 'onset', z, 'offset', z, ...
               'amplitude', z, 'z', z, 'channelZ', z, 'riseMs', z, 'eventId', z, ...
               'commonMode', false(0,1));
end

function r = emptyRejected()
    r = struct('amplitude', 0, 'edge', 0);
end

function e = emptyEvents()
    z = zeros(0,1);
    e = struct('id', z, 'nChannels', z, 'lo', z, 'hi', z, 'commonMode', false(0,1), ...
               'sameSignFraction', z);
end

function d = derivedSamples(opt, srate)
% Sample counts used everywhere else; W and smooth are forced odd.
    ms2n = @(ms) max(1, round(ms * srate / 1000));
    d.W = ms2n(opt.WindowMs);      if mod(d.W, 2) == 0;      d.W = d.W + 1;           end
    d.smooth = ms2n(opt.SmoothMs); if mod(d.smooth, 2) == 0; d.smooth = d.smooth + 1; end
    d.lag   = ms2n(opt.TransitionMs);      % Dslope's lag in samples, exactly
    d.gPre  = floor(d.lag / 2);
    d.gPost = d.lag - d.gPre;              % a one-sample lag reads mc(t+1) - mc(t)
    d.h      = floor(d.W / 2);
    d.refr   = ms2n(opt.RefractoryMs);
    d.cmGap  = ms2n(opt.CommonModeGapMs);
    d.baseline = ms2n(opt.SlopeBaselineMs);
    if mod(d.baseline, 2) == 0; d.baseline = d.baseline + 1; end
end

function valid = validMask(segments, n, W, gPre, gPost)
% true where both level windows fit inside the segment holding t
    valid = false(1, n);
    for s = 1:size(segments,1)
        lo = segments(s,1) + gPre + W - 1;
        hi = segments(s,2) - gPost - W + 1;
        if hi >= lo; valid(lo:hi) = true; end
    end
end

function [D, Dslope, z, s] = levelStatistic(x, W, gPre, gPost, valid, baselineN, noiseFloor, spreadQ)
% D(t): the settled level change across a TransitionMs + WindowMs lag, which sets the sign and
% meets the absolute floor. Dslope(t): the same centred median filter differenced across the
% transition allowance alone, so a level change is scored by how fast it happened. b(t): the
% local linear trend in Dslope (a SlopeBaselineMs-long signed median) - monotonic drift sets
% its own Dslope trend and b tracks it, so a step riding on steady drift is scored against
% that drift's slope, not against zero. s(t): the local robust spread of Dslope around that
% trend - oscillating drift is not removed by b (its median slope is near zero) but is
% absorbed into s, which widens to match its swing. z(t): Dslope's departure from the local
% trend b, in units of the local spread s. NaN outside valid.
    n = numel(x);
    h = floor(W / 2);
    mPre  = movmedian(x, [W-1 0], 'omitnan');       % median of x(t-W+1 .. t)
    mPost = movmedian(x, [0 W-1], 'omitnan');       % median of x(t .. t+W-1)
    mc    = [mPre(h+1:end) nan(1, h)];              % centred median filter: mc(t) = mPre(t+h)
    D = nan(1, n); Dslope = nan(1, n); z = nan(1, n); s = nan(1, n);
    t = find(valid);
    D(t)      = mPost(t + gPost) - mPre(t - gPre);
    Dslope(t) = mc(t + gPost) - mc(t - gPre);
    b = movmedian(Dslope, baselineN, 'omitnan');
    sAll = movingQuantileSd(Dslope - b, baselineN, spreadQ);
    s(t) = sAll(t);
    z(t) = abs(Dslope(t) - b(t)) ./ max(s(t), noiseFloor);
end

function [cand, chanZ, scaleSlope] = gateCandidates(D, Dslope, z, valid, opt)
    loc = valid & z > opt.ZThreshold & abs(D) >= opt.MinAmplitudeMicrovolts;
    Dq = Dslope(valid & ~loc);
    Dq = Dq(isfinite(Dq));
    if isempty(Dq)
        scaleSlope = NaN;
    else
        scaleSlope = quantileSd(Dq - median(Dq), opt.SlopeSpreadQuantile);
    end
    if ~(scaleSlope > 0)
        % A channel whose quiet samples have no spread: the channel gate has nothing to say.
        scaleSlope = eps;
    end
    chanZ = abs(Dslope) / scaleSlope;
    cand = loc & chanZ > opt.ChannelZThreshold;
end

function [rs, re] = signedRuns(cand, D, sgn, refr)
% Start/end indices of runs of candidate samples with sign(D) == sgn; gaps shorter than
% refr samples are closed.
    c = cand & sign(D) == sgn;
    dd = diff([0 c 0]);
    rs = find(dd == 1);
    re = find(dd == -1) - 1;
    if numel(rs) < 2; return; end
    keep = true(size(rs));
    for k = 2:numel(rs)
        prev = find(keep(1:k-1), 1, 'last');
        if rs(k) - re(prev) - 1 < refr
            re(prev) = re(k);
            keep(k) = false;
        end
    end
    rs = rs(keep); re = re(keep);
end

function [row, why] = refineCandidate(x, xs, z, chanZ, rs, re, a, b, opt, srate)
% Bracket the transition around one candidate run, once. Returns [] and a reason when the run
% is not a settled level change.
    d = opt.derived;
    row = []; why = '';

    % the run is a plateau of the gate statistic, so its geometric centre says nothing about
    % where the edge is; its strongest sample does
    [~, kk] = max(z(rs:re));
    mid0 = rs + kk - 1;
    searchLo = max(a, rs - (d.gPre + d.h));
    searchHi = min(b, re + (d.gPost + d.h));

    % initial levels from the windows the statistic itself used at the run midpoint
    Lpre  = median(x(mid0 - d.gPre - d.W + 1 : mid0 - d.gPre), 'omitnan');
    Lpost = median(x(mid0 + d.gPost : mid0 + d.gPost + d.W - 1), 'omitnan');

    A = Lpost - Lpre;
    if ~(abs(A) >= opt.MinAmplitudeMicrovolts); why = 'amplitude'; return; end
    half = (Lpre + Lpost) / 2;

    % mid-level crossing nearest the run midpoint; nothing found falls back to the strongest
    % sample rather than rejecting
    seg = xs(searchLo:searchHi) - half;
    sgnA = sign(A);
    crossed = find(sgnA * seg(2:end) >= 0 & sgnA * seg(1:end-1) < 0) + searchLo; % first post-side sample
    if isempty(crossed)
        tmid = mid0;                      % nothing to lock onto: stay at the strongest sample
    else
        [~, kk] = min(abs(crossed - mid0));
        tmid = crossed(kk);
    end

    % settle bands: a fraction of the amplitude, floored by the smoothed noise
    preWin  = mid0 - d.gPre - d.W + 1 : mid0 - d.gPre;
    postWin = mid0 + d.gPost : mid0 + d.gPost + d.W - 1;
    bandPre  = max(opt.BoundaryFraction * abs(A), opt.BoundaryNoiseMultiple * robustSd(xs(preWin)));
    bandPost = max(opt.BoundaryFraction * abs(A), opt.BoundaryNoiseMultiple * robustSd(xs(postWin)));

    % the walk reaches the segment edge and falls back to the run's own end rather than
    % rejecting: every candidate the gate accepted survives
    k = find(abs(xs(a:max(a, tmid-1)) - Lpre) <= bandPre, 1, 'last');
    if isempty(k)
        onset = min(rs, tmid);
    else
        onset = a + k - 1;
    end
    k = find(abs(xs(tmid:b) - Lpost) <= bandPost, 1, 'first');
    if isempty(k)
        offset = max(re, tmid);
    else
        offset = tmid + k - 1;
    end

    % re-measure the levels from windows adjacent to the bracket
    preWin  = max(a, onset - d.W + 1) : onset;
    postWin = offset : min(b, offset + d.W - 1);
    if numel(preWin) < ceil(d.W/2) || numel(postWin) < ceil(d.W/2); why = 'edge'; return; end
    Lpre  = median(x(preWin), 'omitnan');
    Lpost = median(x(postWin), 'omitnan');

    A = Lpost - Lpre;
    if ~(abs(A) >= opt.MinAmplitudeMicrovolts); why = 'amplitude'; return; end

    row = struct('channel', 0, 'sample', tmid, 'onset', onset, 'offset', offset, ...
                 'amplitude', A, 'z', max(z(rs:re)), 'channelZ', max(chanZ(rs:re)), ...
                 'riseMs', (offset - onset - 1) / srate * 1000);
end

function [row, why] = rawCandidate(~, D, z, chanZ, rs, re, a, b, opt, srate)
% The candidate run itself, with no bracketing: onset and offset are the run's own ends and the
% step is placed at its strongest sample. Bracketing measures settlement against a level taken
% from a window outside the bracket, so any baseline drift adjacent to a step is absorbed into
% the transition and the bracket walks outward until it settles - a fast step on a sloping
% baseline ends up jammed against one edge of a bracket many times its length. This mode gives
% up the transition estimate entirely rather than report a contaminated one: the interval is
% the span over which the gate stayed lit, which for an ideal step is about TransitionMs wide
% and centred on the edge.
    d = opt.derived;
    row = []; why = '';
    [~, kk] = max(z(rs:re));
    tmid = rs + kk - 1;
    A = D(tmid);
    if ~(abs(A) >= opt.MinAmplitudeMicrovolts); why = 'amplitude'; return; end

    preWin  = max(a, tmid - d.gPre - d.W + 1) : max(a, tmid - d.gPre);
    postWin = min(b, tmid + d.gPost) : min(b, tmid + d.gPost + d.W - 1);
    if numel(preWin) < ceil(d.W/2) || numel(postWin) < ceil(d.W/2); why = 'edge'; return; end

    row = struct('channel', 0, 'sample', tmid, 'onset', rs, 'offset', re, ...
                 'amplitude', A, 'z', max(z(rs:re)), 'channelZ', max(chanZ(rs:re)), ...
                 'riseMs', (re - rs - 1) / srate * 1000);
end

function m = segmentMedians(x, segments, valid)
% Each valid sample's own segment median, so a DC offset that differs between inter-boundary
% segments does not read as channel noise.
    m = nan(1, numel(x));
    for k = 1:size(segments, 1)
        idx = segments(k,1):segments(k,2);
        m(idx) = median(x(idx(valid(idx))), 'omitnan');
    end
    m = m(valid);
end

function s = quantileSd(v, q)
% Robust sd from the q-quantile of |v|, scaled to the Gaussian: q = 0.5 is 1.4826 * MAD. A high q
% stays finite where more than half of v is exactly zero, which a short-lag median-filter
% difference is (neighbouring centred medians are equal whenever the entering and leaving samples
% fall on the same side of the median).
    v = abs(v(isfinite(v)));
    if isempty(v); s = NaN; return; end
    s = quantile(v, q) / sqrt(2) / erfinv(q);
end

function s = movingQuantileSd(v, N, q)
% quantileSd over a centred N-sample window at every sample, NaN ignored. Evaluated every
% stride = N/20 samples and held in between: a spread over N samples moves little within N/20.
    n = numel(v);
    stride = max(1, floor(N / 20));
    h = floor(N / 2);
    centres = 1:stride:n;
    a = abs(v);
    vals = nan(1, numel(centres));
    offs = (-h:h)';
    chunk = max(1, floor(2e7 / numel(offs)));
    for c0 = 1:chunk:numel(centres)
        cc = centres(c0 : min(numel(centres), c0 + chunk - 1));
        idx = offs + cc;
        inside = idx >= 1 & idx <= n;
        idx(~inside) = 1;
        M = a(idx);
        M(~inside) = NaN;
        M = sort(M, 1);
        cnt = sum(~isnan(M), 1);
        k = max(1, ceil(q * cnt));
        pick = M(sub2ind(size(M), k, 1:numel(cc)));
        pick(cnt == 0) = NaN;
        vals(c0 : c0 + numel(cc) - 1) = pick;
    end
    held = vals(min(numel(centres), floor(((1:n) - 1) / stride) + 1));
    s = held / sqrt(2) / erfinv(q);
end

function s = robustSd(v)
    v = v(isfinite(v));
    if isempty(v); s = 0; return; end
    s = 1.4826 * median(abs(v - median(v)));
end

function [eventId, commonMode, events] = groupEvents(steps, nb, gapSamples, frac)
% Steps within gapSamples of the previous step (by .sample, across all channels) chain into
% one event; a gap larger than gapSamples starts a new one. eventId is the event's rank in
% time order (1..nEv, by first sample).
    [sortedSample, ord] = sort(steps.sample);
    newEvent = [true; diff(sortedSample) > gapSamples];
    sortedEventId = cumsum(newEvent);
    eventId = zeros(size(steps.sample));
    eventId(ord) = sortedEventId;
    nEv = max(eventId);
    events = emptyEvents();
    events.id = (1:nEv)';
    events.nChannels = zeros(nEv,1); events.lo = zeros(nEv,1); events.hi = zeros(nEv,1);
    events.sameSignFraction = zeros(nEv,1);
    for e = 1:nEv
        m = eventId == e;
        events.nChannels(e) = numel(unique(steps.channel(m)));
        events.lo(e) = min(steps.onset(m));
        events.hi(e) = max(steps.offset(m));
        pos = mean(steps.amplitude(m) > 0);
        events.sameSignFraction(e) = max(pos, 1 - pos);
    end
    % at least two channels, whatever the fraction says on a tiny montage
    events.commonMode = events.nChannels > max(frac * nb, 1);
    commonMode = events.commonMode(eventId);
end
