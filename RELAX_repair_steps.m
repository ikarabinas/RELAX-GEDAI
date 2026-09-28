% RELAX_repair_steps() - Repair electrode-pop transitions: detrend between them, bridge them.
%
% Usage:
%   >> [EEG, info] = RELAX_repair_steps(EEG, steps);
%   >> [EEG, info] = RELAX_repair_steps(EEG, steps, 'FillMethod', 'linear');
%
% Three passes:
%
%  1. DETREND. Each channel is partitioned at the union of {inter-boundary segment starts},
%     {its own transition intervals [onset, offset]} and {the common-mode excision
%     intervals, for every channel}. Every sub-segment between those cuts is detrended in
%     complete isolation with the vendored NoiseTools nt_detrend (no padding, no overlap
%     across a cut). nt_regw's returned trend CARRIES THE WEIGHTED MEAN, so `x - trend` has
%     zero weighted mean over the fit interval: detrending each inter-transition sub-segment
%     removes the level offset as a side effect, and no transition sample ever sits inside a
%     fit window.
%
%  2. FILL. Each transition interval is replaced. The base is a LINEAR BRIDGE between the
%     mean of the last AnchorMs of the detrended left neighbour and the mean of the first
%     AnchorMs of the detrended right neighbour (a constant hold when the interval touches a
%     segment edge). With FillMethod 'spatial' (default) the neighbours' DETAIL is added on
%     top: every basis channel j gets its own chord between the same two anchor windows, and
%     the spherical-spline estimate (RELAX_spline_weights, identical to EEGLAB eeg_interp)
%     of sum_j w_j * (y_j - chord_j) is added to the bridge. Interpolating the
%     anchor-relative detail rather than the raw neighbours keeps the raw DC offsets out of
%     the fill and keeps the joins at anchor-noise level; when the neighbours carry nothing
%     the fill degrades to the plain bridge. Channels in transition themselves over the
%     interval (co-popping neighbours), channels without coordinates, and channels in
%     BasisExcludeLabels are left out of the basis; with fewer than MinBasisChannels left
%     the interval falls back to the linear bridge (counted in info.nFillsLinear).
%
%     What the bridge does and does not fix: it removes the level DISCONTINUITY that a
%     zero-phase highpass would ring on for seconds; the two neighbours remain
%     independently zero-meaned, so continuity of the local signal across the interval is
%     only as good as the spatial detail estimate.
%
%  3. EXCISE. Steps flagged commonMode (a whole-head shift on more than a set fraction of
%     channels) are not filled and, by default (ExciseCommonMode = false), are left in
%     place entirely: they are excluded from detrending's partition points as well as from
%     the fill, so a channel with only common-mode steps is repaired identically to one with
%     no steps at all. With ExciseCommonMode true the union of each such event's intervals
%     is instead deleted from ALL channels through eeg_eegrej, which inserts a 'boundary'
%     event of that duration; events longer than MaxExciseMs are neither excised nor filled
%     even then (counted in info.nCommonModeTooLong) and are left to the later extreme-period
%     marking. Steps flagged eyeMovement (RELAX_flag_eye_steps) are never touched, even
%     when their event is also commonMode: a rim-wide shift that coincides with a blink
%     or saccade is EOG the later ICA needs, not head movement.
%
% Inputs:
%   EEG   - EEGLAB EEG structure (uses .data, .nbchan, .pnts, .srate, .event, .chanlocs)
%   steps - struct from RELAX_detect_steps (.channel .onset .offset, optional .commonMode
%           .eyeMovement .eventId). A legacy struct with only .sample is accepted: each
%           step becomes a zero-length cut at .sample (the old single-sample partition,
%           nothing filled; info.legacyIntervals = true). Pass one with empty columns to
%           apply drift-only detrending with no cuts.
%   Name-value pairs (all optional):
%     MaxPolyOrder      - highest polynomial order used (default 3)
%     MaxCutoffHz       - the invariant: this stage never attenuates above this frequency
%                          (default 1.0, just under the pipeline's own 1 Hz highpass)
%     MaxWindowSeconds  - nt_detrend's overlap-add window for long sub-segments (default 2.0)
%     DetrendThresh     - nt_detrend outlier threshold in sd (default 3)
%     DetrendIterations - nt_detrend reweighting iterations (default 3)
%     MergeSamples      - intervals separated by at most this many samples are merged; an
%                          interval (or legacy cut) within this many samples of a segment
%                          edge is extended into the edge (dropped, for a cut) rather than
%                          leaving a sliver sub-segment (default 20)
%     AnchorMs          - length of the anchor windows either side of a fill (default 3;
%                          the mean of 3 samples, not a 5 ms median, avoids a slope bias
%                          of 1-5 uV at the join)
%     FillMethod        - 'spatial' (default) or 'linear'
%     BasisExcludeLabels - labels never used as a spatial basis (default {}; the wrapper
%                          passes degenerate and over-impedance channels)
%     MinBasisChannels  - below this many usable basis channels a fill is linear (default 16)
%     ExciseCommonMode  - excise common-mode events (default false; a common-mode step is
%                          left in place otherwise - neither excised, filled, nor used as a
%                          partition point)
%     MaxExciseMs       - longest event that is excised (default 500)
%     MinSegmentSeconds - a common-mode excision is extended to a segment edge, or merged with
%                         the next one, when it would otherwise leave a shorter remnant
%                         (default 1.0; 0 = exact cuts)
%
% Outputs:
%   EEG  - same structure with .data repaired; when anything was excised, .pnts / .event /
%          .xmax are as eeg_eegrej leaves them, and every sample index in `steps` refers
%          to the PRE-excision timeline (info.excised lists what was removed).
%   info - .available, .status, .params, .nChannelsRepaired, .nSubSegments,
%          .subSegmentsPerChannel, .treatmentCounts (order3oa/order2oa/order1oa/order3/
%          order2/order1/meanOnly/constant/fallback), .minSubSegmentDuration,
%          .medianSubSegmentDuration, .fractionSamplesMeanOnly, .lowFreqVarianceRemoved,
%          .nIntervals, .intervalsPerChannel, .fractionSamplesInterpolated,
%          .fractionInterpolatedPerChannel, .maxChannelFractionInterpolated,
%          .medianTransitionMs, .nFillsSpatial, .nFillsLinear, .medianBasisSize,
%          .nCommonModeSkipped, .nEyeMovementSkipped, .nCommonModeTooLong, .excised
%          ([k x 2] pre-excision sample ranges), .nExcised, .excisedSeconds,
%          .legacyIntervals
%
% TREATMENT LADDER (unchanged). For a sub-segment of duration T, the largest order
% n <= MaxPolyOrder satisfying c(n)*n/(2T) <= MaxCutoffHz is used, with c = [1.69 1.24 1.09]
% MEASURED for orders 1-3. At the defaults:
%
%     T > 2.000 s           order 3, wsize = 2.0 s overlap-add     0.84 Hz
%     1.635 <= T <= 2.000   order 3, single window                 0.82-1.00 Hz
%     1.240 <= T < 1.635    order 2, single window                 0.76-1.00 Hz
%     0.845 <= T < 1.240    order 1, single window                 0.68-1.00 Hz
%     T < 0.845 s           robust mean subtraction only           DC only
%
% Cuts at boundaries, transitions and excisions are HARD; subdividing a long sub-segment to
% respect MaxWindowSeconds is SOFT (nt_detrend's own overlap-add). Channels with no steps are
% still partitioned by boundaries and excisions and still detrended, so every channel receives
% the same drift treatment.

function [EEG, info] = RELAX_repair_steps(EEG, steps, varargin)

    p = inputParser;
    p.addParameter('MaxPolyOrder', 3);
    p.addParameter('MaxCutoffHz', 1.0);
    p.addParameter('MaxWindowSeconds', 2.0);
    p.addParameter('DetrendThresh', 3);
    p.addParameter('DetrendIterations', 3);
    p.addParameter('MergeSamples', 20);
    p.addParameter('AnchorMs', 3);
    p.addParameter('FillMethod', 'spatial');
    p.addParameter('BasisExcludeLabels', {});
    p.addParameter('MinBasisChannels', 16);
    p.addParameter('ExciseCommonMode', false);
    p.addParameter('MaxExciseMs', 500);
    p.addParameter('MinSegmentSeconds', 1.0);
    p.parse(varargin{:});
    opt = p.Results;

    C = [1.69 1.24 1.09];   % measured f_3dB = C(n)*n/(2T); see the spec's response tables

    if opt.MaxPolyOrder > numel(C)
        warning('RELAX_repair_steps:maxPolyOrderClamped', ...
            ['MaxPolyOrder=%d exceeds the %d measured response coefficients in C, so the ' ...
             'effective maximum order is %d. Extend C with a measured coefficient if a ' ...
             'higher order is really wanted.'], opt.MaxPolyOrder, numel(C), numel(C));
    end

    nb = EEG.nbchan;
    counts = struct('order3oa', 0, 'order2oa', 0, 'order1oa', 0, ...
                    'order3', 0, 'order2', 0, 'order1', 0, 'meanOnly', 0, 'constant', 0, 'fallback', 0);
    info = struct('available', false, 'status', '', 'params', opt, ...
                  'nChannelsRepaired', 0, 'nSubSegments', 0, ...
                  'subSegmentsPerChannel', zeros(nb,1), 'treatmentCounts', counts, ...
                  'minSubSegmentDuration', NaN, 'medianSubSegmentDuration', NaN, ...
                  'fractionSamplesMeanOnly', 0, 'lowFreqVarianceRemoved', zeros(nb,1), ...
                  'nIntervals', 0, 'intervalsPerChannel', zeros(nb,1), ...
                  'fractionSamplesInterpolated', 0, 'fractionInterpolatedPerChannel', zeros(nb,1), ...
                  'maxChannelFractionInterpolated', 0, 'medianTransitionMs', NaN, ...
                  'nFillsSpatial', 0, 'nFillsLinear', 0, 'medianBasisSize', NaN, ...
                  'nCommonModeSkipped', 0, 'nEyeMovementSkipped', 0, 'nCommonModeTooLong', 0, ...
                  'excised', zeros(0,2), 'nExcised', 0, 'excisedSeconds', 0, ...
                  'legacyIntervals', false);

    if isempty(which('nt_detrend'))
        info.status = 'nt_detrend not on the path - call RELAX_add_vendored_noisetools()';
        return
    end

    segments  = RELAX_get_segments(EEG);
    srate     = EEG.srate;
    pnts      = EEG.pnts;
    wsize     = 2 * floor(round(opt.MaxWindowSeconds * srate) / 2);  % samples, even
    anchorN   = max(1, round(opt.AnchorMs * srate / 1000));
    origClass = class(EEG.data);
    data      = double(EEG.data);
    Y         = data;

    % ---- what to fill, what to skip, what to excise ----
    [ints, isLegacy, nCM, nEye] = stepIntervals(steps, nb);
    info.legacyIntervals = isLegacy;
    canExcise = opt.ExciseCommonMode && ~isempty(which('eeg_eegrej'));
    if opt.ExciseCommonMode && ~canExcise
        warning('RELAX_repair_steps:noEegrej', ...
            'eeg_eegrej is not on the path; common-mode events are skipped instead of excised');
    end
    [excise, nTooLong, nCMskipped] = exciseIntervals(steps, segments, canExcise, ...
        round(opt.MaxExciseMs * srate / 1000), round(opt.MinSegmentSeconds * srate));
    info.nCommonModeSkipped  = nCM * ~canExcise + nCMskipped * canExcise;
    info.nEyeMovementSkipped = nEye;
    info.nCommonModeTooLong  = nTooLong;

    % ---- pass 1: partition and detrend ----
    % Every channel's sub-segments are detrended independently of every other channel's, so
    % this loop is parfor'd: the dominant cost is the nt_detrend calls inside
    % detrendSubSegment, called once per sub-segment per channel. counts/allDurations can't
    % be accumulated directly in a parfor body (a struct field update and a growing array are
    % both unsupported reduction forms), so each worker fills a per-channel slot and the
    % results are merged once the loop finishes.
    fills           = cell(nb, 1);   % per channel [lo hi] to bridge
    subSegCounts    = zeros(nb, 1);
    durationsPerCh  = cell(nb, 1);
    meanOnlyPerCh   = zeros(nb, 1);
    countsPerCh     = cell(nb, 1);
    parfor ch = 1:nb
        chCounts = counts;
        [subs, fillsCh] = partitionChannel(ints{ch}, excise, segments, opt.MergeSamples);
        fills{ch} = fillsCh;
        subSegCounts(ch) = size(subs, 1);
        durs = zeros(size(subs, 1), 1);
        meanOnly = 0;
        yRow = data(ch, :);          % parfor only allows sliced Y(ch,:) with a full colon,
        for k = 1:size(subs, 1)      % so accumulate into a local row and write it back once
            a = subs(k,1); b = subs(k,2);
            x = yRow(a:b)';
            T = numel(x) / srate;
            durs(k) = T;
            [y, chCounts, nMean] = detrendSubSegment(x, T, opt, C, wsize, chCounts);
            meanOnly = meanOnly + nMean;
            yRow(a:b) = y';
        end
        Y(ch, :) = yRow;
        durationsPerCh{ch} = durs;
        meanOnlyPerCh(ch)  = meanOnly;
        countsPerCh{ch}    = chCounts;
    end
    info.subSegmentsPerChannel = subSegCounts;
    allDurations  = vertcat(durationsPerCh{:});
    meanOnlySamps = sum(meanOnlyPerCh);
    counts = mergeCounts(countsPerCh, counts);

    % ---- pass 2: fill the transition intervals ----
    XYZ = channelCoordinates(EEG.chanlocs, nb);
    hasXYZ = all(isfinite(XYZ), 2);
    labels = {EEG.chanlocs.labels};
    basisOk = hasXYZ & ~ismember(labels(:), opt.BasisExcludeLabels(:));
    spatial = strcmpi(opt.FillMethod, 'spatial') && any(hasXYZ);
    % every channel's busy intervals (fills) plus the global excisions, for the co-popping test
    busy = zeros(0, 3);            % [channel lo hi]
    for ch = 1:nb
        if ~isempty(fills{ch}); busy = [busy; repmat(ch, size(fills{ch},1), 1) fills{ch}]; end %#ok<AGROW>
    end
    weightCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    basisSizes = [];
    fillLens   = [];
    for ch = 1:nb
        F = fills{ch};
        info.intervalsPerChannel(ch) = size(F, 1);
        for k = 1:size(F, 1)
            lo = F(k,1); hi = F(k,2);
            fillLens(end+1,1) = hi - lo + 1; %#ok<AGROW>
            s = find(segments(:,1) <= lo & segments(:,2) >= hi, 1);
            a = segments(s,1); b = segments(s,2);
            kL = min(anchorN, lo - a);          % anchor samples available on each side
            kR = min(anchorN, b - hi);

            basis = [];
            if spatial && hasXYZ(ch)
                inTransition = busy(busy(:,1) ~= ch & busy(:,2) <= hi + anchorN & busy(:,3) >= lo - anchorN, 1);
                cand = find(basisOk);
                cand = cand(cand ~= ch & ~ismember(cand, inTransition));
                if numel(cand) >= opt.MinBasisChannels; basis = cand(:)'; end
            end

            if isempty(basis)
                Y(ch, lo:hi) = chord(Y(ch, :), lo, hi, kL, kR);
                info.nFillsLinear = info.nFillsLinear + 1;
            else
                key = sprintf('%d,', basis);
                if isKey(weightCache, key)
                    w = weightCache(key);
                else
                    w = RELAX_spline_weights(XYZ(basis, :), XYZ(ch, :));
                    weightCache(key) = w;
                end
                detail = zeros(numel(basis), hi - lo + 1);
                for j = 1:numel(basis)
                    detail(j, :) = Y(basis(j), lo:hi) - chord(Y(basis(j), :), lo, hi, kL, kR);
                end
                detail(~isfinite(detail)) = 0;
                Y(ch, lo:hi) = chord(Y(ch, :), lo, hi, kL, kR) + w * detail;
                info.nFillsSpatial = info.nFillsSpatial + 1;
                basisSizes(end+1,1) = numel(basis); %#ok<AGROW>
            end
        end
    end

    % ---- covariates ----
    removed = data - Y;
    for ch = 1:nb
        v = var(data(ch, ~isnan(data(ch,:))));
        if v > 0
            info.lowFreqVarianceRemoved(ch) = var(removed(ch, ~isnan(removed(ch,:)))) / v;
        end
    end
    info.nIntervals = sum(info.intervalsPerChannel);
    info.fractionInterpolatedPerChannel = cellfun(@(F) sum(F(:,2) - F(:,1) + 1), fills) / pnts;
    info.fractionSamplesInterpolated = sum(info.fractionInterpolatedPerChannel) / nb;
    info.maxChannelFractionInterpolated = max([0; info.fractionInterpolatedPerChannel]);
    if ~isempty(fillLens); info.medianTransitionMs = median(fillLens) / srate * 1000; end
    if ~isempty(basisSizes); info.medianBasisSize = median(basisSizes); end
    info.nChannelsRepaired        = nb;
    info.nSubSegments             = sum(info.subSegmentsPerChannel);
    info.treatmentCounts          = counts;
    info.minSubSegmentDuration    = min(allDurations);
    info.medianSubSegmentDuration = median(allDurations);
    info.fractionSamplesMeanOnly  = meanOnlySamps / (nb * pnts);

    EEG.data = cast(Y, origClass);

    % ---- pass 3: excise common-mode events ----
    if ~isempty(excise)
        EEG = eeg_eegrej(EEG, excise);
        info.excised = excise;
        info.nExcised = size(excise, 1);
        info.excisedSeconds = sum(excise(:,2) - excise(:,1) + 1) / srate;
    end

    info.available = true;
    info.status    = 'ok';
end

% ==========================================================================================
function [ints, isLegacy, nCM, nEye] = stepIntervals(steps, nb)
% Per-channel [onset offset] rows to repair. Skipped steps are counted, not returned.
    ints = repmat({zeros(0,2)}, nb, 1);
    nCM = 0; nEye = 0;
    isLegacy = ~(isfield(steps, 'onset') && isfield(steps, 'offset'));
    if isempty(steps.channel); return; end
    ch = steps.channel(:);
    if isLegacy
        lo = steps.sample(:); hi = lo - 1;
    else
        lo = steps.onset(:); hi = steps.offset(:);
        if numel(lo) ~= numel(ch) || numel(hi) ~= numel(ch)
            error('RELAX_repair_steps:badSteps', 'steps.channel/.onset/.offset must be parallel columns');
        end
    end
    % Eye movement takes precedence over common mode: a rim-wide event that coincides
    % with a blink / saccade is EOG, which the later ICA needs intact, so it is neither
    % filled nor excised.
    keep = true(numel(ch), 1);
    if isfield(steps, 'eyeMovement') && ~isempty(steps.eyeMovement)
        em = logical(steps.eyeMovement(:)); nEye = sum(em); keep = keep & ~em;
    end
    if isfield(steps, 'commonMode') && ~isempty(steps.commonMode)
        cm = logical(steps.commonMode(:)); nCM = sum(cm & keep); keep = keep & ~cm;
    end
    keep = keep & ch >= 1 & ch <= nb;
    for c = unique(ch(keep))'
        m = keep & ch == c;
        ints{c} = [lo(m) hi(m)];
    end
end

function [excise, nTooLong, nSkipped] = exciseIntervals(steps, segments, canExcise, maxLen, minSeg)
% Union of each common-mode event's intervals, clipped to segments and merged. A cut that
% would leave less than minSeg samples between itself and a segment edge, or between itself
% and the next cut, is extended / merged: a sub-second remnant is orphaned by the _boundary
% epoching anyway, and a ~1-sample one crashes zapline's pwelch (m292 day 3).
    excise = zeros(0, 2); nTooLong = 0; nSkipped = 0;
    if ~canExcise || ~isfield(steps, 'commonMode') || isempty(steps.commonMode); return; end
    cm = logical(steps.commonMode(:));
    if isfield(steps, 'eyeMovement') && ~isempty(steps.eyeMovement)
        cm = cm & ~logical(steps.eyeMovement(:));     % eye-flagged events stay in place
    end
    if ~any(cm); return; end
    if isfield(steps, 'eventId'); ev = steps.eventId(:); else; ev = ones(numel(cm), 1); end
    rows = zeros(0, 2);
    for e = unique(ev(cm))'
        m = cm & ev == e;
        lo = min(steps.onset(m)); hi = max(steps.offset(m));
        if hi - lo + 1 > maxLen
            nTooLong = nTooLong + 1; nSkipped = nSkipped + sum(m);
        else
            rows(end+1, :) = [lo hi]; %#ok<AGROW>
        end
    end
    % clip to segments (a piece per segment), extend into edges that would be left with a
    % remnant shorter than minSeg, then merge cuts that overlap or leave such a remnant
    clipped = zeros(0, 2);
    for r = 1:size(rows, 1)
        for s = 1:size(segments, 1)
            a = segments(s,1); b = segments(s,2);
            lo = max(rows(r,1), a); hi = min(rows(r,2), b);
            if hi < lo; continue; end
            if lo - a < minSeg; lo = a; end
            if b - hi < minSeg; hi = b; end
            clipped(end+1, :) = [lo hi]; %#ok<AGROW>
        end
    end
    excise = mergeRows(sortrows(clipped), max(0, minSeg - 1));
end

function [subs, fills] = partitionChannel(own, excise, segments, mergeN)
% Sub-segments to detrend and intervals to fill for one channel.
% own: [lo hi] rows; a row with hi == lo-1 is a legacy zero-length cut.
    subs = zeros(0, 2); fills = zeros(0, 2);
    for s = 1:size(segments, 1)
        a = segments(s,1); b = segments(s,2);
        % this segment's rows: type 1 own fill, 0 legacy cut, 2 excision
        R = zeros(0, 3);
        for r = 1:size(own, 1)
            lo = own(r,1); hi = own(r,2);
            if hi < lo                          % zero-length cut at lo
                if lo > a && lo <= b && (lo - a) > mergeN && (b + 1 - lo) > mergeN
                    R(end+1, :) = [lo lo-1 0]; %#ok<AGROW>
                end
            else
                lo = max(lo, a); hi = min(hi, b);
                if hi >= lo; R(end+1, :) = [lo hi 1]; end %#ok<AGROW>
            end
        end
        for r = 1:size(excise, 1)
            lo = max(excise(r,1), a); hi = min(excise(r,2), b);
            if hi >= lo; R(end+1, :) = [lo hi 2]; end %#ok<AGROW>
        end
        R = sortrows(R);
        % merge: a non-empty row absorbs anything within mergeN; cuts never merge with cuts
        M = zeros(0, 3);
        for r = 1:size(R, 1)
            if ~isempty(M) && R(r,1) <= M(end,2) + 1 + mergeN && ~(M(end,3) == 0 && R(r,3) == 0)
                M(end,2) = max(M(end,2), R(r,2));
                M(end,3) = max(M(end,3), R(r,3));
            else
                M(end+1, :) = R(r, :); %#ok<AGROW>
            end
        end
        % extend fills that nearly touch a segment edge into it
        for r = 1:size(M, 1)
            if M(r,3) == 1
                if M(r,1) - a <= mergeN; M(r,1) = a; end
                if b - M(r,2) <= mergeN; M(r,2) = b; end
            end
        end
        % a fill covering the whole segment has no anchors: treat it as a sub-segment
        M(M(:,3) == 1 & M(:,1) == a & M(:,2) == b, :) = [];
        % complement -> sub-segments
        cur = a;
        for r = 1:size(M, 1)
            if M(r,1) > cur; subs(end+1, :) = [cur M(r,1)-1]; end %#ok<AGROW>
            cur = M(r,2) + 1;
        end
        if cur <= b; subs(end+1, :) = [cur b]; end %#ok<AGROW>
        fills = [fills; M(M(:,3) == 1, 1:2)]; %#ok<AGROW>
    end
end

function total = mergeCounts(countsPerCh, total)
% Sum the treatment-ladder struct returned by each parfor worker into one struct.
    fn = fieldnames(total);
    for c = 1:numel(countsPerCh)
        cc = countsPerCh{c};
        for i = 1:numel(fn)
            total.(fn{i}) = total.(fn{i}) + cc.(fn{i});
        end
    end
end

function M = mergeRows(R, gap)
    M = zeros(0, 2);
    for r = 1:size(R, 1)
        if ~isempty(M) && R(r,1) <= M(end,2) + 1 + gap
            M(end,2) = max(M(end,2), R(r,2));
        else
            M(end+1, :) = R(r, :); %#ok<AGROW>
        end
    end
end

function v = chord(y, lo, hi, kL, kR)
% Linear bridge across [lo, hi] between the anchor means either side; constant hold when
% only one side exists; zeros when neither does.
    L = NaN; R = NaN;
    if kL > 0; L = mean(y(lo-kL:lo-1), 'omitnan'); end
    if kR > 0; R = mean(y(hi+1:hi+kR), 'omitnan'); end
    if isnan(L) && isnan(R); L = 0; R = 0; end
    if isnan(L); L = R; end
    if isnan(R); R = L; end
    v = linspace(L, R, hi - lo + 1);
end

function XYZ = channelCoordinates(chanlocs, nb)
    XYZ = nan(nb, 3);
    if ~all(isfield(chanlocs, {'X', 'Y', 'Z'})); return; end
    for c = 1:nb
        x = chanlocs(c).X; y = chanlocs(c).Y; z = chanlocs(c).Z;
        if isscalar(x) && isscalar(y) && isscalar(z); XYZ(c, :) = [x y z]; end
    end
end

function [y, counts, nMean] = detrendSubSegment(x, T, opt, C, wsize, counts)
% The treatment ladder for one sub-segment (unchanged from the single-sample repair).
    nMean = 0;
    order = chooseOrder(T, C, opt.MaxPolyOrder, opt.MaxCutoffHz);
    if (max(x) - min(x)) == 0
        % A constant stretch has no trend and no scale for the reweighting. It is returned
        % as it is rather than demeaned to exact zero: a channel that is constant in some
        % segments (a railed contact, m544 day 2) must not come out globally flat, which
        % crashes pop_jointprob downstream. One that is constant in EVERY segment is
        % degenerate to the detector and rejected before it gets here.
        y = x;
        counts.constant = counts.constant + 1;
        return
    end
    if any(isnan(x))
        % A NaN would propagate through the whole fit.
        order = 0;
    end
    % Compare against the ACHIEVABLE ceiling: chooseOrder saturates at numel(C).
    useOverlapAdd = (order == min(opt.MaxPolyOrder, numel(C))) && (numel(x) > wsize) && (wsize > 0);
    done = false;
    while ~done
        if order == 0
            y = robustDemean(x, opt.DetrendThresh, opt.DetrendIterations);
            counts.meanOnly = counts.meanOnly + 1;
            nMean = numel(x);
            done = true;
        else
            try
                if useOverlapAdd
                    y = nt_detrend(x, order, [], [], opt.DetrendThresh, opt.DetrendIterations, wsize);
                    oaField = sprintf('order%doa', order);
                    counts.(oaField) = counts.(oaField) + 1;
                else
                    y = nt_detrend(x, order, [], [], opt.DetrendThresh, opt.DetrendIterations);
                    switch order
                        case 3, counts.order3 = counts.order3 + 1;
                        case 2, counts.order2 = counts.order2 + 1;
                        case 1, counts.order1 = counts.order1 + 1;
                    end
                end
                if any(isnan(y)); error('RELAX:nanFromDetrend', 'nt_detrend returned NaN'); end
                done = true;
            catch
                % nt_detrend's degenerate-window failure is a bare error('!'). Step DOWN
                % the ladder rather than killing a 257-channel SLURM task hours into a run.
                counts.fallback = counts.fallback + 1;
                useOverlapAdd = false;
                order = order - 1;
                if order < 1; order = 0; end
            end
        end
    end
end

function order = chooseOrder(T, C, maxOrder, maxCutoffHz)
% Largest order n <= maxOrder whose measured -3 dB point stays at or below maxCutoffHz.
% 0 means "no polynomial fits in this duration" - robust mean subtraction only. The few-ULP
% tolerance rescues exact rung boundaries (T = 1.635 s evaluates to 1.0 + 1 ULP) without
% shifting any other rung; keep it relative to maxCutoffHz.
    order = 0;
    tol = 4 * eps(maxCutoffHz);
    for n = 1:min(maxOrder, numel(C))
        if C(n) * n / (2*T) <= maxCutoffHz + tol
            order = n;
        end
    end
end

function y = robustDemean(x, thresh, niter)
% Iteratively reweighted mean subtraction: the order-0 rung, and the fallback when
% nt_detrend cannot fit (order 0 builds an empty basis nt_regw cannot regress on).
    finite = isfinite(x);
    if ~any(finite)
        y = x;
        return
    end
    w = double(finite);
    for k = 1:niter
        m = sum(w(finite) .* x(finite)) / sum(w(finite));
        r = x - m;
        s = sqrt(sum(w(finite) .* r(finite).^2) / sum(w(finite)));
        if ~isfinite(s) || s == 0; break; end
        wNew = double(abs(r) < thresh * s) .* finite;
        if ~any(wNew); break; end
        w = wNew;
    end
    m = sum(w(finite) .* x(finite)) / sum(w(finite));
    y = x - m;
end
