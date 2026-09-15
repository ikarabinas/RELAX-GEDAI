% RELAX_repair_steps() - Repair electrode-pop steps by per-sub-segment robust detrending.
%
% Usage:
%   >> [EEG, info] = RELAX_repair_steps(EEG, steps);
%   >> [EEG, info] = RELAX_repair_steps(EEG, steps, 'MaxCutoffHz', 1.0);
%
% Each channel is partitioned at the union of {inter-boundary segment starts} and {that
% channel's own step samples}. Every sub-segment is then detrended in complete isolation -
% no padding, no overlap across a cut.
%
% That isolation is the whole mechanism. nt_regw's returned trend CARRIES THE WEIGHTED
% MEAN (nt_regw.m lines 31 and 43/58: `mn = y - nt_demean(y,w)` is added back into `z`),
% so `x - z` has zero weighted mean over the fit interval. Detrending each inter-step
% sub-segment therefore removes that sub-segment's offset as a side effect, and a separate
% step-subtraction pass (nt_destep) would find nothing left to do.
%
% The step itself is never fitted: it sits ON a partition boundary, not inside a fit
% window.
%
% Inputs:
%   EEG   - EEGLAB EEG structure (uses .data, .nbchan, .pnts, .srate, .event)
%   steps - struct from RELAX_detect_steps. Pass one with empty .sample to apply
%            drift-only detrending with no step cuts.
%   Name-value pairs (all optional):
%     MaxPolyOrder      - highest polynomial order used (default 3)
%     MaxCutoffHz       - the invariant: this stage never attenuates above this frequency
%                          (default 1.0, just under the pipeline's own 1 Hz highpass, so
%                          the stage can only remove signal the highpass was going to
%                          remove anyway)
%     MaxWindowSeconds  - nt_detrend's overlap-add window for long sub-segments (default
%                          2.0, measured -3 dB 0.841 Hz at order 3). Keeps the cutoff from
%                          sliding below 0.5 Hz on a long segment.
%     DetrendThresh     - nt_detrend outlier threshold in sd (default 3)
%     DetrendIterations - nt_detrend reweighting iterations (default 3)
%     MergeSamples      - a step within this many samples of a segment edge is merged into
%                          it rather than creating a sliver sub-segment (default 20)
%
% Outputs:
%   EEG  - same structure with .data repaired
%   info - .available, .status, .params, .nChannelsRepaired, .nSubSegments,
%          .subSegmentsPerChannel, .treatmentCounts (order3oa/order2oa/order1oa/
%          order3/order2/order1/meanOnly/fallback; the *oa fields count overlap-add
%          windowed fits and name the order actually used), .minSubSegmentDuration,
%          .medianSubSegmentDuration, .fractionSamplesMeanOnly, .lowFreqVarianceRemoved
%
% TREATMENT LADDER. For a sub-segment of duration T, the largest order n <= MaxPolyOrder
% satisfying c(n)*n/(2T) <= MaxCutoffHz is used, with c = [1.69 1.24 1.09] MEASURED for
% orders 1-3. (The familiar n/(2T) rule of thumb underestimates the true -3 dB point by up
% to 70% at order 1, which is why the constants are here rather than the heuristic.) At the
% defaults that is:
%
%     T > 2.000 s           order 3, wsize = 2.0 s overlap-add     0.84 Hz
%     1.635 <= T <= 2.000   order 3, single window                 0.82-1.00 Hz
%     1.240 <= T < 1.635    order 2, single window                 0.76-1.00 Hz
%     0.845 <= T < 1.240    order 1, single window                 0.68-1.00 Hz
%     T < 0.845 s           robust mean subtraction only           DC only
%
% HARD vs SOFT SUBDIVISION. Cuts at boundaries and steps are HARD: the signal genuinely
% discontinues there. Subdividing a long sub-segment to respect MaxWindowSeconds is SOFT,
% and goes through nt_detrend's own wsize overlap-add (50% overlap, triangular weighting,
% continuous reconstructed trend). Hard-cutting a smooth stretch every 2 s would MANUFACTURE
% step discontinuities at 2 s intervals - exactly the input this stage exists to remove.
%
% Channels with no detected steps are still partitioned by boundaries and still detrended,
% so every channel receives the same drift treatment and the stage introduces no systematic
% difference between repaired and unrepaired channels.

function [EEG, info] = RELAX_repair_steps(EEG, steps, varargin)

    p = inputParser;
    p.addParameter('MaxPolyOrder', 3);
    p.addParameter('MaxCutoffHz', 1.0);
    p.addParameter('MaxWindowSeconds', 2.0);
    p.addParameter('DetrendThresh', 3);
    p.addParameter('DetrendIterations', 3);
    p.addParameter('MergeSamples', 20);
    p.parse(varargin{:});
    opt = p.Results;

    C = [1.69 1.24 1.09];   % measured f_3dB = C(n)*n/(2T); see the spec's response tables

    % chooseOrder can never return more than numel(C), so a MaxPolyOrder above that is
    % not merely ignored - before this was guarded it silently disabled overlap-add
    % entirely (the ceiling and the configured value stopped coinciding), collapsing
    % every long sub-segment to a single global fit whose cutoff sits far below the
    % 0.5 Hz floor the ladder is documented to guarantee. Say so rather than clamp
    % quietly: the caller asked for a fit this file cannot deliver.
    if opt.MaxPolyOrder > numel(C)
        warning('RELAX_repair_steps:maxPolyOrderClamped', ...
            ['MaxPolyOrder=%d exceeds the %d measured response coefficients in C, so the ' ...
             'effective maximum order is %d. Extend C with a measured coefficient if a ' ...
             'higher order is really wanted.'], opt.MaxPolyOrder, numel(C), numel(C));
    end

    % order1oa/order2oa exist so the overlap-add counter can name the order actually
    % used. They stay zero in any shipped configuration (MaxPolyOrder=3), but a
    % rederivation pass that lowers MaxPolyOrder would otherwise have its windowed
    % sub-segments misreported as order-3.
    counts = struct('order3oa', 0, 'order2oa', 0, 'order1oa', 0, ...
                    'order3', 0, 'order2', 0, 'order1', 0, ...
                    'meanOnly', 0, 'fallback', 0);
    info = struct('available', false, 'status', '', 'params', opt, ...
                  'nChannelsRepaired', 0, 'nSubSegments', 0, ...
                  'subSegmentsPerChannel', zeros(EEG.nbchan,1), ...
                  'treatmentCounts', counts, ...
                  'minSubSegmentDuration', NaN, 'medianSubSegmentDuration', NaN, ...
                  'fractionSamplesMeanOnly', 0, ...
                  'lowFreqVarianceRemoved', zeros(EEG.nbchan,1));

    if isempty(which('nt_detrend'))
        info.status = 'nt_detrend not on the path - call RELAX_add_vendored_noisetools()';
        return
    end

    segments   = RELAX_get_segments(EEG);
    wsize      = 2 * floor(round(opt.MaxWindowSeconds * EEG.srate) / 2);  % samples, even
    origClass  = class(EEG.data);
    data       = double(EEG.data);
    repaired   = data;

    allDurations  = [];
    meanOnlySamps = 0;

    for ch = 1:EEG.nbchan
        % --- build this channel's partition ---
        chSteps = steps.sample(steps.channel == ch);
        cuts    = segments(:,1);
        for k = 1:numel(chSteps)
            % merge a step that sits on top of a segment edge rather than slivering it
            if all(abs(chSteps(k) - segments(:,1)) > opt.MergeSamples) && ...
               all(abs(chSteps(k) - (segments(:,2)+1)) > opt.MergeSamples)
                cuts(end+1) = chSteps(k); %#ok<AGROW>
            end
        end
        % cuts(:) guards against MATLAB's ambiguous-orientation growth: when segments
        % has exactly one row (the common no-boundary case), segments(:,1) is a 1x1
        % array, and cuts(end+1)=... on a 1x1 silently grows it into a ROW vector
        % instead of a column, which then breaks the column-vector vertcat below.
        cuts = unique(cuts(:));

        subStarts = cuts;
        subEnds   = [cuts(2:end) - 1; EEG.pnts];
        % a cut may fall inside a later boundary segment; clip each sub-segment to its own
        for k = 1:numel(subStarts)
            s = find(segments(:,1) <= subStarts(k) & segments(:,2) >= subStarts(k), 1);
            subEnds(k) = min(subEnds(k), segments(s,2));
        end
        keep = subEnds >= subStarts;
        subStarts = subStarts(keep);
        subEnds   = subEnds(keep);

        info.subSegmentsPerChannel(ch) = numel(subStarts);

        for k = 1:numel(subStarts)
            a = subStarts(k); b = subEnds(k);
            x = data(ch, a:b)';                 % nt_detrend wants time on dim 1
            T = numel(x) / EEG.srate;
            allDurations(end+1,1) = T;          %#ok<AGROW>

            order = chooseOrder(T, C, opt.MaxPolyOrder, opt.MaxCutoffHz);

            if any(isnan(x)) || (max(x) - min(x)) == 0
                % A NaN would propagate through the whole fit; a constant stretch has no
                % trend to find and a zero scale for the reweighting.
                order = 0;
            end

            % Compare against the ACHIEVABLE ceiling, not the raw configured value.
            % chooseOrder saturates at numel(C), so gating on opt.MaxPolyOrder alone
            % makes this condition unsatisfiable whenever MaxPolyOrder exceeds numel(C),
            % silently switching overlap-add off for every long sub-segment.
            useOverlapAdd = (order == min(opt.MaxPolyOrder, numel(C))) ...
                            && (numel(x) > wsize) && (wsize > 0);

            done = false;
            while ~done
                if order == 0
                    y = robustDemean(x, opt.DetrendThresh, opt.DetrendIterations);
                    counts.meanOnly = counts.meanOnly + 1;
                    meanOnlySamps = meanOnlySamps + numel(x);
                    done = true;
                else
                    try
                        if useOverlapAdd
                            y = nt_detrend(x, order, [], [], opt.DetrendThresh, ...
                                           opt.DetrendIterations, wsize);
                            % Name the order actually used. order is 1..numel(C) here:
                            % the order==0 case is handled by the meanOnly branch above,
                            % and chooseOrder cannot exceed numel(C).
                            oaField = sprintf('order%doa', order);
                            counts.(oaField) = counts.(oaField) + 1;
                        else
                            y = nt_detrend(x, order, [], [], opt.DetrendThresh, ...
                                           opt.DetrendIterations);
                            switch order
                                case 3, counts.order3 = counts.order3 + 1;
                                case 2, counts.order2 = counts.order2 + 1;
                                case 1, counts.order1 = counts.order1 + 1;
                            end
                        end
                        if any(isnan(y)); error('RELAX:nanFromDetrend', 'nt_detrend returned NaN'); end
                        done = true;
                    catch
                        % nt_detrend's degenerate-window failure is a bare error('!') with
                        % no identifier. Step DOWN the ladder rather than killing a
                        % 257-channel SLURM task hours into a run.
                        counts.fallback = counts.fallback + 1;
                        useOverlapAdd = false;
                        order = order - 1;
                        if order < 1; order = 0; end
                    end
                end
            end

            repaired(ch, a:b) = y';
        end
    end

    removed = data - repaired;
    for ch = 1:EEG.nbchan
        v = var(data(ch, ~isnan(data(ch,:))));
        if v > 0
            info.lowFreqVarianceRemoved(ch) = var(removed(ch, ~isnan(removed(ch,:)))) / v;
        end
    end

    EEG.data = cast(repaired, origClass);

    info.available               = true;
    info.status                  = 'ok';
    info.nChannelsRepaired       = EEG.nbchan;
    info.nSubSegments            = sum(info.subSegmentsPerChannel);
    info.treatmentCounts         = counts;
    info.minSubSegmentDuration   = min(allDurations);
    info.medianSubSegmentDuration = median(allDurations);
    info.fractionSamplesMeanOnly = meanOnlySamps / (EEG.nbchan * EEG.pnts);
end

function order = chooseOrder(T, C, maxOrder, maxCutoffHz)
% Largest order n <= maxOrder whose measured -3 dB point stays at or below maxCutoffHz.
% 0 means "no polynomial fits in this duration" - robust mean subtraction only.
    order = 0;
    % C(n)*n/(2*T) and maxCutoffHz are independently-rounded floating-point values that
    % need not cancel exactly even when the duration T is chosen to land exactly on a
    % documented rung boundary: at T = 1.635 s, n = 3, C(3) = 1.09, the comparison
    % evaluates to 1.0 + 1 ULP and fails a strict <=, silently dropping the sub-segment
    % to order 2 even though the docstring table documents that boundary as order 3.
    % A tolerance of a few ULPs of maxCutoffHz rescues exact-boundary cases without
    % shifting any other rung (verified: this tolerance does not change the chosen order
    % at any of the non-boundary durations exercised by the test suite). Do not replace
    % this with a hardcoded absolute constant (e.g. 1e-9) - it must stay correct if
    % MaxCutoffHz is configured far from 1.0.
    tol = 4 * eps(maxCutoffHz);
    for n = 1:min(maxOrder, numel(C))
        if C(n) * n / (2*T) <= maxCutoffHz + tol
            order = n;
        end
    end
end

function y = robustDemean(x, thresh, niter)
% Iteratively reweighted mean subtraction: the order-0 rung of the ladder, and the
% fallback when nt_detrend cannot fit. Written out rather than calling nt_detrend(x,0,...)
% because order 0 builds an empty polynomial basis, which nt_regw cannot regress on.
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
        if ~any(wNew); break; end   % never leave an all-zero weight vector behind
        w = wNew;
    end
    m = sum(w(finite) .* x(finite)) / sum(w(finite));
    y = x - m;
end
