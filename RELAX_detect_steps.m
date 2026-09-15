% RELAX_detect_steps() - Detect electrode-pop step artifacts in continuous EEG.
%
% Usage:
%   >> [steps, info] = RELAX_detect_steps(EEG);
%   >> [steps, info] = RELAX_detect_steps(EEG, 'ZThreshold', 20, 'RefractoryMs', 20);
%
% An electrode pop is a near-instantaneous change in a single channel's DC baseline,
% caused by a momentary change in the electrode-electrolyte contact. It is a STEP, not a
% transient: the baseline moves and stays moved. Detection follows FieldTrip's
% ft_artifact_jump recipe (median filter, absolute derivative, z threshold), with two
% additions this data needed: an absolute microvolt floor, and a rise-time criterion
% without which smooth ramps are flagged as steps.
%
% Run this BEFORE any filtering. A zero-phase highpass converts a step into a bipolar
% excursion spread over seconds in BOTH time directions, after which no step remains to
% find.
%
% Inputs:
%   EEG - EEGLAB EEG structure (uses .data, .nbchan, .pnts, .srate, .chanlocs, .event)
%   Name-value pairs (all optional):
%     MedianFilterOrder    - median filter length in samples, suppresses single-sample
%                             spikes without smearing a genuine step edge (default 9,
%                             FieldTrip's value)
%     ZThreshold           - robust z threshold on the absolute derivative (default 20,
%                             FieldTrip's value)
%     MinAmplitudeMicrovolts - absolute floor on the level change (default 10). Stops a
%                             very quiet channel generating detections from its own noise
%                             floor, where a tiny step is enormous in z.
%     RefractoryMs         - minimum spacing between steps on one channel (default 20)
%     SustainMs            - window either side used to measure the pre- and post-step
%                             level (default 100). The level change must persist across
%                             it - that is what makes this a step detector rather than a
%                             spike detector.
%     MinRiseFraction      - fraction of the total level change that must occur in the
%                             single largest derivative sample (default 0.5). This is the
%                             rise-time criterion: a genuine pop delivers most of its
%                             change in one sample, a ramp spreads it over thousands.
%     StepsPerMinuteToRejectChannel - per-channel step rate above which the channel is a
%                             chronically bad contact rather than a repairable one
%                             (default 30). Derived from survey data: the clearest bad
%                             contact ran at ~72 steps/min while sporadically popping
%                             channels averaged well under 1/min, so 30 sits in the gap.
%     CommonModeBinMs      - bin width for the reference-pop co-occurrence count (default 25)
%     CommonModeChannelFraction - fraction of channels in one bin above which the bin is
%                             counted as a common-mode (reference) step (default 0.25).
%                             REPORTING ONLY - such steps are still repaired per channel.
%
% Outputs:
%   steps - struct of parallel [N x 1] columns, one entry per step, sorted by channel then
%           sample: .channel (row index into EEG.data), .label (channel label, so the
%           struct survives a pop_select), .sample (FIRST sample of the new level),
%           .amplitude (signed uV, post-level minus pre-level), .z. An empty result has
%           0x1 / {} fields, never a bare [], so numel(steps.sample) is always valid.
%   info  - .available, .status, .params, .nSteps, .nChannelsAffected, .stepsPerChannel,
%           .stepRatePerMinute, .badChannels, .badChannelReason, .degenerateChannels,
%           .nCommonModeSteps, .commonModeSamples
%           .badChannels/.badChannelReason are ordered worst-first (highest step rate
%           first), matching the ordering RELAX_channel_budget's 'spend' requires of any
%           list it truncates.
%
% Detection runs inside each inter-boundary segment separately (via RELAX_get_segments),
% so the splice discontinuity at a boundary is never mistaken for a pop. The robust
% z-scale, however, is computed per channel across the WHOLE recording: an 8 s segment is
% too short for a stable scale estimate, and MAD is insensitive to the pops themselves.

function [steps, info] = RELAX_detect_steps(EEG, varargin)

    p = inputParser;
    p.addParameter('MedianFilterOrder', 9);
    p.addParameter('ZThreshold', 20);
    p.addParameter('MinAmplitudeMicrovolts', 10);
    p.addParameter('RefractoryMs', 20);
    p.addParameter('SustainMs', 100);
    p.addParameter('MinRiseFraction', 0.5);
    p.addParameter('StepsPerMinuteToRejectChannel', 30);
    p.addParameter('CommonModeBinMs', 25);
    p.addParameter('CommonModeChannelFraction', 0.25);
    p.parse(varargin{:});
    opt = p.Results;

    steps = emptySteps();
    info  = struct('available', false, 'status', '', 'params', opt, ...
                   'nSteps', 0, 'nChannelsAffected', 0, ...
                   'stepsPerChannel', zeros(EEG.nbchan,1), ...
                   'stepRatePerMinute', zeros(EEG.nbchan,1), ...
                   'badChannels', {{}}, 'badChannelReason', {{}}, ...
                   'degenerateChannels', {{}}, ...
                   'nCommonModeSteps', 0, 'commonModeSamples', []);

    refractory = max(1, round(opt.RefractoryMs * EEG.srate / 1000));
    sustain    = max(2, round(opt.SustainMs    * EEG.srate / 1000));

    segments = RELAX_get_segments(EEG);
    % A segment must hold a pre-window, a post-window and the step itself to be usable.
    minSegLen = 2*sustain + 2;
    usable = (segments(:,2) - segments(:,1) + 1) >= minSegLen;
    if ~any(usable)
        info.status = sprintf(['no inter-boundary segment reaches the %d samples needed for a ' ...
            '%g ms sustain window either side'], minSegLen, opt.SustainMs);
        return
    end
    segments = segments(usable, :);

    labels = {EEG.chanlocs.labels};
    durationMinutes = (EEG.pnts / EEG.srate) / 60;

    allCh = []; allSamp = []; allAmp = []; allZ = [];
    degenerate = {};

    for ch = 1:EEG.nbchan
        x = double(EEG.data(ch, :));

        % 1. median-filter and differentiate, inside each segment only
        d = nan(1, EEG.pnts - 1);   % d(i) is the transition between sample i and i+1
        for s = 1:size(segments,1)
            a = segments(s,1); b = segments(s,2);
            xm = movmedian(x(a:b), opt.MedianFilterOrder);
            d(a : b-1) = abs(diff(xm));
        end
        dValid = d(~isnan(d));

        % 2. robust scale across the whole recording
        med = median(dValid);
        mad = median(abs(dValid - med));
        if mad <= 0
            % A flat, dead or heavily quantised channel. No scale exists, so no z exists.
            degenerate{end+1} = labels{ch}; %#ok<AGROW>
            continue
        end
        z = (d - med) / (1.4826 * mad);

        % 3. candidates over threshold, strongest first, refractory enforced greedily
        cand = find(z > opt.ZThreshold);
        if isempty(cand); continue; end
        [~, order] = sort(z(cand), 'descend');
        cand = cand(order);

        taken = false(1, EEG.pnts);
        for k = 1:numel(cand)
            i = cand(k);                       % transition index
            lo = max(1, i - refractory);
            hi = min(EEG.pnts, i + refractory);
            if any(taken(lo:hi)); continue; end

            stepSample = i + 1;                % first sample of the NEW level

            % the sustain windows must fit inside THIS segment
            s = find(segments(:,1) <= i & segments(:,2) >= stepSample, 1);
            if isempty(s); continue; end
            a = segments(s,1); b = segments(s,2);
            if (i - a + 1) < sustain || (b - stepSample + 1) < sustain; continue; end

            preLevel  = median(x(i - sustain + 1 : i));
            postLevel = median(x(stepSample : stepSample + sustain - 1));
            amplitude = postLevel - preLevel;

            % 4. absolute floor
            if abs(amplitude) < opt.MinAmplitudeMicrovolts; continue; end

            % 5. rise-time criterion: most of the change in one sample.
            %    A ramp delivers its change over thousands of samples and fails here.
            if abs(x(stepSample) - x(i)) < opt.MinRiseFraction * abs(amplitude); continue; end

            taken(lo:hi) = true;
            allCh(end+1,1)   = ch;          %#ok<AGROW>
            allSamp(end+1,1) = stepSample;  %#ok<AGROW>
            allAmp(end+1,1)  = amplitude;   %#ok<AGROW>
            allZ(end+1,1)    = z(i);        %#ok<AGROW>
        end
    end

    info.degenerateChannels = degenerate;
    info.available = true;
    info.status    = 'ok';

    if isempty(allCh)
        return
    end

    % sort by channel then sample
    [~, ord] = sortrows([allCh allSamp]);
    steps.channel   = allCh(ord);
    steps.sample    = allSamp(ord);
    steps.amplitude = allAmp(ord);
    steps.z         = allZ(ord);
    steps.label     = reshape(labels(steps.channel), [], 1);

    info.nSteps            = numel(steps.sample);
    info.stepsPerChannel   = accumarray(steps.channel, 1, [EEG.nbchan 1]);
    info.stepRatePerMinute = info.stepsPerChannel / max(durationMinutes, eps);
    info.nChannelsAffected = sum(info.stepsPerChannel > 0);

    % Population A: chronically bad contacts, ordered worst-first for the shared budget
    badIdx = find(info.stepRatePerMinute > opt.StepsPerMinuteToRejectChannel);
    [~, worstFirst] = sort(info.stepRatePerMinute(badIdx), 'descend');
    badIdx = badIdx(worstFirst);
    info.badChannels = reshape(labels(badIdx), [], 1);
    info.badChannelReason = reshape(arrayfun(@(k) sprintf('%.1f steps/min exceeds %g', ...
        info.stepRatePerMinute(k), opt.StepsPerMinuteToRejectChannel), ...
        badIdx, 'UniformOutput', false), [], 1);

    % Common-mode (reference) steps: reporting only, per design decision 10
    binSamples = max(1, round(opt.CommonModeBinMs * EEG.srate / 1000));
    binIdx = floor((steps.sample - 1) / binSamples) + 1;
    uniqueBins = unique(binIdx);
    nPerBin = arrayfun(@(b) numel(unique(steps.channel(binIdx == b))), uniqueBins);
    commonBins = uniqueBins(nPerBin > opt.CommonModeChannelFraction * EEG.nbchan);
    info.nCommonModeSteps  = numel(commonBins);
    info.commonModeSamples = (commonBins - 1) * binSamples + 1;
end

function s = emptySteps()
    s = struct('channel', zeros(0,1), 'label', {cell(0,1)}, 'sample', zeros(0,1), ...
               'amplitude', zeros(0,1), 'z', zeros(0,1));
end
