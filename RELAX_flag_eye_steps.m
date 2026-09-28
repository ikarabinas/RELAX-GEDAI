% RELAX_flag_eye_steps() - Mark detected steps that are eye movements rather than pops.
%
% Usage:
%   >> [steps, info] = RELAX_flag_eye_steps(steps, eyeMask, blinkLabels);
%   >> [steps, info] = RELAX_flag_eye_steps(steps, eyeMask, blinkLabels, 'MaxAmplitudeMicrovolts', 500);
%
% A horizontal saccade produces a 20-60 ms, 50-300 uV same-sign level shift on a handful of
% frontal / peri-ocular channels, which the two-window step detector cannot tell from a
% local cluster of electrode pops. Repairing it would remove EOG structure the later ICA
% relies on to isolate eye components. This function flags, for the repair stage to skip,
% every member of a co-occurrence event (steps.eventId) that
%   (1) spans at least two channels,
%   (2) includes at least two channels listed in blinkLabels (the peri-ocular electrodes,
%       RELAX_cfg.BlinkElectrodes plus the HEOG pair) - or at least one when the event is
%       common-mode (steps.commonMode), where an eye electrode is one of many,
%   (3) has the mid-level crossing (steps.sample) of at least one of those eye electrodes
%       inside an eye-movement period in eyeMask, and
%   (4) has no member whose |amplitude| exceeds MaxAmplitudeMicrovolts.
% A single-channel step on a blink electrode is never flagged: blinks and saccades are
% shared across neighbours, so a lone step on a face electrode is a pop. Common-mode
% events are eligible too - RELAX_repair_steps then leaves them in place rather than
% excising them, since a rim-wide shift inside an eye period is EOG, not head movement.
%
% Rules (2)-(4) were tightened on 2026-09-18 after an audit of the five validation files:
% with the blink mask covering 25-30 % of a file, 2-5-channel events containing one eye
% electrode were being flagged at 1.2-1.7x the chance-overlap rate, i.e. a pop pair with
% one member on E18 was usually left unrepaired; and on m155 every eye-electrode step (all
% 1.5-2.6 mV, 3-channel clusters) was flagged. Saccades and blinks stay well under 500 uV
% on this montage; a transition interval [onset, offset] is wide enough to touch a +-400 ms
% blink window by chance far more often than its mid-crossing is.
%
% Inputs:
%   steps       - struct from RELAX_detect_steps (needs .channel .label .sample .amplitude
%                 .eventId .commonMode)
%   eyeMask     - logical [1 x EEG.pnts] at the SAME sample rate and timeline as steps, true
%                 during eye-movement periods (e.g. the blink mask produced by
%                 RELAX_blinks_IQR_method on a temporary copy, OR'd with the horizontal
%                 eye-movement mask when available)
%   blinkLabels - cell array of channel labels counted as eye electrodes
%   Name-value pairs (optional):
%     MaxAmplitudeMicrovolts - an event with any member over this |amplitude| is a pop
%                              cluster, not eye movement (default 500)
%
% Outputs:
%   steps - same struct with a logical [N x 1] .eyeMovement column added
%   info  - .nEyeSteps, .nEyeEvents, .nEventsOverAmplitude (eligible events that failed
%           rule 4 alone)

function [steps, info] = RELAX_flag_eye_steps(steps, eyeMask, blinkLabels, varargin)

    p = inputParser;
    p.addParameter('MaxAmplitudeMicrovolts', 500, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.parse(varargin{:});
    maxAmplitude = p.Results.MaxAmplitudeMicrovolts;

    n = numel(steps.sample);
    steps.eyeMovement = false(n, 1);
    info = struct('nEyeSteps', 0, 'nEyeEvents', 0, 'nEventsOverAmplitude', 0);
    if n == 0; return; end

    eyeMask = logical(eyeMask(:)');
    isBlinkCh = ismember(steps.label(:), blinkLabels(:));
    sampleInMask = false(n, 1);
    inRange = steps.sample(:) >= 1 & steps.sample(:) <= numel(eyeMask);
    sampleInMask(inRange) = eyeMask(steps.sample(inRange));
    ids = unique(steps.eventId(:))';
    nEv = 0;
    for e = ids
        m = steps.eventId == e;
        if numel(unique(steps.channel(m))) < 2; continue; end
        nEyeChannels = numel(unique(steps.channel(m & isBlinkCh)));
        if any(steps.commonMode(m)); needed = 1; else; needed = 2; end
        if nEyeChannels < needed; continue; end
        if ~any(sampleInMask(m & isBlinkCh)); continue; end
        if any(abs(steps.amplitude(m)) > maxAmplitude)
            info.nEventsOverAmplitude = info.nEventsOverAmplitude + 1;
            continue;
        end
        steps.eyeMovement(m) = true;
        nEv = nEv + 1;
    end
    info.nEyeSteps  = sum(steps.eyeMovement);
    info.nEyeEvents = nEv;
end
