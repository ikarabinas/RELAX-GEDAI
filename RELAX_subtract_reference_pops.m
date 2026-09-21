% RELAX_subtract_reference_pops() - Remove detected reference-electrode pops from every
%                                   channel.
%
% Usage:
%   >> [EEG, info] = RELAX_subtract_reference_pops(EEG, refSteps, refTrace);
%   >> [EEG, info] = RELAX_subtract_reference_pops(EEG, refSteps, refTrace, 'WindowMs', 150);
%
% Every channel is x_i(t) = s_i(t) - r(t). For each reference pop [onset, offset] found on
% the cross-channel median refTrace (RELAX_detect_reference_pops), the correction is
%     0                         for t <= onset
%     refTrace(t) - Lpre        for onset < t < offset     (the transition itself)
%     amplitude                 for offset <= t <= segment end
% with Lpre the median of refTrace over the WindowMs samples ending at onset, and it is
% subtracted from every channel that is not flat (RELAX_flat_channels). Several pops in one
% segment stack because each is relative to its own pre-level. Boundaries are respected:
% a correction never crosses into the next inter-boundary segment. Apply the same call to
% the auxiliary (keep-but-don't-clean) electrodes: they share the reference.
%
% Inputs:
%   EEG      - EEGLAB EEG structure (.data, .nbchan, .pnts, .srate, .chanlocs, .event)
%   refSteps - struct with N x 1 .onset .offset .amplitude (RELAX_detect_reference_pops)
%   refTrace - 1 x pnts virtual channel (refInfo.refChannel.data)
%   Name-value pairs (optional):
%     WindowMs - pre-level window (default 150; use the detector's WindowMs)
%
% Outputs:
%   EEG  - .data corrected (class preserved); nothing else changed
%   info - .available .status .nSteps .nChannelsCorrected .skippedLabels (flat channels)
%          .transitionSeconds .medianAmplitude (of |amplitude|; NaN when empty)
%          .correction (1 x pnts single: what was subtracted, for auditing)

function [EEG, info] = RELAX_subtract_reference_pops(EEG, refSteps, refTrace, varargin)
    p = inputParser;
    p.addParameter('WindowMs', 150);
    p.parse(varargin{:});
    W = max(1, round(p.Results.WindowMs * EEG.srate / 1000));

    n = EEG.pnts;
    labels = reshape({EEG.chanlocs.labels}, [], 1);
    segments = RELAX_get_segments(EEG);
    flat = RELAX_flat_channels(EEG.data, segments);
    nSteps = numel(refSteps.onset);

    info = struct('available', false, 'status', '', 'nSteps', nSteps, ...
                  'nChannelsCorrected', sum(~flat), 'skippedLabels', {labels(flat)}, ...
                  'transitionSeconds', 0, 'medianAmplitude', NaN, 'correction', zeros(1, n, 'single'));
    if numel(refTrace) ~= n
        info.status = sprintf('refTrace has %d samples, EEG has %d', numel(refTrace), n);
        return
    end
    refTrace = double(refTrace(:))';

    correction = zeros(1, n);
    for k = 1:nSteps
        onset = refSteps.onset(k); offset = refSteps.offset(k);
        seg = find(segments(:,1) <= onset & segments(:,2) >= offset, 1);
        if isempty(seg); continue; end
        a = segments(seg,1); b = segments(seg,2);
        Lpre = median(refTrace(max(a, onset - W + 1) : onset), 'omitnan');
        inside = onset + 1 : offset - 1;
        correction(inside) = correction(inside) + (refTrace(inside) - Lpre);
        correction(offset : b) = correction(offset : b) + refSteps.amplitude(k);
        info.transitionSeconds = info.transitionSeconds + numel(inside) / EEG.srate;
    end

    apply = find(~flat);
    EEG.data(apply, :) = EEG.data(apply, :) - cast(correction, 'like', EEG.data);
    info.correction = single(correction);
    if nSteps > 0; info.medianAmplitude = median(abs(refSteps.amplitude)); end
    info.available = true;
    info.status = 'ok';
end