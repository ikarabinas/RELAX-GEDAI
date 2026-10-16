% RELAX_QC_blink_peaks() - Samples of the pipeline's own EyeBlinkMax markers that can be
%                          epoched without spanning a boundary.
%
% Usage:
%   >> peaks = RELAX_QC_blink_peaks(EEG, halfWindowSamples);
%
% Outputs:
%   peaks - [1 x K] sorted sample indices whose [peak - halfWindowSamples,
%           peak + halfWindowSamples] lies wholly inside one RELAX_get_segments segment

function peaks = RELAX_QC_blink_peaks(EEG, halfWindowSamples)

    peaks = zeros(1, 0);
    if ~isfield(EEG, 'event') || isempty(EEG.event), return, end

    isBlink = arrayfun(@(e) (ischar(e.type) || isstring(e.type)) && strcmp(e.type, 'EyeBlinkMax'), EEG.event);
    segments = RELAX_get_segments(EEG);
    for p = round([EEG.event(isBlink).latency])
        s = find(segments(:, 1) <= p & segments(:, 2) >= p, 1);
        if ~isempty(s) && p - halfWindowSamples >= segments(s, 1) && p + halfWindowSamples <= segments(s, 2)
            peaks(end+1) = p; %#ok<AGROW>
        end
    end
    peaks = sort(peaks);
end
