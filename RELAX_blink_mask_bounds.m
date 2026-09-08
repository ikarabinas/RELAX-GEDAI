% RELAX_blink_mask_bounds() - Blink mask window, truncated to the segment containing
%                             the blink.
%
% Usage:
%   >> [maskStart, maskEnd] = RELAX_blink_mask_bounds(latency, halfWidth, segments, nPoints);
%
% Inputs:
%   latency   - sample index of the blink peak (may be non-integer; it is rounded)
%   halfWidth - how many samples to mark either side of the peak
%   segments  - [nSegments x 2] matrix from RELAX_get_segments
%   nPoints   - total samples in the record, used as the fallback range
%
% Outputs:
%   maskStart, maskEnd - inclusive integer bounds of the blink window, truncated to
%                        the segment containing `latency`. maskStart <= maskEnd always.
%
% Every blink gets marked; a blink close to a discontinuity simply gets a shorter
% window rather than being discarded. Truncating rather than dropping matters here
% because the segments between boundary events are only ~8 s long, so a large share of
% blinks sit within 400 ms of a segment edge.

function [maskStart, maskEnd] = RELAX_blink_mask_bounds(latency, halfWidth, segments, nPoints)

    segStart = 1;
    segEnd   = nPoints;

    match = find(segments(:,1) <= latency & segments(:,2) >= latency, 1, 'first');
    if ~isempty(match)
        segStart = segments(match,1);
        segEnd   = segments(match,2);
    end

    % Clamp the peak into its own segment first. A real blink peak is always inside
    % its segment already, so this is a no-op in practice; doing it up front means
    % maskStart <= maskEnd holds structurally even in the fallback case of a latency
    % that matched no segment.
    peak = min(max(round(latency), segStart), segEnd);

    maskStart = max(peak - round(halfWidth), segStart);
    maskEnd   = min(peak + round(halfWidth), segEnd);
end
