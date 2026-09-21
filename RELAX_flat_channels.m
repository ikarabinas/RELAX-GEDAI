% RELAX_flat_channels() - Which channels carry nothing: constant (or non-finite) inside
%                         EVERY segment.
%
% Usage:
%   >> flat = RELAX_flat_channels(EEG.data, RELAX_get_segments(EEG));
%
% A globally constant channel, a railed one that sits at a different constant level in
% each inter-boundary segment (m544 day 2, E154: 0 / +333333 uV), and an all-NaN channel
% are all flat. One segment with any variation keeps the channel.
%
% Inputs:
%   data     - nbchan x pnts
%   segments - [nSeg x 2] 1-based inclusive sample ranges (RELAX_get_segments)
%
% Outputs:
%   flat - nbchan x 1 logical

function flat = RELAX_flat_channels(data, segments)
    nb = size(data, 1);
    flat = true(nb, 1);
    for ch = 1:nb
        for s = 1:size(segments, 1)
            x = double(data(ch, segments(s,1):segments(s,2)));
            x = x(isfinite(x));
            if ~isempty(x) && max(x) - min(x) > 0
                flat(ch) = false;
                break
            end
        end
    end
end
