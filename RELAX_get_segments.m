% RELAX_get_segments() - Return the sample ranges of contiguous data between
%                        boundary events.
%
% Usage:
%   >> segments = RELAX_get_segments(EEG);
%
% Inputs:
%   EEG - EEGLAB EEG structure. Requires .pnts; uses .event if present.
%
% Outputs:
%   segments - [nSegments x 2] matrix of [startSample endSample], 1-based and
%              inclusive, one row per stretch of contiguous data. Data with no
%              boundary events returns a single row [1 EEG.pnts], so callers
%              degrade to whole-record behaviour on ordinary continuous data.
%
% Boundary latencies are converted with the integer-safe form
%   starts = [1, floor(b)+1],  ends = [floor(b), pnts]
% which handles both EEGLAB's convention of placing boundaries at x.5 (between
% samples) and tools that write integer boundary latencies, without ever
% assigning the same sample to two segments.

function segments = RELAX_get_segments(EEG)

    nPoints = EEG.pnts;

    boundaries = [];
    if isfield(EEG,'event') && ~isempty(EEG.event)
        isBoundary = false(1,numel(EEG.event));
        for e = 1:numel(EEG.event)
            thisType = EEG.event(e).type;
            if (ischar(thisType) || isstring(thisType)) && strcmp(thisType,'boundary')
                isBoundary(e) = true;
            end
        end
        boundaries = sort([EEG.event(isBoundary).latency]);
    end

    % Ignore boundaries outside the data, and boundaries at the very first or last
    % sample, which do not split the record into two usable pieces:
    boundaries = boundaries(boundaries > 1 & boundaries < nPoints);

    starts = [1, floor(boundaries) + 1];
    ends   = [floor(boundaries), nPoints];

    valid = starts <= ends;
    segments = [starts(valid)', ends(valid)'];
end
