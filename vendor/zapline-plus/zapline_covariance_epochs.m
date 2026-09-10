function [epochs, epochSegment] = zapline_covariance_epochs(segments, epochSamples, wholeSegmentSamples, minEpochSamples)
% ZAPLINE_COVARIANCE_EPOCHS - epochs for adaptive chunk detection that respect boundaries.
%
% Upstream zapline tiles fixed-length covariance epochs flatly from sample 1, so an
% epoch can straddle a discontinuity and a chunk edge derived from it can slice a
% stretch of continuous data into a useless fragment. This tiles inside each segment
% instead.
%
%   segments            - [nSegments x 2] from zapline_boundary_segments
%   epochSamples        - nominal epoch length in samples (segmentLength*srate)
%   wholeSegmentSamples - segments shorter than this become ONE epoch spanning the
%                         whole segment, so they can never be subdivided
%   minEpochSamples     - segments shorter than this emit NO epoch. A covariance from
%                         too few samples is rank-starved and reads as a large distance
%                         from its neighbours, manufacturing a spurious chunk edge.
%
%   epochs              - [nEpochs x 2] inclusive [start end]. Every epoch lies wholly
%                         inside one segment and is at least epochSamples long.
%   epochSegment        - [nEpochs x 1] row index into segments for each epoch.
%
% RELAX-GEDAI local addition to zapline-plus. See vendor/zapline-plus/VENDORING.md.

    % Both thresholds are floored at epochSamples so "every epoch is at least
    % segmentLength long" holds however the caller parameterises this.
    minEpochSamples     = max(minEpochSamples, epochSamples);
    wholeSegmentSamples = max(wholeSegmentSamples, epochSamples);

    epochs       = zeros(0, 2);
    epochSegment = zeros(0, 1);

    for k = 1:size(segments, 1)
        segStart = segments(k, 1);
        segEnd   = segments(k, 2);
        segLen   = segEnd - segStart + 1;

        if segLen < minEpochSamples
            continue;
        end

        if segLen < wholeSegmentSamples
            epochs(end+1, :)       = [segStart segEnd];       %#ok<AGROW>
            epochSegment(end+1, 1) = k;                       %#ok<AGROW>
            continue;
        end

        nEpochs = floor(segLen / epochSamples);
        for j = 1:nEpochs
            thisStart = segStart + (j-1)*epochSamples;
            if j < nEpochs
                thisEnd = segStart + j*epochSamples - 1;
            else
                % The last epoch absorbs the remainder so it ends exactly on the
                % boundary and is never shorter than epochSamples.
                thisEnd = segEnd;
            end
            epochs(end+1, :)       = [thisStart thisEnd];     %#ok<AGROW>
            epochSegment(end+1, 1) = k;                       %#ok<AGROW>
        end
    end
end
