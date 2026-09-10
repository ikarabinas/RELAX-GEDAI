function chunkIndices = zapline_chunk_indices(nSamples, epochs, epochSegment, segments, peakLocs, peakHeights, minChunkSamples, edgeMarginSamples)
% ZAPLINE_CHUNK_INDICES - chunk edges that never fragment a stretch of continuous data.
%
%   nSamples          - number of samples in the record
%   epochs, epochSegment - from zapline_covariance_epochs
%   segments          - from zapline_boundary_segments
%   peakLocs          - findpeaks locations in the DISTANCE vector: peakLocs(i)==k means
%                       the covariance jumped between epoch k and epoch k+1
%   peakHeights       - the corresponding peak values, used to resolve spacing conflicts
%   minChunkSamples   - minimum chunk length in samples (minChunkLength*srate)
%   edgeMarginSamples - an edge inside a segment must leave at least this many samples
%                       of that segment on both sides
%
%   chunkIndices      - column vector, [1; edges; nSamples+1], consumed as
%                       data(chunkIndices(i):chunkIndices(i+1)-1,:)
%
% RELAX-GEDAI local addition to zapline-plus. See vendor/zapline-plus/VENDORING.md.

    candidates = zeros(0, 1);
    heights    = zeros(0, 1);

    for i = 1:numel(peakLocs)
        nextEpoch = peakLocs(i) + 1;
        if nextEpoch > size(epochs, 1)
            continue;
        end

        % The edge is the FIRST sample of the following epoch. Upstream used the last
        % sample of the current one, which under boundary alignment lands one sample
        % before a boundary and leaves a 1-sample fragment behind.
        edge = epochs(nextEpoch, 1);
        seg  = segments(epochSegment(nextEpoch), :);

        isSegmentStart = (edge == seg(1));
        leavesRoom     = (edge - seg(1) >= edgeMarginSamples) && ...
                         (seg(2) - edge + 1 >= edgeMarginSamples);

        if isSegmentStart || leavesRoom
            candidates(end+1, 1) = edge;             %#ok<AGROW>
            heights(end+1, 1)    = peakHeights(i);   %#ok<AGROW>
        end
    end

    % Reproduce findpeaks' MinPeakDistance rule in the SAMPLE domain: keep the tallest
    % peak, then drop every candidate within minChunkSamples of one already kept.
    % findpeaks itself cannot do this any more, because its distance is measured in
    % epoch indices and boundary-aligned epochs are not all the same length.
    [~, order] = sort(heights, 'descend');
    accepted = zeros(0, 1);
    for i = reshape(order, 1, [])
        if isempty(accepted) || all(abs(accepted - candidates(i)) >= minChunkSamples)
            accepted(end+1, 1) = candidates(i);      %#ok<AGROW>
        end
    end
    accepted = unique(accepted);

    chunkIndices = [1; accepted(:); nSamples + 1];

    % Upstream's guards: the first and last chunks must also reach the minimum length.
    if numel(chunkIndices) >= 3 && chunkIndices(2) - chunkIndices(1) < minChunkSamples
        chunkIndices(2) = [];
    end
    if numel(chunkIndices) >= 3 && chunkIndices(end) - chunkIndices(end-1) < minChunkSamples
        chunkIndices(end-1) = [];
    end
end
