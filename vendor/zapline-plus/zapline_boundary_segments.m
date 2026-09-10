function segments = zapline_boundary_segments(nSamples, boundaries)
% ZAPLINE_BOUNDARY_SEGMENTS - stretches of continuous data between discontinuities.
%
%   nSamples   - number of samples in the record
%   boundaries - sample indices of discontinuities. Each value is the FIRST sample
%                of a new segment: the eeglab wrapper writes round(EEG.event.latency),
%                and EEGLAB places boundaries at x.5, so round() lands on the first
%                post-gap sample. Same convention as pwelch_segmented.
%
%   segments   - [nSegments x 2] inclusive 1-based [start end] rows. Contiguous and
%                exhaustive: segments(1,1) == 1, segments(end,2) == nSamples, and each
%                row starts one sample after the previous row ends. Empty or unusable
%                boundaries degrade to a single whole-record segment, so continuous
%                data behaves exactly as it did before. A boundary at exactly nSamples
%                produces a legal 1-sample trailing segment; dropping such slivers is
%                deliberately left to the epoch layer (zapline_covariance_epochs drops
%                any segment below minEpochSamples), not this function's job.
%
% RELAX-GEDAI local addition to zapline-plus. See vendor/zapline-plus/VENDORING.md.

    if isempty(boundaries)
        segments = [1 nSamples];
        return;
    end

    b = unique(round(boundaries(:)'));
    b = b(b > 1 & b <= nSamples);

    if isempty(b)
        segments = [1 nSamples];
        return;
    end

    starts = [1, b];
    ends   = [b - 1, nSamples];
    segments = [starts(:), ends(:)];
end
