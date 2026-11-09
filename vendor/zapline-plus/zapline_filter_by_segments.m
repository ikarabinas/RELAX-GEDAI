function out = zapline_filter_by_segments(data, segments, fn, padWidth, minSegmentSamples)
% ZAPLINE_FILTER_BY_SEGMENTS - apply a filter within each segment, never across one.
%
%   data              - [nSamples x nChannels], zapline's own orientation. The full
%                       array is never transposed: it can be 257 x 245000 doubles.
%   segments          - [nSegments x 2] from zapline_boundary_segments
%   fn                - function handle taking and returning [nSamples x nChannels]
%   padWidth          - reflection pad in samples applied to each end of a segment,
%                       clamped to segLen-1. One filter order is enough: an FIR of
%                       order N cannot see further than N samples.
%   minSegmentSamples - segments shorter than this are copied through UNCHANGED with a
%                       warning. Callers must ensure such segments are excluded
%                       downstream; zapline_covariance_epochs drops them.
%
%   out               - [nSamples x nChannels], same size as data.
%
% Reflection is delegated to RELAX_reflect_pad, which works in [nChannels x nSamples],
% so only one segment at a time is transposed. The vendored tree is only ever reached
% through RELAX_add_vendored_zapline(), so that function is always on the path.
%
% RELAX-GEDAI local addition to zapline-plus. See vendor/zapline-plus/VENDORING.md.

    out = data;

    for k = 1:size(segments, 1)
        segStart = segments(k, 1);
        segEnd   = segments(k, 2);
        segLen   = segEnd - segStart + 1;

        if segLen < minSegmentSamples
            warning('zapline_filter_by_segments:shortSegment', ...
                'Segment %d (%d samples) is shorter than the filter needs; left unfiltered.', ...
                k, segLen);
            continue;
        end

        p        = min(padWidth, segLen - 1);
        padded   = RELAX_reflect_pad(data(segStart:segEnd, :).', p).';
        filtered = fn(padded);
        out(segStart:segEnd, :) = filtered(p+1 : p+segLen, :);
    end
end
