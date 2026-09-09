% RELAX_apply_by_segments() - Apply a function to each segment of data independently,
%                             with reflection padding at each segment edge.
%
% Usage:
%   >> out = RELAX_apply_by_segments(data, segments, fn, p);
%   >> out = RELAX_apply_by_segments(data, segments, fn, p, minSegmentLength);
%
% Inputs:
%   data             - [nChannels x nSamples] data
%   segments         - [nSegments x 2] matrix of [startSample endSample], as returned
%                       by RELAX_get_segments. Must tile the record exactly: segment 1
%                       starts at sample 1, each segment starts one sample after the
%                       previous one ends, and the last segment ends at the last sample.
%   fn               - function handle taking [nChannels x N] and returning
%                       [nChannels x N]
%   p                - reflection pad width in samples, clamped per segment to
%                       segLen-1
%   minSegmentLength - (optional, default 1) segments shorter than this are copied
%                       through UNCHANGED rather than passed to fn, with a warning
%                       naming the segment index and length. The default of 1 means
%                       fn is applied to EVERY segment - this helper is generic (fn
%                       can be anything from a pointwise mask op to a filtfilt call),
%                       so it imposes no length floor of its own. Callers whose fn
%                       has a real minimum-length requirement (anything built on
%                       filtfilt, or a wavelet transform) must pass that floor
%                       explicitly - e.g. an order-4 filtfilt needs a signal longer
%                       than 3x its order, and coif5's level-5 support is ~480
%                       samples, so a sliver segment handed to either can crash or
%                       corrupt it. A caller with no such requirement (e.g. a plain
%                       mask/pointwise op) should leave this at the default so it
%                       keeps applying fn to every segment, including short ones.
%
% Outputs:
%   out - [nChannels x nSamples], the same size as data. Any samples not covered by
%         a segment are passed through unchanged.
%
% Each segment is padded, transformed and trimmed on its own, so no segment can
% influence any other. Callers must derive `segments` from the CURRENT EEG.event:
% eeg_eegrej adds boundaries part-way through the pipeline, so a segment list must
% never be cached across stages.

function out = RELAX_apply_by_segments(data, segments, fn, p, minSegmentLength)

    if nargin < 5 || isempty(minSegmentLength)
        minSegmentLength = 1; % apply fn to every segment unless a caller asks otherwise
    end

    assert(~isempty(segments) && segments(1,1) == 1 && segments(end,2) == size(data,2), ...
        'RELAX_apply_by_segments:segmentsMustTile', ...
        'segments must tile the record: first segment must start at sample 1 and last must end at sample %d.', ...
        size(data,2));
    assert(all(segments(2:end,1) == segments(1:end-1,2) + 1), ...
        'RELAX_apply_by_segments:segmentsMustTile', ...
        'segments must tile the record: each segment must start one sample after the previous one ends.');

    out = data;

    for s = 1:size(segments,1)
        segStart = segments(s,1);
        segEnd   = segments(s,2);
        segData  = data(:, segStart:segEnd);
        segLen   = size(segData,2);

        if segLen < minSegmentLength
            warning('RELAX_apply_by_segments:segmentTooShort', ...
                'Segment %d (samples %d:%d, length %d) is shorter than minSegmentLength=%d; passing it through unchanged.', ...
                s, segStart, segEnd, segLen, minSegmentLength);
            continue % out(:, segStart:segEnd) already holds the unchanged data
        end

        thisPad  = min(p, segLen-1);
        if thisPad < 0
            thisPad = 0;
        end

        padded   = RELAX_reflect_pad(segData, thisPad);
        filtered = fn(padded);

        out(:, segStart:segEnd) = filtered(:, thisPad+1 : thisPad+segLen);
    end
end
