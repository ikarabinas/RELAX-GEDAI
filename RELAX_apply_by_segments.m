% RELAX_apply_by_segments() - Apply a function to each segment of data independently,
%                             with reflection padding at each segment edge.
%
% Usage:
%   >> out = RELAX_apply_by_segments(data, segments, fn, p);
%
% Inputs:
%   data     - [nChannels x nSamples] data
%   segments - [nSegments x 2] matrix of [startSample endSample], as returned by
%              RELAX_get_segments
%   fn       - function handle taking [nChannels x N] and returning [nChannels x N]
%   p        - reflection pad width in samples, clamped per segment to segLen-1
%
% Outputs:
%   out - [nChannels x nSamples], the same size as data. Any samples not covered by
%         a segment are passed through unchanged.
%
% Each segment is padded, transformed and trimmed on its own, so no segment can
% influence any other. Callers must derive `segments` from the CURRENT EEG.event:
% eeg_eegrej adds boundaries part-way through the pipeline, so a segment list must
% never be cached across stages.

function out = RELAX_apply_by_segments(data, segments, fn, p)

    out = data;

    for s = 1:size(segments,1)
        segStart = segments(s,1);
        segEnd   = segments(s,2);
        segData  = data(:, segStart:segEnd);

        thisPad  = min(p, size(segData,2)-1);
        if thisPad < 0
            thisPad = 0;
        end

        padded   = RELAX_reflect_pad(segData, thisPad);
        filtered = fn(padded);

        out(:, segStart:segEnd) = filtered(:, thisPad+1 : thisPad+size(segData,2));
    end
end
