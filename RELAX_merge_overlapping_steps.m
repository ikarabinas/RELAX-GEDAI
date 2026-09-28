% RELAX_merge_overlapping_steps() - Collapse steps whose [onset, offset] intervals overlap
%                                   or coincide into one step per run.
%
% Usage:
%   >> [merged, groupOf] = RELAX_merge_overlapping_steps(steps, x, segments, W, srate);
%
% RELAX_detect_steps can bracket one transition more than once: a staircase or a long
% transition carries several local-z peaks, and the same bracket is sometimes reported at
% two z values. Subtracting each bracket on its own would correct the transition twice
% (RELAX_subtract_reference_pops stacks corrections), so every run of overlapping steps
% becomes a single step: the union interval, with the amplitude re-measured on x exactly
% as the detector measures it - the median of the W samples from the offset minus the
% median of the W samples up to the onset, both clipped to the step's segment. Overlapping
% intervals share at least one sample, so they always lie in one segment.
%
% Inputs:
%   steps    - struct of parallel N x 1 numeric columns with at least .sample .onset
%              .offset .amplitude .riseMs .z .channelZ; any other columns are carried
%              through from the member with the largest |z|
%   x        - 1 x pnts trace the steps were detected on
%   segments - S x 2 [start end] inter-boundary segments (RELAX_get_segments)
%   W        - detector window in samples (detectInfo.params.derived.W)
%   srate    - sampling rate, for riseMs of a merged interval
%
% Outputs:
%   merged   - the same columns, one row per run, in onset order
%   groupOf  - N x 1 row of merged that each input step went to

function [merged, groupOf] = RELAX_merge_overlapping_steps(steps, x, segments, W, srate)
    n = numel(steps.onset);
    x = double(x(:))';
    [~, order] = sort(steps.onset);
    groupOf = zeros(n, 1);
    nGroups = 0; groupEnd = -Inf;
    for i = order'
        if steps.onset(i) > groupEnd; nGroups = nGroups + 1; end
        groupOf(i) = nGroups;
        groupEnd = max(groupEnd, steps.offset(i));
    end

    rows = cell(nGroups, 1);
    for g = 1:nGroups
        members = find(groupOf == g);
        [~, leadPos] = max(abs(steps.z(members)));
        row = structfun(@(f) f(members(leadPos)), steps, 'UniformOutput', false);
        if numel(members) > 1
            onset = min(steps.onset(members)); offset = max(steps.offset(members));
            seg = find(segments(:,1) <= onset & segments(:,2) >= offset, 1);
            preWin  = max(segments(seg,1), onset - W + 1) : onset;
            postWin = offset : min(segments(seg,2), offset + W - 1);
            row.onset = onset; row.offset = offset;
            row.amplitude = median(x(postWin), 'omitnan') - median(x(preWin), 'omitnan');
            row.riseMs = (offset - onset - 1) / srate * 1000;
            row.z = max(steps.z(members)); row.channelZ = max(steps.channelZ(members));
        end
        rows{g} = row;
    end

    merged = structfun(@(f) f(zeros(0, 1)), steps, 'UniformOutput', false);
    if nGroups == 0; return; end
    rows = [rows{:}];
    for f = fieldnames(steps)'
        merged.(f{1}) = reshape([rows.(f{1})], [], 1);
    end
end
