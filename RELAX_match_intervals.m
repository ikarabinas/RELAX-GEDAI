% RELAX_match_intervals() - Which intervals in one list overlap any interval in another.
%
% Usage:
%   >> [aHasB, bHasA] = RELAX_match_intervals(aLo, aHi, bLo, bHi, slack);
%
% Inclusive 1-based sample intervals; a is widened by slack samples on both sides before
% the overlap test (so slack = 50 ms in samples pairs a reference pop with a per-channel
% common-mode event that was binned a few samples away).
%
% Outputs: aHasB (numel(aLo) x 1 logical), bHasA (numel(bLo) x 1 logical).

function [aHasB, bHasA] = RELAX_match_intervals(aLo, aHi, bLo, bHi, slack)
    aLo = aLo(:) - slack; aHi = aHi(:) + slack;
    bLo = bLo(:); bHi = bHi(:);
    overlap = aLo <= bHi' & aHi >= bLo';      % numel(a) x numel(b)
    aHasB = reshape(any(overlap, 2), [], 1);
    bHasA = reshape(any(overlap, 1), [], 1);
    if isempty(bLo); aHasB = false(numel(aLo), 1); end
    if isempty(aLo); bHasA = false(numel(bLo), 1); end
end
