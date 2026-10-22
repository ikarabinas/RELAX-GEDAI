% RELAX_screen_top_separated() - Highest-scoring times that are at least minSeparation apart.
%
% Usage:
%   >> centres = RELAX_screen_top_separated(times, scores, k, minSeparation);
%
% Inputs:
%   times         - [1 x n] or [n x 1] time of each score
%   scores        - vector of the same length as times; NaN scores are ignored
%   k             - most centres to return
%   minSeparation - least distance allowed between any two returned centres, in the units of times
%
% Outputs:
%   centres - [1 x <=k] times in descending score order, chosen greedily so the highest score
%             always wins. Fewer than k when too few candidates are far enough apart

function centres = RELAX_screen_top_separated(times, scores, k, minSeparation)

    scored = ~isnan(scores(:));
    [~, order] = sort(scores(scored), 'descend');
    candidateTimes = times(scored);
    candidateTimes = candidateTimes(order);

    centres = zeros(1, 0);
    while ~isempty(candidateTimes) && numel(centres) < k
        centres(end+1) = candidateTimes(1); %#ok<AGROW>
        candidateTimes = candidateTimes(abs(candidateTimes - centres(end)) >= minSeparation);
    end
end
