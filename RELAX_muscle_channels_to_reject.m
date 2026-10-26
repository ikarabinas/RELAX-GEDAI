% RELAX_muscle_channels_to_reject() - Channels to reject for muscle contamination, worst first.
%
% Usage:
%   >> rejected = RELAX_muscle_channels_to_reject(labels, proportions, threshold, maxRejectable);
%
% Inputs:
%   labels        - cell of channel labels
%   proportions   - per-channel proportion of epochs whose muscle slope exceeds
%                   RELAX_cfg.MuscleSlopeThreshold, same order as labels (row or column)
%   threshold     - RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel;
%                   a channel must be strictly over it
%   maxRejectable - channel-rejection budget still unspent; <= 0 rejects nothing
%
% Outputs:
%   rejected - [1 x k] cell of labels, most contaminated first, at most maxRejectable long
%
% Unlike the epoching stage's version of this screen, a tie at the budget's edge is broken by
% channel order rather than by rejecting none of the tied channels.

function rejected = RELAX_muscle_channels_to_reject(labels, proportions, threshold, maxRejectable)

    proportions = proportions(:)';
    overThreshold = find(proportions > threshold);
    [~, worstFirst] = sort(proportions(overThreshold), 'descend');
    overThreshold = overThreshold(worstFirst);
    rejected = labels(overThreshold(1:min(numel(overThreshold), max(maxRejectable, 0))));
    rejected = reshape(rejected, 1, []);
end
