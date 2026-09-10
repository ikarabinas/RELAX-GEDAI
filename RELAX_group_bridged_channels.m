% RELAX_group_bridged_channels() - Group eBridge-flagged channel pairs into
%                                   connected clusters of electrically bridged channels.
%
% Usage:
%   >> [clusters, pairLabels] = RELAX_group_bridged_channels(labels, pairPairs);
%
% Inputs:
%   labels     - cell array of channel labels, in EEG.chanlocs order
%   pairPairs  - eBridge's EB.Bridged.Pairs: 2 x nPairs cell array of channel label
%                pairs flagged as bridged. A pair naming a label not present in
%                `labels` (e.g. a channel already removed upstream) is silently
%                dropped, matching the source pipeline's behaviour.
%
% Outputs:
%   clusters   - 1 x nClusters cell array; each element is a cell array of channel
%                labels electrically shorted together (connected components of the
%                flagged-pair graph, so bridge chains of more than 2 channels are
%                grouped rather than split into separate pairs)
%   pairLabels - 1 x nMatchedPairs cell array of 'A-B' strings, one per matched pair
%                (pairs with an unmatched label excluded), for QC reporting
%
% Bridges chain: electrolyte spreading across a patch commonly shorts 3-5 electrodes
% into one cluster, and handling flagged pairs independently would miss that they are
% one measurement, not several unrelated ones. Split out of RELAX_bridged_channels so
% this grouping logic can be unit-tested without the eBridge plugin installed.

function [clusters, pairLabels] = RELAX_group_bridged_channels(labels, pairPairs)

    clusters = {};
    pairLabels = {};

    numChans = numel(labels);
    numPairs = size(pairPairs, 2);
    if numChans == 0 || numPairs == 0
        return;
    end

    adjacency = false(numChans);
    pairLabels = cell(1, numPairs);
    for k = 1:numPairs
        i = find(strcmp(labels, pairPairs{1,k}), 1);
        j = find(strcmp(labels, pairPairs{2,k}), 1);
        if isempty(i) || isempty(j) || i == j
            continue; % a self-pair would set a nonzero adjacency diagonal, which graph() rejects
        end
        adjacency(i,j) = true;
        adjacency(j,i) = true;
        pairLabels{k} = sprintf('%s-%s', labels{i}, labels{j});
    end
    pairLabels = pairLabels(~cellfun('isempty', pairLabels));

    comp = conncomp(graph(adjacency));
    for c = 1:max(comp)
        members = find(comp == c);
        if numel(members) < 2
            continue; % an unbridged channel sitting in its own component
        end
        clusters{end+1} = labels(members); %#ok<AGROW>
    end
end
