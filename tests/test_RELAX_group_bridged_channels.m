function test_RELAX_group_bridged_channels()
% Connected-components grouping of eBridge-flagged pairs.

    % --- A single simple pair becomes a 2-member cluster ---
    labels = {'E1','E2','E3','E4'};
    pairs = {'E1'; 'E2'}; % 2x1: one pair
    [clusters, pairLabels] = RELAX_group_bridged_channels(labels, pairs);
    assert(numel(clusters) == 1 && isequal(sort(clusters{1}), {'E1','E2'}), ...
        'a single flagged pair must become one 2-member cluster');
    assert(isequal(pairLabels, {'E1-E2'}), 'pair label string uses channel order');

    % --- A chain of pairs merges into one cluster, not several ---
    labels = {'E1','E2','E3','E4','E5'};
    pairs = [{'E1','E2','E3','E4'}; {'E2','E3','E4','E5'}]; % chain E1-E2-E3-E4-E5
    [clusters, ~] = RELAX_group_bridged_channels(labels, pairs);
    assert(numel(clusters) == 1, 'a chain of pairs must merge into a single cluster');
    assert(isequal(sort(clusters{1}), {'E1','E2','E3','E4','E5'}), ...
        'the merged cluster must contain every chained channel');

    % --- Two disjoint pairs stay as two separate clusters ---
    labels = {'E1','E2','E3','E4','E5','E6'};
    pairs = [{'E1','E4'}; {'E2','E5'}]; % E1-E2 and E4-E5, unrelated
    [clusters, ~] = RELAX_group_bridged_channels(labels, pairs);
    assert(numel(clusters) == 2, 'unrelated pairs must not be merged');
    sizes = cellfun(@numel, clusters);
    assert(all(sizes == 2), 'each disjoint pair is its own 2-member cluster');

    % --- A channel with no flagged pair never appears in any cluster ---
    labels = {'E1','E2','E3'};
    pairs = {'E1'; 'E2'};
    [clusters, ~] = RELAX_group_bridged_channels(labels, pairs);
    allMembers = [clusters{:}];
    assert(~ismember('E3', allMembers), 'an unpaired channel must not appear in any cluster');

    % --- A pair naming an unknown label is silently skipped ---
    labels = {'E1','E2'};
    pairs = [{'E1','E1'}; {'E2','E99'}]; % second pair references E99, not in labels
    [clusters, pairLabels] = RELAX_group_bridged_channels(labels, pairs);
    assert(numel(clusters) == 1 && isequal(sort(clusters{1}), {'E1','E2'}), ...
        'the valid pair must still be grouped');
    assert(isequal(pairLabels, {'E1-E2'}), 'the pair naming an unknown label must be dropped from pairLabels');

    % --- No pairs -> no clusters ---
    [clusters, pairLabels] = RELAX_group_bridged_channels({'E1','E2'}, cell(2,0));
    assert(isempty(clusters) && isempty(pairLabels), 'empty pair list produces no clusters');

    % --- Empty label list -> no clusters, no error ---
    [clusters, pairLabels] = RELAX_group_bridged_channels({}, cell(2,0));
    assert(isempty(clusters) && isempty(pairLabels), 'empty channel list tolerated');

    fprintf('test_RELAX_group_bridged_channels: PASS\n');
end
