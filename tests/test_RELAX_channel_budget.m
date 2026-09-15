function test_RELAX_channel_budget()
% Shared running channel-rejection budget across screens.

    mkEEG = @(n) struct('nbchan', n, 'chanlocs', struct('labels', ...
        arrayfun(@(k) sprintf('E%d', k), 1:n, 'UniformOutput', false)));

    % init records the montage size and the derived cap
    EEG = RELAX_channel_budget(mkEEG(100), 'init', 0.10);
    b = EEG.RELAX.ChannelBudget;
    assert(b.nOriginalChannels == 100, 'nOriginalChannels');
    assert(b.maxRejectable == 10, 'maxRejectable should be floor(0.10*100)');
    assert(b.nSpent == 0, 'nothing spent yet');
    assert(~b.exhausted, 'not exhausted yet');

    % A request inside the budget passes through untouched
    [fits, EEG] = RELAX_channel_budget(EEG, 'spend', 'pop', {'E1','E2','E3'});
    assert(isequal(fits, {'E1','E2','E3'}), 'in-budget request should pass through whole');
    assert(EEG.RELAX.ChannelBudget.nSpent == 3, 'nSpent should be 3');
    assert(isequal(EEG.RELAX.ChannelBudget.spentBy.pop, {'E1','E2','E3'}), 'attribution by stage');

    % A second screen sees the REDUCED remainder, and truncation keeps the worst first
    [fits, EEG] = RELAX_channel_budget(EEG, 'spend', 'impedance', ...
        {'E10','E11','E12','E13','E14','E15','E16','E17','E18','E19'});
    assert(numel(fits) == 7, 'only 7 of the 10 should fit into the remaining budget');
    assert(isequal(fits, {'E10','E11','E12','E13','E14','E15','E16'}), ...
        'truncation must keep the worst-first prefix');
    assert(EEG.RELAX.ChannelBudget.nSpent == 10, 'budget fully spent');
    assert(EEG.RELAX.ChannelBudget.exhausted, 'exhausted flag must be set');

    % A third screen gets nothing, and gets an empty cell rather than []
    [fits, EEG] = RELAX_channel_budget(EEG, 'spend', 'PREP', {'E20'});
    assert(iscell(fits) && isempty(fits), 'exhausted budget must return an empty cell');
    assert(EEG.RELAX.ChannelBudget.nSpent == 10, 'nSpent must not grow past the cap');

    % Spending nothing is legal and does not set the exhausted flag
    EEG2 = RELAX_channel_budget(mkEEG(100), 'init', 0.10);
    [fits, EEG2] = RELAX_channel_budget(EEG2, 'spend', 'pop', {});
    assert(iscell(fits) && isempty(fits), 'empty request returns empty cell');
    assert(~EEG2.RELAX.ChannelBudget.exhausted, 'empty request must not exhaust the budget');

    % A budget that rounds down to zero rejects everything requested
    EEG3 = RELAX_channel_budget(mkEEG(5), 'init', 0.10);
    assert(EEG3.RELAX.ChannelBudget.maxRejectable == 0, 'floor(0.10*5) is 0');
    [fits, ~] = RELAX_channel_budget(EEG3, 'spend', 'pop', {'E1'});
    assert(isempty(fits), 'a zero budget grants nothing');

    % Init must not disturb an existing EEG.RELAX
    EEG4 = mkEEG(50);
    EEG4.RELAX.SomethingElse = 42;
    EEG4 = RELAX_channel_budget(EEG4, 'init', 0.10);
    assert(EEG4.RELAX.SomethingElse == 42, 'init must not clobber other EEG.RELAX fields');

    % --- the impedance screen's backward compatibility ---
    % RELAX_impedance_bad_channels reads kOhm from chanlocs(k).impedance, and returns
    % a COLUMN cell of labels (see tests/test_RELAX_impedance_bad_channels.m).
    chanlocs = struct('labels', {'E1','E2','E3','E4','E5'}, ...
                      'impedance', {900, 800, 700, 10, 10});

    % 3-argument form: cap is floor(0.10*5) = 0, so nothing is rejected
    rej3 = RELAX_impedance_bad_channels(chanlocs, 500, 0.10);
    assert(isempty(rej3), '3-arg form must keep its existing floor(maxProportion*n) cap');

    % 4-argument form: the explicit count overrides, worst impedance first
    rej4 = RELAX_impedance_bad_channels(chanlocs, 500, 0.10, 2);
    assert(numel(rej4) == 2, '4-arg form should honour the explicit count');
    assert(all(ismember({'E1','E2'}, rej4)), 'the two worst impedances must be the ones kept');

    % An explicit count large enough to not bind returns everything over threshold
    rej5 = RELAX_impedance_bad_channels(chanlocs, 500, 0.10, 99);
    assert(numel(rej5) == 3, 'a non-binding count should return all 3 over-threshold channels');

    fprintf('test_RELAX_channel_budget: PASS\n');
end
