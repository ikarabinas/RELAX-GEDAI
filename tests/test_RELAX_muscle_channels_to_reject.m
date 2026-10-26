function test_RELAX_muscle_channels_to_reject()
% Choosing which muscle-contaminated channels to reject, worst first, within a budget.

    labels = {'E1','E2','E3','E4','E5'};
    proportions = [0.10 0.80 0.55 0.50 0.95];

    % Strictly over the threshold, worst first
    rejected = RELAX_muscle_channels_to_reject(labels, proportions, 0.50, 10);
    assert(isequal(rejected, {'E5','E2','E3'}), 'over-threshold channels, worst first');

    % The budget keeps only the worst
    rejected = RELAX_muscle_channels_to_reject(labels, proportions, 0.50, 2);
    assert(isequal(rejected, {'E5','E2'}), 'budget keeps the worst two');

    % No budget, nothing rejected, and an empty cell rather than []
    rejected = RELAX_muscle_channels_to_reject(labels, proportions, 0.50, 0);
    assert(iscell(rejected) && isempty(rejected), 'zero budget rejects nothing');
    rejected = RELAX_muscle_channels_to_reject(labels, proportions, 0.50, -3);
    assert(iscell(rejected) && isempty(rejected), 'an overspent budget rejects nothing');

    % A threshold of 1 is the configured off switch
    rejected = RELAX_muscle_channels_to_reject(labels, [1 1 1 1 1], 1, 10);
    assert(isempty(rejected), 'threshold 1 never rejects');

    % NaN proportions (no usable epochs) are never rejected
    rejected = RELAX_muscle_channels_to_reject(labels, [NaN 0.9 NaN 0 0], 0.50, 10);
    assert(isequal(rejected, {'E2'}), 'NaN proportions are skipped');

    % Column proportions, as RELAX_metrics_muscle returns them, give the same answer
    rejected = RELAX_muscle_channels_to_reject(labels, proportions', 0.50, 10);
    assert(isequal(rejected, {'E5','E2','E3'}), 'column proportions accepted');

    fprintf('test_RELAX_muscle_channels_to_reject: PASS\n');
end
