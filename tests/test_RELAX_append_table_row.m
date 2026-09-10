function test_RELAX_append_table_row()
% Schema-reconciling row assignment into the aggregate tables.

    oldRow = @(name, chans, retained) struct2table(struct( ...
        'aFileName', {{name}}, 'PREPBasedChannelToReject', {{chans}}, ...
        'RetainedChannelCount', retained), 'AsArray', true);
    % Note the new column sits in the MIDDLE, so a name-set fix that ignores order fails
    newRow = @(name, chans, imp, retained) struct2table(struct( ...
        'aFileName', {{name}}, 'PREPBasedChannelToReject', {{chans}}, ...
        'ImpedanceBasedChannelToReject', {{imp}}, ...
        'RetainedChannelCount', retained), 'AsArray', true);

    % --- Empty start: [] becomes the row ---
    T = RELAX_append_table_row([], 1, newRow('f1', {'E1'}, {'E9'}, 236));
    assert(height(T) == 1 && strcmp(T.aFileName{1}, 'f1'), 'first row on an empty table');

    % --- Empty start at a row other than 1 keeps the caller's indexing ---
    T = RELAX_append_table_row([], 3, newRow('f3', {'E1'}, {'E9'}, 236));
    assert(height(T) == 3, 'assigning row 3 into an empty table gives height 3');
    assert(strcmp(T.aFileName{3}, 'f3'), 'the row lands at the requested index');

    % --- Matching schema behaves like plain assignment ---
    T = RELAX_append_table_row([], 1, newRow('f1', {'E1'}, {'E9'}, 236));
    T = RELAX_append_table_row(T, 2, newRow('f2', {'E2'}, {'E8'}, 240));
    assert(height(T) == 2 && strcmp(T.aFileName{2}, 'f2'), 'matching schema appends');
    assert(T.RetainedChannelCount(2) == 240, 'numeric column intact');

    % --- Old table + new row: the table gains the column, existing rows survive ---
    T = RELAX_append_table_row([], 1, oldRow('f1', {'E1'}, 250));
    T = RELAX_append_table_row(T, 2, newRow('f2', {'E2'}, {'E9';'E10'}, 236));
    assert(height(T) == 2, 'both rows kept');
    assert(strcmp(T.aFileName{1}, 'f1'), 'the pre-existing row is preserved');
    assert(ismember('ImpedanceBasedChannelToReject', T.Properties.VariableNames), 'column added');
    assert(isequal(T.ImpedanceBasedChannelToReject{2}, {'E9';'E10'}), 'new row written correctly');
    assert(isempty(T.ImpedanceBasedChannelToReject{1}), 'the old row is backfilled empty');
    % The regression this guard exists for: a positional mismatch would put the cell
    % column into RetainedChannelCount and either error or corrupt it
    assert(isnumeric(T.RetainedChannelCount), 'numeric column must not become a cell');
    assert(T.RetainedChannelCount(1) == 250 && T.RetainedChannelCount(2) == 236, ...
        'columns must be matched by name, not position');

    % --- New table + old row: the row is backfilled instead of blocking assignment ---
    T = RELAX_append_table_row([], 1, newRow('f1', {'E1'}, {'E9'}, 236));
    T = RELAX_append_table_row(T, 2, oldRow('f2', {'E2'}, 250));
    assert(height(T) == 2, 'both rows kept in the reverse direction');
    assert(isempty(T.ImpedanceBasedChannelToReject{2}), 'the older row backfills empty');
    assert(T.RetainedChannelCount(2) == 250, 'numeric column still matched by name');

    % --- Reordered but otherwise identical schema is matched by name ---
    reordered = struct2table(struct('RetainedChannelCount', 231, ...
        'ImpedanceBasedChannelToReject', {{{'E7'}}}, 'PREPBasedChannelToReject', {{{'E3'}}}, ...
        'aFileName', {{'f3'}}), 'AsArray', true);
    T = RELAX_append_table_row(T, 3, reordered);
    assert(strcmp(T.aFileName{3}, 'f3'), 'reordered row lands in the right columns');
    assert(T.RetainedChannelCount(3) == 231, 'reordered numeric column matched by name');
    assert(isequal(T.ImpedanceBasedChannelToReject{3}, {'E7'}), 'reordered cell column matched by name');

    fprintf('test_RELAX_append_table_row: PASS\n');
end
