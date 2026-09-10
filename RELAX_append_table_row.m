% RELAX_append_table_row() - Assign a row into an aggregate table, reconciling a
%                            schema that has changed between versions.
%
% Usage:
%   >> T = RELAX_append_table_row(T, RowNumber, NewRow);
%
% Inputs:
%   T         - Existing table, or [] / an empty table on the first row of a run.
%   RowNumber - Row index to assign into.
%   NewRow    - 1-row table, typically struct2table(...,'AsArray',true).
%
% Outputs:
%   T - The table with NewRow assigned at RowNumber.
%
% The wrapper reloads per-participant aggregate tables checkpointed by earlier runs so a
% run can resume. A table written before a column was added carries a different variable
% set, and plain row assignment errors on that. Rather than discard rows the table already
% holds, this fills the missing columns on whichever side lacks them.
%
% The ordering step is the subtle part: table row assignment matches variables by POSITION,
% not by name. Widening both sides to the same SET of names is not enough - if the orders
% differ, MATLAB lines a cell column up against a numeric one and either errors with a
% confusing "Conversion to double from cell" or, worse, writes the wrong column. The row is
% therefore reordered to the table's variable order before assignment.
%
% Backfilled values are {} for cell columns and NaN for everything else, so "this run did
% not record that" is distinguishable from a recorded zero.

function T = RELAX_append_table_row(T, RowNumber, NewRow)

    if isempty(T) || ~istable(T)
        T = NewRow;
        if RowNumber ~= 1
            % Preserve the caller's indexing when the first file processed is not file 1
            T = NewRow([], :);
            T(RowNumber, :) = NewRow;
        end
        return;
    end

    ColumnsMissingFromTable = setdiff(NewRow.Properties.VariableNames, ...
        T.Properties.VariableNames, 'stable');
    for MissingColumn = 1:numel(ColumnsMissingFromTable)
        ColumnName = ColumnsMissingFromTable{MissingColumn};
        if iscell(NewRow.(ColumnName))
            T.(ColumnName) = cell(height(T), 1);
        else
            T.(ColumnName) = NaN(height(T), 1);
        end
    end

    ColumnsMissingFromRow = setdiff(T.Properties.VariableNames, ...
        NewRow.Properties.VariableNames, 'stable');
    for MissingColumn = 1:numel(ColumnsMissingFromRow)
        ColumnName = ColumnsMissingFromRow{MissingColumn};
        if iscell(T.(ColumnName))
            NewRow.(ColumnName) = cell(1, 1);
        else
            NewRow.(ColumnName) = NaN;
        end
    end

    % Match by name, not by position - see the note above:
    NewRow = NewRow(:, T.Properties.VariableNames);

    T(RowNumber, :) = NewRow;
end
