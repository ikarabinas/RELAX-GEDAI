%% find_missing_files.m
% Compares input .set files against RELAXProcessed output files to identify
% which input files were not successfully processed.
% Returns: "missing_files_report.csv" and "missing_files_only.csv"

%% ---- Paths ----
data_dir = '/athena/grosenicklab/scratch/imk2003/acc_tmseeg/eeg_data/RELAX_GEDAI/RELAX_twICA_GEDAI-dmpfc';
output_dir = fullfile(data_dir, 'RELAXProcessed', 'Cleaned_Data');
csv_outpath = fullfile(data_dir, 'RELAXProcessed', 'missing_files_report.csv');

%% ---- Load input files ----
cd(data_dir);
inputDirList = dir('*.set');
inputFiles = {inputDirList.name};
fprintf('Total input .set files found: %d\n', numel(inputFiles));

%% ---- Load output (cleaned) files ----
if ~exist(output_dir, 'dir')
    error('Output directory not found: %s', output_dir);
end
outputDirList = dir(fullfile(output_dir, '*RELAX*.set'));
outputFiles = {outputDirList.name};
fprintf('Total output (cleaned) .set files found: %d\n\n', numel(outputFiles));

%% ---- Extract participant IDs from input files (for reporting) ----
participantIDs = cell(size(inputFiles));
for f = 1:numel(inputFiles)
    tok = regexp(inputFiles{f}, '^([A-Za-z]\d+)', 'tokens');
    if isempty(tok)
        participantIDs{f} = 'UNPARSED';
    else
        participantIDs{f} = tok{1}{1};
    end
end

%% ---- Match each input file to output files by basename containment ----
% RELAX appends suffixes to the original filename when saving cleaned
% output (e.g. "<basename>_RELAX.set" or similar), so matching is done by
% checking whether the input file's basename appears within any output
% filename - mirroring the same 'contains' logic used in the pipeline's
% own skip-logic.

isProcessed = false(size(inputFiles));
matchedOutputFile = cell(size(inputFiles));

for f = 1:numel(inputFiles)
    [~, basename, ~] = fileparts(inputFiles{f});
    matchIdx = find(contains(outputFiles, basename), 1);
    if ~isempty(matchIdx)
        isProcessed(f) = true;
        matchedOutputFile{f} = outputFiles{matchIdx};
    else
        matchedOutputFile{f} = '';
    end
end

missingFiles = inputFiles(~isProcessed);
missingParticipantIDs = participantIDs(~isProcessed);

fprintf('Input files with a matching cleaned output: %d\n', sum(isProcessed));
fprintf('Input files MISSING from output (failed/unprocessed): %d\n\n', numel(missingFiles));

%% ---- Summarize missing files by participant ----
if ~isempty(missingParticipantIDs)
    uniqueMissingParticipants = unique(missingParticipantIDs, 'stable');
    fprintf('Missing files span %d unique participant(s):\n', numel(uniqueMissingParticipants));
    for p = 1:numel(uniqueMissingParticipants)
        pid = uniqueMissingParticipants{p};
        nMissing = sum(strcmp(missingParticipantIDs, pid));
        fprintf('  %s: %d missing file(s)\n', pid, nMissing);
    end
end

%% ---- Write CSV report ----
% Report includes ALL input files, with a flag column for whether each
% was found in the output directory, so the CSV is a complete audit trail
% (not just the missing ones) - easier to sanity check totals against.

InputFile = inputFiles(:);
ParticipantID = participantIDs(:);
Processed = isProcessed(:);
MatchedOutputFile = matchedOutputFile(:);

resultsTable = table(InputFile, ParticipantID, Processed, MatchedOutputFile);

% Sort so missing files appear first for easy review
resultsTable = sortrows(resultsTable, 'Processed', 'ascend');

writetable(resultsTable, csv_outpath);
fprintf('\nFull report (all %d input files) written to:\n  %s\n', numel(inputFiles), csv_outpath);

% Also write a second CSV containing ONLY the missing files, for convenience
missing_only_csv = fullfile(data_dir, 'RELAXProcessed', 'missing_files_only.csv');
missingTable = resultsTable(~resultsTable.Processed, {'InputFile','ParticipantID'});
writetable(missingTable, missing_only_csv);
fprintf('Missing-only report (%d files) written to:\n  %s\n', height(missingTable), missing_only_csv);