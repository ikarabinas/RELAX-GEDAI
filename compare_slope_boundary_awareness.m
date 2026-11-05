% compare_slope_boundary_awareness() - How many ICA components targeted wICA flags as muscle when the
%                                      slope is fitted on the whole record versus within each
%                                      stretch of contiguous data, for finished wICA runs.
%
% Usage:
%   >> compare_slope_boundary_awareness(runDir, fileIndices, outCsv);
%
% Inputs:
%   runDir      - RELAXProcessed output folder of a targeted wICA run
%   fileIndices - which files, counting runDir's Cleaned_Data alphabetically
%   outCsv      - one row per file: boundaries, components, share flagged as muscle each way
%
% The wICA input is rebuilt from Extremes_Rejected as review_wICA_components.m does. The
% whole-record share is checked against the one the run stored.

function compare_slope_boundary_awareness(runDir, fileIndices, outCsv)

    cleanedFiles = dir(fullfile(runDir, 'Cleaned_Data', '*.set'));
    rows = table();
    for fileIndex = fileIndices(:)'
        cleanedName = cleanedFiles(fileIndex).name;
        fileName = extractBefore(cleanedName, '_RELAX');
        loaded = load(fullfile(runDir, 'RELAX_cfg', [fileName '.mat']));
        cfgField = fieldnames(loaded);
        RELAX_cfg = loaded.(cfgField{1});

        cleaned = pop_loadset('filepath', fullfile(runDir, 'Cleaned_Data'), 'filename', cleanedName, 'loadmode', 'info');
        wICAinput = pop_loadset('filepath', fullfile(runDir, 'Extremes_Rejected'), 'filename', [fileName '_Extremes_Rejected.set']);
        if ~isempty(wICAinput.RELAX.ExtremelyBadPeriodsForDeletion)
            wICAinput = eeg_eegrej(wICAinput, wICAinput.RELAX.ExtremelyBadPeriodsForDeletion);
        end
        wICAinput = RELAX_average_rereference(wICAinput);
        icaLabels = {cleaned.chanlocs(cleaned.icachansind).labels};
        [inSavedICA, columns] = ismember({wICAinput.chanlocs.labels}, icaLabels);
        wICAinput = pop_select(wICAinput, 'channel', find(inSavedICA));
        unmixing = cleaned.icaweights * cleaned.icasphere;
        unmixing = unmixing(:, columns(inSavedICA));
        activations = unmixing * double(wICAinput.data);

        segments = RELAX_get_segments(wICAinput);
        wholeRecord = RELAX_IC_muscle_slopes(activations, wICAinput.srate, RELAX_cfg.LineNoiseFrequency);
        bySegment = RELAX_IC_muscle_slopes(activations, wICAinput.srate, RELAX_cfg.LineNoiseFrequency, segments);
        threshold = RELAX_cfg.MuscleSlopeThreshold;

        rows = [rows; table({fileName}, size(segments,1) - 1, size(activations,1), ...
            mean(wholeRecord >= threshold), mean(bySegment >= threshold), ...
            sum(wholeRecord < threshold & bySegment >= threshold), sum(wholeRecord >= threshold & bySegment < threshold), ...
            median(bySegment - wholeRecord), cleaned.RELAXProcessing_wICA.Proportion_ICs_muscle_by_slope, ...
            'VariableNames', {'File','Boundaries','Components','MuscleShareWholeRecord','MuscleShareBySegment', ...
            'GainedMuscle','LostMuscle','MedianSlopeShift','StoredMuscleShare'})]; %#ok<AGROW>
        disp(rows(end,:))
    end
    writetable(rows, outCsv);
end
