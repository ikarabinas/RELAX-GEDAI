% review_wICA_components() - Per-component record of how targeted wICA judged each ICA component
%                            in finished RELAX_Wrapper_slurm.m runs: ICLabel probabilities and
%                            thresholded label, the muscle spectral slope, and icablinkmetrics.
%
% Usage:
%   >> review_wICA_components('extract', runDirs, runNames, outDir, fileIndex);
%   >> review_wICA_components('combine', runDirs, runNames, outDir);
%
% Inputs:
%   runDirs   - cellstr of RELAXProcessed output folders of targeted wICA runs over the same files
%   runNames  - cellstr of short names for the runs, same order
%   outDir    - folder for per_file/*.mat, wICA_components_<run>.csv, summary.md and figures
%   fileIndex - 'extract' only: which file, counting the first run's Cleaned_Data alphabetically
%
% The ICLabel probabilities are the ones the run stored on the cleaned file. The slopes and
% icablinkmetrics are recomputed on the wICA input, rebuilt from Extremes_Rejected as the wrapper
% built it (excise ExtremelyBadPeriodsForDeletion, average re-reference; no interpolation, so
% ICA's own channels), using the saved unmixing matrix. The recomputed muscle-by-slope share is
% checked against the one the run stored, so a mismatch shows the rebuild is wrong. Channels the
% deferred muscle screen rejected after wICA are no longer in the saved ICA; they are dropped from
% the rebuilt input too, the column ChannelsDroppedFromICA counts them, and those files' slopes and
% blink metrics are approximate.
%
% icablinkmetrics did not run in runs whose MATLAB path lacked it (the wrapper skips it
% silently); then the blink columns say what it would have flagged, and BlinkMetricsRanInRun is
% false. Each run's RELAX_cfg/<file>.mat supplies the thresholds.

function review_wICA_components(mode, runDirs, runNames, outDir, fileIndex)

    perFileDir = fullfile(outDir, 'per_file');
    if ~exist(perFileDir, 'dir'); mkdir(perFileDir); end

    switch mode
        case 'extract'
            cleanedFiles = dir(fullfile(runDirs{1}, 'Cleaned_Data', '*.set'));
            cleanedName = cleanedFiles(fileIndex).name;
            fileName = extractBefore(cleanedName, '_RELAX');
            components = table();
            for r = 1:numel(runDirs)
                components = [components; reviewFile(runDirs{r}, runNames{r}, fileName, cleanedName)]; %#ok<AGROW>
            end
            save(fullfile(perFileDir, [fileName '.mat']), 'components');
        case 'combine'
            combine(runNames, outDir, perFileDir);
    end
end

function components = reviewFile(runDir, runName, fileName, cleanedName)

    loaded = load(fullfile(runDir, 'RELAX_cfg', [fileName '.mat']));
    cfgField = fieldnames(loaded);
    RELAX_cfg = loaded.(cfgField{1});

    cleaned = pop_loadset('filepath', fullfile(runDir, 'Cleaned_Data'), 'filename', cleanedName, 'loadmode', 'info');
    wICAinput = pop_loadset('filepath', fullfile(runDir, 'Extremes_Rejected'), ...
        'filename', [fileName '_Extremes_Rejected.set']);
    if ~isempty(wICAinput.RELAX.ExtremelyBadPeriodsForDeletion)
        wICAinput = eeg_eegrej(wICAinput, wICAinput.RELAX.ExtremelyBadPeriodsForDeletion);
    end
    wICAinput = RELAX_average_rereference(wICAinput);

    icaLabels = {cleaned.chanlocs(cleaned.icachansind).labels};
    [inSavedICA, columns] = ismember({wICAinput.chanlocs.labels}, icaLabels);
    channelsDroppedFromICA = sum(~inSavedICA);
    wICAinput = pop_select(wICAinput, 'channel', find(inSavedICA));
    assert(wICAinput.nbchan == numel(icaLabels), 'channels of the saved ICA missing from the wICA input: %s', ...
        strjoin(setdiff(icaLabels, {wICAinput.chanlocs.labels}), ', '));
    unmixing = cleaned.icaweights * cleaned.icasphere;
    unmixing = unmixing(:, columns(inSavedICA));
    mixing = pinv(unmixing);

    activations = unmixing * double(wICAinput.data);
    nComponents = size(activations, 1);
    slopes = RELAX_IC_muscle_slopes(activations, wICAinput.srate, RELAX_cfg.LineNoiseFrequency);
    muscleBySlope = slopes >= RELAX_cfg.MuscleSlopeThreshold;

    probabilities = cleaned.etc.ic_classification.ICLabel.classifications;
    classes = cleaned.etc.ic_classification.ICLabel.classes;
    [~, unthresholdedClass] = max(probabilities, [], 2);
    thresholded = probabilities;
    thresholded(thresholded < RELAX_cfg.ICLabel_thresholds) = 0;
    [~, thresholdedClass] = max(thresholded, [], 2);

    inputPower = sum(double(wICAinput.data(:)).^2);
    varianceShare = sum(mixing.^2, 1)' .* sum(activations.^2, 2) / inputPower;

    blink = icablinkmetricsColumns(wICAinput, unmixing, mixing, activations, RELAX_cfg);

    storedShare = cleaned.RELAXProcessing_wICA.Proportion_ICs_muscle_by_slope;

    components = table();
    components.Run = repmat(string(runName), nComponents, 1);
    components.File = repmat(string(fileName), nComponents, 1);
    components.Component = (1:nComponents)';
    components.VarianceShare = varianceShare;
    for c = 1:numel(classes)
        components.(['P_' strrep(classes{c}, ' ', '')]) = probabilities(:, c);
    end
    components.ICLabelTop = string(classes(unthresholdedClass))';
    components.ICLabelThresholded = string(classes(thresholdedClass))';
    components.MuscleSlope = slopes';
    components.MuscleSlopeThreshold = repmat(RELAX_cfg.MuscleSlopeThreshold, nComponents, 1);
    components.MuscleBySlope = muscleBySlope';
    components = [components, blink];
    components.BlinkMetricsRanInRun = repmat(isfield(cleaned.RELAXProcessing_wICA, 'BlinkMetricsRan') ...
        && cleaned.RELAXProcessing_wICA.BlinkMetricsRan, nComponents, 1);
    isEye = components.ICLabelThresholded == "Eye";
    isICLabelMuscle = components.ICLabelThresholded == "Muscle";
    if strcmp(RELAX_cfg.Clean_other_comps, 'yes')
        isWaveletCleaned = components.ICLabelThresholded ~= "Brain";
    else
        isWaveletCleaned = isEye | isICLabelMuscle;
    end
    components.WaveletCleanedBlinkPeriodsOnly = isEye;
    components.WaveletCleanedWholeRecording = isWaveletCleaned & ~isEye;
    components.HighFrequencySubtracted = components.MuscleBySlope;
    components.StoredMuscleBySlopeShare = repmat(storedShare, nComponents, 1);
    components.ChannelsDroppedFromICA = repmat(channelsDroppedFromICA, nComponents, 1);
end

function blink = icablinkmetricsColumns(wICAinput, unmixing, mixing, activations, RELAX_cfg)

    nComponents = size(activations, 1);
    blink = table(nan(nComponents, 1), nan(nComponents, 1), nan(nComponents, 1), ...
        nan(nComponents, 1), nan(nComponents, 1), nan(nComponents, 1), false(nComponents, 1), ...
        repmat("", nComponents, 1), 'VariableNames', {'Blink_Correlation', 'Blink_CorrelationP', ...
        'Blink_Convolution', 'Blink_ConvolutionP', 'Blink_PercentReduction', 'Blink_PercentReductionP', ...
        'BlinkMetricsIdentified', 'BlinkMetricsError'});

    withICA = wICAinput;
    withICA.icaweights = unmixing;
    withICA.icasphere = eye(wICAinput.nbchan);
    withICA.icawinv = mixing;
    withICA.icachansind = 1:wICAinput.nbchan;
    withICA.icaact = activations;
    allElectrodes = pop_interp(wICAinput, wICAinput.allchan, 'spherical');
    blinkChannel = strcmpi({allElectrodes.chanlocs.labels}, RELAX_cfg.BlinkElectrodes(1,1));
    try
        out = icablinkmetrics(withICA, 'ArtifactChannel', allElectrodes.data(blinkChannel, :), ...
            'Alpha', 0.001, 'VisualizeData', 'False');
        blink.Blink_Correlation = [out.metrics.correlation]';
        blink.Blink_CorrelationP = [out.metrics.corr_Pvalue]';
        blink.Blink_Convolution = [out.metrics.convolution]';
        blink.Blink_ConvolutionP = [out.metrics.conv_Pvalue]';
        blink.Blink_PercentReduction = [out.metrics.percentreduction]';
        blink.Blink_PercentReductionP = [out.metrics.perc_Pvalue]';
        identified = out.identifiedcomponents(out.identifiedcomponents > 0);
        blink.BlinkMetricsIdentified(identified) = true;
    catch err
        blink.BlinkMetricsError(:) = string(err.message);
    end
end

function combine(runNames, outDir, perFileDir)

    perFile = dir(fullfile(perFileDir, '*.mat'));
    all = cell(numel(perFile), 1);
    for f = 1:numel(perFile)
        loaded = load(fullfile(perFileDir, perFile(f).name));
        all{f} = loaded.components;
    end
    all = vertcat(all{:});

    lines = {'# Per-component review of targeted wICA', '', ...
        sprintf('%d files. One row per component in `wICA_components_<run>.csv`.', numel(perFile)), ''};
    for r = 1:numel(runNames)
        run = all(all.Run == runNames{r}, :);
        writetable(run, fullfile(outDir, ['wICA_components_' runNames{r} '.csv']));
        lines = [lines, runSummary(run, runNames{r})]; %#ok<AGROW>
        drawSlopeFigure(run, runNames{r}, outDir);
    end
    writelines(string(lines), fullfile(outDir, 'summary.md'));
end

function lines = runSummary(run, runName)

    [files, ~, fileOfRow] = unique(run.File);
    recomputedShare = accumarray(fileOfRow, run.MuscleBySlope, [], @mean);
    storedShare = accumarray(fileOfRow, run.StoredMuscleBySlopeShare, [], @mean);
    mismatch = abs(recomputedShare - storedShare) > 1e-9;
    droppedChannels = accumarray(fileOfRow, run.ChannelsDroppedFromICA, [], @max);
    approximate = droppedChannels > 0;

    lines = {['## ' runName], '', ...
        sprintf('- Rebuild check: recomputed muscle-by-slope share matches the stored one in %d of %d files%s', ...
            sum(~mismatch), numel(files), mismatchList(files(mismatch))), ...
        sprintf(['- %d files had channels rejected after wICA (the deferred muscle screen removes them from the saved ICA); ' ...
            'their slopes and blink metrics use the saved unmixing on the remaining channels, so they are approximate, ' ...
            'and %d of them differ from the stored muscle-by-slope share'], sum(approximate), sum(approximate & mismatch)), ...
        sprintf('- icablinkmetrics ran in the run for %d of %d files; recomputed here it failed for %d files', ...
            numel(unique(run.File(run.BlinkMetricsRanInRun))), numel(files), ...
            numel(unique(run.File(run.BlinkMetricsError ~= "")))), '', ...
        '| ICLabel label (thresholded) | Components | Median slope | Muscle by slope | Blink-flagged by icablinkmetrics | Median variance share | Variance in muscle-by-slope components (sum per file, median) |', ...
        '|---|---|---|---|---|---|---|'};
    for label = unique(run.ICLabelThresholded)'
        rows = run.ICLabelThresholded == label;
        subtractedShare = accumarray(fileOfRow(rows), run.VarianceShare(rows) .* run.MuscleBySlope(rows), [numel(files) 1]);
        lines{end+1} = sprintf('| %s | %d | %.2f | %.0f%% | %d | %.4f | %.3f |', label, sum(rows), ...
            median(run.MuscleSlope(rows)), 100 * mean(run.MuscleBySlope(rows)), sum(run.BlinkMetricsIdentified(rows)), ...
            median(run.VarianceShare(rows)), median(subtractedShare)); %#ok<AGROW>
    end
    eye = run.ICLabelThresholded == "Eye";
    lines = [lines, {'', sprintf(['- Eye components: %d by ICLabel, %d by icablinkmetrics, %d by both, ' ...
        '%d icablinkmetrics-only (these would join blink-period wavelet cleaning)'], sum(eye), ...
        sum(run.BlinkMetricsIdentified), sum(eye & run.BlinkMetricsIdentified), ...
        sum(~eye & run.BlinkMetricsIdentified)), ''}];
end

function text = mismatchList(files)
    if isempty(files)
        text = '';
    else
        text = [' (mismatched: ' strjoin(cellstr(files), ', ') ')'];
    end
end

function drawSlopeFigure(run, runName, outDir)

    figure('visible', 'off', 'Position', [0 0 900 500]);
    hold on
    labels = unique(run.ICLabelThresholded)';
    edges = -4:0.1:2;
    for label = labels
        histogram(run.MuscleSlope(run.ICLabelThresholded == label), edges, 'DisplayStyle', 'stairs', 'LineWidth', 1.5);
    end
    xline(run.MuscleSlopeThreshold(1), 'k--', sprintf('threshold %.2f', run.MuscleSlopeThreshold(1)));
    legend([labels, "threshold"], 'Location', 'northwest');
    xlabel('Component log-log slope, 12-70 Hz');
    ylabel('Components');
    title(sprintf('%s: muscle slope by ICLabel label', runName));
    print(gcf, fullfile(outDir, ['slopes_by_label_' runName '.png']), '-dpng', '-r150');
    close(gcf);
end
