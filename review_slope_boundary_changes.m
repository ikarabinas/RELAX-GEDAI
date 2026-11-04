% review_slope_boundary_changes() - Components whose muscle call changes when the slope is fitted
%                                    within each segment instead of over the whole record, for one
%                                    finished wICA run file: a table, and spectra and topographies.
%
% Usage:
%   >> review_slope_boundary_changes(runDir, cleanedName, outDir);
%
% Inputs:
%   runDir      - RELAXProcessed output folder of a targeted wICA run
%   cleanedName - file name in runDir/Cleaned_Data
%   outDir      - folder for components.csv and components_page<N>.png
%
% Each plotted row is one component whose call changes: its topography, then its binned
% spectrum fitted over the whole record and within segments, over 12-70 Hz excluding line noise.
% The wICA input is rebuilt as review_wICA_components.m does.

function review_slope_boundary_changes(runDir, cleanedName, outDir)

    rowsPerPage = 6;
    if ~exist(outDir, 'dir'); mkdir(outDir); end
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
    mixing = pinv(unmixing);
    activations = unmixing * double(wICAinput.data);

    segments = RELAX_get_segments(wICAinput);
    [wholeRecord, binsWhole, freq] = RELAX_IC_muscle_slopes(activations, wICAinput.srate, RELAX_cfg.LineNoiseFrequency);
    [bySegment, binsSegment] = RELAX_IC_muscle_slopes(activations, wICAinput.srate, RELAX_cfg.LineNoiseFrequency, segments);
    threshold = RELAX_cfg.MuscleSlopeThreshold;
    gained = wholeRecord < threshold & bySegment >= threshold;
    lost = wholeRecord >= threshold & bySegment < threshold;
    changed = find(gained | lost);

    probabilities = cleaned.etc.ic_classification.ICLabel.classifications;
    classes = cleaned.etc.ic_classification.ICLabel.classes;
    [topProbability, topClass] = max(probabilities, [], 2);
    inputPower = sum(double(wICAinput.data(:)).^2);
    varianceShare = sum(mixing.^2, 1)' .* sum(activations.^2, 2) / inputPower;
    segmentLengths = segments(:,2) - segments(:,1) + 1;

    components = table(changed', repmat({''}, numel(changed), 1), wholeRecord(changed)', bySegment(changed)', ...
        classes(topClass(changed))', topProbability(changed), varianceShare(changed), ...
        'VariableNames', {'Component','Change','SlopeWholeRecord','SlopeBySegment','ICLabelClass','ICLabelProbability','VarianceShare'});
    components.Change(gained(changed)) = {'gained'};
    components.Change(lost(changed)) = {'lost'};
    writetable(components, fullfile(outDir, 'components.csv'));
    fprintf('%s: %d boundaries, segments %.1f-%.1f s (median %.1f s), %d components, %d gained, %d lost\n', ...
        fileName, size(segments,1) - 1, min(segmentLengths)/wICAinput.srate, max(segmentLengths)/wICAinput.srate, ...
        median(segmentLengths)/wICAinput.srate, size(activations,1), sum(gained), sum(lost));
    disp(components)

    fitBand = freq >= 12 & freq <= 70 & (freq < RELAX_cfg.LineNoiseFrequency - 6 | freq > RELAX_cfg.LineNoiseFrequency + 6);
    for page = 1:ceil(numel(changed) / rowsPerPage)
        figure('Visible', 'off', 'Position', [0 0 900 160 * rowsPerPage]);
        for r = 1:min(rowsPerPage, numel(changed) - (page - 1) * rowsPerPage)
            c = changed((page - 1) * rowsPerPage + r);
            subplot(rowsPerPage, 4, (r - 1) * 4 + 1);
            topoplot(mixing(:, c), wICAinput.chanlocs, 'electrodes', 'off', 'numcontour', 0);
            title(sprintf('IC%d %s (%s %.2f)', c, components.Change{components.Component == c}, classes{topClass(c)}, topProbability(c)), 'FontSize', 8);
            subplot(rowsPerPage, 4, (r - 1) * 4 + (2:4));
            loglog(freq, binsWhole(c,:), 'Color', [0.6 0.6 0.6]); hold on
            loglog(freq, binsSegment(c,:), 'k');
            for fitted = {{binsWhole, wholeRecord, [0.6 0.6 0.6]}, {binsSegment, bySegment, 'k'}}
                bins = fitted{1}{1}; slope = fitted{1}{2}(c);
                p = polyfit(log(freq(fitBand)), log(bins(c, fitBand)), 1);
                loglog(freq(fitBand), exp(polyval([slope p(2)], log(freq(fitBand)))), '--', 'Color', fitted{1}{3}, 'LineWidth', 1.5);
            end
            xlim([1 100]);
            title(sprintf('slope whole %.2f (grey), by segment %.2f (black); threshold %.2f', wholeRecord(c), bySegment(c), threshold), 'FontSize', 8);
        end
        exportgraphics(gcf, fullfile(outDir, sprintf('components_page%d.png', page)), 'Resolution', 110);
        close(gcf);
    end
end
