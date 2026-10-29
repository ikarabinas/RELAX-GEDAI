% compare_wICA_runs() - Paired comparison of what targeted wICA removed in two or more runs of
%                       RELAX_Wrapper_slurm.m over the same files.
%
% Usage:
%   >> summary = compare_wICA_runs(runDirs, runNames, outDir);
%
% Inputs:
%   runDirs  - cellstr of RELAXProcessed output folders, each run with
%              Compute_wICA_removal_metrics=1; the first is the baseline the others are compared to
%   runNames - cellstr of short names for the runs, same order
%   outDir   - folder for wICA_comparison.csv, summary.md and the figures
%
% The yardstick is RELAX_wICA_removal_metrics: the share of input power wICA removed in epochs
% without muscle (brain signal lost, lower is better, chiefly alpha and beta) and in muscle epochs
% (artifact taken, high gamma); blink epochs are reported separately. Only files every run processed are compared. Each comparison row
% gives the per-run medians, the median within-file difference from the baseline, the number of
% files where the run removed less than the baseline, and a Wilcoxon signed-rank p.
%
% Output:
%   summary - table, one row per metric and compared run

function summary = compare_wICA_runs(runDirs, runNames, outDir)

    if ~exist(outDir, 'dir'); mkdir(outDir); end

    metricNames = {'NumberOfICs', 'Proportion_ICs_muscle_by_slope', 'ProportionICs_was_Muscle', ...
        'ProportionICs_was_Eye', 'Proportion_power_removed_by_muscle_subtraction', ...
        'ProportionChannelEpochsMuscle', 'ProportionEpochsBlink'};
    for kind = {'Clean', 'Muscle', 'Blink'}
        for band = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'HighGamma', 'Broadband'}
            metricNames{end+1} = ['RemovedFraction_' kind{1} '_' band{1}]; %#ok<AGROW>
        end
    end

    runs = cellfun(@loadRunTable, runDirs, 'UniformOutput', false);
    files = runs{1}.aFileName;
    for r = 2:numel(runs)
        files = intersect(files, runs{r}.aFileName, 'stable');
    end
    fprintf('%d files processed by every run\n', numel(files));

    values = nan(numel(files), numel(metricNames), numel(runs));
    for r = 1:numel(runs)
        [~, rows] = ismember(files, runs{r}.aFileName);
        for m = 1:numel(metricNames)
            if ismember(metricNames{m}, runs{r}.Properties.VariableNames)
                values(:, m, r) = double(runs{r}.(metricNames{m})(rows));
            end
        end
    end

    perFile = table(files, 'VariableNames', {'File'});
    for r = 1:numel(runs)
        for m = 1:numel(metricNames)
            perFile.([runNames{r} '_' metricNames{m}]) = values(:, m, r);
        end
    end
    writetable(perFile, fullfile(outDir, 'wICA_comparison_per_file.csv'));

    summary = table();
    for r = 2:numel(runs)
        for m = 1:numel(metricNames)
            baseline = values(:, m, 1);
            other = values(:, m, r);
            paired = ~isnan(baseline) & ~isnan(other);
            p = NaN;
            if nnz(paired) >= 5 && any(other(paired) ~= baseline(paired))
                p = signrank(other(paired), baseline(paired));
            end
            summary = [summary; table(string(runNames{r}), string(metricNames{m}), nnz(paired), ...
                median(baseline(paired)), median(other(paired)), median(other(paired) - baseline(paired)), ...
                nnz(other(paired) < baseline(paired)), p, ...
                'VariableNames', {'Run', 'Metric', 'Files', 'BaselineMedian', 'RunMedian', ...
                'MedianDifference', 'FilesWhereRunLower', 'SignrankP'})]; %#ok<AGROW>
        end
    end
    writetable(summary, fullfile(outDir, 'wICA_comparison.csv'));

    drawPairedScatter(values, metricNames, runNames, outDir);
    drawTopographies(runDirs, runNames, files, outDir);
    writeSummaryMarkdown(summary, runDirs, runNames, numel(files), outDir);
end

function runTable = loadRunTable(runDir)
    statsFiles = dir(fullfile(runDir, 'PerParticipantStats', 'ProcessingStatistics_wICA_*.mat'));
    runTable = [];
    for f = 1:numel(statsFiles)
        s = load(fullfile(statsFiles(f).folder, statsFiles(f).name), 'RELAXProcessing_wICA_AllParticipants');
        rows = s.RELAXProcessing_wICA_AllParticipants;
        rows = rows(~cellfun(@isempty, rows.aFileName), :);
        for k = 1:height(rows)
            runTable = RELAX_append_table_row(runTable, height(runTable) + 1, rows(k, :));
        end
    end
    runTable.aFileName = string(runTable.aFileName);
end

function drawPairedScatter(values, metricNames, runNames, outDir)
    shown = {'RemovedFraction_Clean_Alpha', 'RemovedFraction_Clean_Beta', 'RemovedFraction_Clean_Broadband', ...
        'RemovedFraction_Muscle_HighGamma', 'Proportion_ICs_muscle_by_slope', 'NumberOfICs'};
    for r = 2:size(values, 3)
        figure('Visible', 'off', 'Position', [0 0 1200 750]);
        tiledlayout(2, 3, 'TileSpacing', 'compact');
        for k = 1:numel(shown)
            m = strcmp(metricNames, shown{k});
            nexttile;
            x = values(:, m, 1); y = values(:, m, r);
            scatter(x, y, 18, 'filled', 'MarkerFaceAlpha', 0.6); hold on
            limits = [min([x; y], [], 'omitnan') max([x; y], [], 'omitnan')];
            if all(isfinite(limits)) && diff(limits) > 0
                plot(limits, limits, 'k:'); xlim(limits); ylim(limits);
            end
            axis square
            xlabel(runNames{1}, 'Interpreter', 'none'); ylabel(runNames{r}, 'Interpreter', 'none');
            title(strrep(shown{k}, '_', ' '));
        end
        exportgraphics(gcf, fullfile(outDir, sprintf('paired_%s_vs_%s.png', runNames{r}, runNames{1})), 'Resolution', 150);
        close(gcf);
    end
end

function drawTopographies(runDirs, runNames, files, outDir)
% Mean over files of each channel's removed share of alpha in its clean (no muscle, no blink) epochs.
    alphaBand = 3;
    figure('Visible', 'off', 'Position', [0 0 450 * numel(runDirs) 420]);
    tiledlayout(1, numel(runDirs), 'TileSpacing', 'compact');
    for r = 1:numel(runDirs)
        sums = containers.Map('KeyType', 'char', 'ValueType', 'double');
        counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
        locs = containers.Map('KeyType', 'char', 'ValueType', 'any');
        for f = 1:numel(files)
            removalFile = fullfile(runDirs{r}, 'wICA_removal', files(f) + "_wICA_removal.mat");
            if ~isfile(removalFile); continue; end
            s = load(removalFile, 'removalPerChannel');
            pc = s.removalPerChannel;
            share = pc.removedPowerClean(:, alphaBand) ./ pc.inputPowerClean(:, alphaBand);
            for c = find(pc.channelEpochsClean(:)' > 0)
                label = pc.labels{c};
                if ~isKey(sums, label); sums(label) = 0; counts(label) = 0; locs(label) = pc.chanlocs(c); end
                sums(label) = sums(label) + share(c);
                counts(label) = counts(label) + 1;
            end
        end
        nexttile;
        if sums.Count == 0; title(runNames{r}); continue; end
        labels = keys(sums);
        chanlocs = cellfun(@(l) locs(l), labels);
        meanShare = cellfun(@(l) sums(l) / counts(l), labels);
        topoplot(meanShare, chanlocs, 'electrodes', 'off', 'maplimits', [0 max(0.05, max(meanShare))]);
        colorbar;
        title(sprintf('%s: alpha removed, clean epochs', runNames{r}), 'Interpreter', 'none');
    end
    exportgraphics(gcf, fullfile(outDir, 'topography_clean_alpha_removed.png'), 'Resolution', 150);
    close(gcf);
end

function writeSummaryMarkdown(summary, runDirs, runNames, nFiles, outDir)
    fid = fopen(fullfile(outDir, 'summary.md'), 'w');
    closeFile = onCleanup(@() fclose(fid));
    fprintf(fid, '# wICA removal comparison\n\n');
    for r = 1:numel(runDirs)
        fprintf(fid, '- `%s`: %s%s\n', runNames{r}, runDirs{r}, repmat(' (baseline)', 1, r == 1));
    end
    fprintf(fid, '\n%d files processed by every run. RemovedFraction = power wICA removed / power going in, ', nFiles);
    fprintf(fid, 'pooled over channels and 1 s epochs. Blink = epochs touching a blink; in the rest each channel is ');
    fprintf(fid, 'Muscle where its own 7-75 Hz slope reaches the muscle threshold, else Clean (brain loss, lower is better).\n\n');
    fprintf(fid, '| Run | Metric | Files | Baseline median | Run median | Median diff | Files run lower | Signrank p |\n');
    fprintf(fid, '|---|---|---|---|---|---|---|---|\n');
    for k = 1:height(summary)
        fprintf(fid, '| %s | %s | %d | %.4g | %.4g | %.4g | %d | %.3g |\n', summary.Run(k), summary.Metric(k), ...
            summary.Files(k), summary.BaselineMedian(k), summary.RunMedian(k), summary.MedianDifference(k), ...
            summary.FilesWhereRunLower(k), summary.SignrankP(k));
    end
end
