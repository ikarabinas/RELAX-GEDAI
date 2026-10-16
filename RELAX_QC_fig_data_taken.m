% RELAX_QC_fig_data_taken() - Figure 4: how much data each cleaning stage took, per file.

function caption = RELAX_QC_fig_data_taken(cache, outDir)

    style = RELAX_QC_style();
    shown = ismember(cache.tables.Outcome, ["cleaned", "excluded", "error"]);
    T = cache.tables(shown, :);
    S = cache.summaries(shown);
    isCleaned = T.Outcome == "cleaned";

    channelColumns = {'NChanImpedance', 'NChanPop', 'NChanBridge', 'NChanExtreme', 'NChanMuscle'};
    channels = zeroed(T{:, channelColumns});
    channels(~isCleaned, :) = 0;
    montage = median(arrayfun(@(s) s.nChannels, S), 'omitnan');
    duration = arrayfun(@(s) s.durationSeconds, S);
    excised = zeroed([T.ExtremeExcisedSeconds, T.PopExcisedSeconds]);
    excisedPercent = 100 * excised ./ (duration + sum(excised, 2));
    interpolatedPercent = 100 * T.PopFractionSamplesInterpolated;
    variancePercent = 100 * arrayfun(@(s) s.varianceRemovedFraction, S);
    icPercent = 100 * T.WicaArtifactICsReduced;

    lost = sum(channels, 2) / montage + sum(excisedPercent, 2, 'omitnan') / 100;
    lost(~isCleaned) = Inf;
    [~, order] = sort(lost);
    [T, channels, excisedPercent, interpolatedPercent, variancePercent, icPercent, isCleaned] = deal( ...
        T(order, :), channels(order, :), excisedPercent(order, :), interpolatedPercent(order), ...
        variancePercent(order), icPercent(order), isCleaned(order));
    x = 1:height(T);

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 max(18, 0.35 * numel(x) + 6) 24], 'Color', 'w');
    layout = tiledlayout(fig, 4, 1, 'TileSpacing', 'loose', 'Padding', 'compact');
    if isempty(x)
        RELAX_QC_draw_note(nexttile(layout, [4 1]), "No processed files");
    else
        ax = nexttile(layout);
        stack(ax, x, channels, style.channelScreens, {'Impedance', 'Pop', 'Bridge', 'Extreme data', 'Muscle'});
        if ~isnan(montage)
            yline(ax, cache.meta.maxChannelProportion * montage, '--', 'Rejection budget', 'Color', style.muted, ...
                  'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'FontSize', style.fontSize - 1);
        end
        ylabel(ax, {'Channels', 'rejected'});
        markNotCleaned(ax, x(~isCleaned), T.Outcome(~isCleaned), style);

        ax = nexttile(layout);
        stack(ax, x, excisedPercent, style.channelScreens([4 2], :), {'Extreme periods', 'Pop excision'});
        ylabel(ax, {'% of recording', 'excised'});

        ax = nexttile(layout);
        bar(ax, x, interpolatedPercent, 0.8, 'FaceColor', style.channelScreens(2, :), 'EdgeColor', 'none');
        ylabel(ax, {'% samples filled', 'by pop repair'});

        ax = nexttile(layout); hold(ax, 'on');
        bar(ax, x, variancePercent, 0.8, 'FaceColor', style.removed, 'EdgeColor', 'none', 'DisplayName', 'Variance removed by wICA');
        plot(ax, x, icPercent, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', style.ink, 'MarkerFaceColor', 'w', ...
             'DisplayName', 'Artifact ICs reduced');
        ylabel(ax, '%'); legend(ax, 'Location', 'northoutside', 'NumColumns', 2, 'Box', 'off');
        xticks(ax, x); xticklabels(ax, T.File); ax.TickLabelInterpreter = 'none'; xtickangle(ax, 90);

        axesList = findobj(fig, 'Type', 'axes');
        set(axesList, 'XLim', [0.4, numel(x) + 0.6], 'FontSize', style.fontSize, 'YGrid', 'on');
        set(axesList(axesList ~= ax), 'XTickLabel', []);
    end
    RELAX_QC_save_figure(fig, outDir, 'fig4_data_taken');

    caption = string(sprintf(['## Figure 4 — Data taken by cleaning, per file\n\n' ...
        '%d files shown (%d cleaned, %d excluded, %d errored), sorted by the fraction of channels rejected plus the ' ...
        'fraction of time excised; files not yet processed are left out. Panels: channels rejected by each screen ' ...
        'against the shared rejection budget; %% of the recording excised; %% of samples interpolated by pop repair; ' ...
        'variance removed by wICA (baseline − cleaned over baseline) with the share of artifact ICs it reduced. ' ...
        'Look for files at the right-hand end that lost far more than the rest.\n'], ...
        height(T), sum(isCleaned), sum(T.Outcome == "excluded"), sum(T.Outcome == "error")));
end

function x = zeroed(x)
    x(isnan(x)) = 0;
end

function stack(ax, x, values, colours, names)
    b = bar(ax, x, values, 0.8, 'stacked', 'EdgeColor', 'w', 'LineWidth', 0.5);
    for k = 1:numel(b)
        b(k).FaceColor = colours(k, :);
        b(k).DisplayName = names{k};
    end
    legend(ax, b, 'Location', 'northoutside', 'NumColumns', numel(b), 'Box', 'off', 'AutoUpdate', 'off');
end

function markNotCleaned(ax, x, outcome, style)
    for k = 1:numel(x)
        text(ax, x(k), 0, " " + outcome(k), 'Rotation', 90, 'Color', style.muted, ...
             'FontSize', style.fontSize - 1, 'VerticalAlignment', 'middle');
    end
end
