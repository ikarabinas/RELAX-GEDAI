% RELAX_QC_fig_data_taken() - Figure 4: how much data each cleaning stage took, per file.

function caption = RELAX_QC_fig_data_taken(cache, outDir)

    style = RELAX_QC_style();
    shown = ismember(cache.tables.Outcome, ["cleaned", "excluded", "error"]);
    T = cache.tables(shown, :);
    S = cache.summaries(shown);
    isCleaned = T.Outcome == "cleaned";

    channelColumns = {'NChanImpedance', 'NChanPop', 'NChanBridge', 'NChanExtreme', 'NChanMuscle'};
    requiredColumns = {'NChanImpedance', 'NChanPop', 'NChanExtreme', 'ExtremeExcisedSeconds', 'PopFractionSamplesInterpolated', ...
                       'WicaArtifactICsReduced'};
    isSkipped = string({S.status})' == "skipped";
    lacksCovariates = isCleaned & ~isSkipped & any(isnan(T{:, requiredColumns}), 2);
    channels = zeroed(T{:, channelColumns});
    channels(~isCleaned, :) = 0;
    montage = median(arrayfun(@(s) s.nChannels, S), 'omitnan');
    duration = arrayfun(@(s) s.durationSeconds, S);
    excised = zeroed([T.ExtremeExcisedSeconds, T.PopExcisedSeconds]);
    excisedPercent = 100 * excised ./ (duration + sum(excised, 2));
    interpolatedPer10k = 1e4 * T.PopFractionSamplesInterpolated;
    variancePercent = 100 * arrayfun(@(s) s.varianceRemovedFraction, S);
    icPercent = 100 * T.WicaArtifactICsReduced;
    interpolatedPer10k(isSkipped) = NaN;
    icPercent(isSkipped) = NaN;

    lost = sum(channels, 2) / montage + sum(excisedPercent, 2, 'omitnan') / 100;
    lost(~isCleaned) = Inf;
    [~, order] = sort(lost);
    [T, channels, excisedPercent, interpolatedPer10k, variancePercent, icPercent, isCleaned, isSkipped, lacksCovariates] = deal( ...
        T(order, :), channels(order, :), excisedPercent(order, :), interpolatedPer10k(order), ...
        variancePercent(order), icPercent(order), isCleaned(order), isSkipped(order), lacksCovariates(order));
    x = 1:height(T);
    budget = floor(cache.meta.maxChannelProportion * montage);
    suffix = commonSuffix(T.File);
    labels = extractBefore(T.File, strlength(T.File) - strlength(suffix) + 1);

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 max(18, 0.35 * numel(x) + 6) 24], 'Color', 'w');
    layout = tiledlayout(fig, 4, 1, 'TileSpacing', 'loose', 'Padding', 'compact');
    if isempty(x)
        RELAX_QC_draw_note(nexttile(layout, [4 1]), "No processed files");
    else
        ax = nexttile(layout);
        stack(ax, x, channels, style.channelScreens, {'Impedance', 'Pop', 'Bridge', 'Extreme data', 'Muscle'});
        limit = max(5, 1.2 * max(sum(channels, 2)));
        ylim(ax, [0, limit]); yticks(ax, 0:max(1, ceil(limit / 6)):limit);
        if budget <= limit
            yline(ax, budget, '--', 'Rejection budget', 'Color', style.muted, ...
                  'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'FontSize', style.fontSize - 1);
        end
        ylabel(ax, {'Channels', 'rejected'});
        markNotCleaned(ax, x(~isCleaned), T.Outcome(~isCleaned), style);
        markNotCleaned(ax, x(lacksCovariates), repmat("n/a", nnz(lacksCovariates), 1), style);

        ax = nexttile(layout);
        stack(ax, x, excisedPercent, style.channelScreens([4 2], :), {'Extreme periods', 'Pop excision'});
        ylabel(ax, {'% of recording', 'excised'});
        headroom(ax, sum(excisedPercent, 2, 'omitnan'));
        markNotCleaned(ax, x(isSkipped), repmat("skipped", nnz(isSkipped), 1), style);

        ax = nexttile(layout);
        bar(ax, x, interpolatedPer10k, 0.8, 'FaceColor', style.channelScreens(2, :), 'EdgeColor', 'none');
        ax.YAxis.Exponent = 0;
        headroom(ax, interpolatedPer10k);
        ylabel(ax, {'Samples filled by', 'pop repair (per 10^4)'});

        ax = nexttile(layout); hold(ax, 'on');
        bar(ax, x, variancePercent, 0.8, 'FaceColor', style.removed, 'EdgeColor', 'none', 'DisplayName', 'Variance removed by wICA + GEDAI');
        plot(ax, x, icPercent, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', style.ink, 'MarkerFaceColor', 'w', ...
             'DisplayName', 'Artifact ICs (% of all ICs)');
        markNotCleaned(ax, x(isSkipped), repmat("skipped", nnz(isSkipped), 1), style);
        ylabel(ax, '%'); legend(ax, 'Location', 'northoutside', 'NumColumns', 2, 'Box', 'off');
        xticks(ax, x); xticklabels(ax, labels); ax.TickLabelInterpreter = 'none'; xtickangle(ax, 90); ax.FontSize = style.fontSize - 1;

        axesList = findobj(fig, 'Type', 'axes');
        set(axesList, 'XLim', [0.4, numel(x) + 0.6], 'YGrid', 'on');
        set(axesList(axesList ~= ax), 'FontSize', style.fontSize);
        set(axesList(axesList ~= ax), 'XTickLabel', []);
    end
    RELAX_QC_save_figure(fig, outDir, 'fig4_data_taken');

    caption = string(sprintf(['## Figure 4 — Data taken by cleaning, per file\n\n' ...
        '%d files shown (%d cleaned, %d skipped, %d excluded, %d errored), sorted by the fraction of channels rejected plus the ' ...
        'fraction of time excised; files not yet processed are left out. Panels: channels rejected by each screen ' ...
        'against the shared rejection budget; %% of the recording excised; samples filled by pop repair per 10,000 samples; ' ...
        'variance removed by wICA + GEDAI, where on (baseline − cleaned over baseline), with markers for the share of all ICs that ICLabel ' ...
        'classed as artifact and wICA reduced. ' ...
        'Look for files at the right-hand end that lost far more than the rest. ' ...
        'The shared rejection budget is %s. Missing values are drawn as zero-height bars: %d cleaned file%s had no covariate ' ...
        'row (marked "n/a"), and %d file%s skipped for alignment: excision and variance bars are absent (marked "skipped").%s\n'], ...
        height(T), sum(isCleaned & ~isSkipped), nnz(isSkipped), sum(T.Outcome == "excluded"), sum(T.Outcome == "error"), ...
        budgetText(cache.meta.maxChannelProportion, budget, montage), nnz(lacksCovariates), pluralSuffix(nnz(lacksCovariates)), ...
        nnz(isSkipped), pluralSuffix(nnz(isSkipped)), suffixNote(suffix)));
end

function s = pluralSuffix(n)
    s = 's';
    if n == 1, s = ''; end
end

function headroom(ax, values)
    top = max(values, [], 'omitnan');
    if top > 0, ylim(ax, [0, 1.1 * top]); end
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

function suffix = commonSuffix(names)
    suffix = "";
    if numel(names) < 2, return, end
    shortest = min(strlength(names));
    for k = 1:shortest - 1
        candidate = extractAfter(names(1), strlength(names(1)) - k);
        if ~all(endsWith(names, candidate)), break, end
        suffix = candidate;
    end
end

function note = suffixNote(suffix)
    note = "";
    if suffix ~= "", note = sprintf(' File labels omit the shared suffix "%s".', suffix); end
end

function text = budgetText(proportion, budget, montage)
    text = sprintf('%.0f%% of the montage', 100 * proportion);
    if ~isnan(montage), text = sprintf('%s (%.0f of %.0f channels)', text, budget, montage); end
end
