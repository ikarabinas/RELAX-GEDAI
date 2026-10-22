% RELAX_QC_fig_blink_topo() - Figure 5: scalp maps at the blink peak before and after
%                             wICA, for one example file and averaged over the dataset,
%                             with the blink-locked time course.

function caption = RELAX_QC_fig_blink_topo(cache, outDir)

    style = RELAX_QC_style();
    ok = cache.summaries([cache.summaries.status] == "ok");
    eligible = ok([ok.blinkN] >= cache.meta.minBlinksForGroup);
    tooFew = numel(ok) - numel(eligible);

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 30 12], 'Color', 'w');
    layout = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'loose');

    [example, exampleFile] = blinkExample(cache, ok);
    if isempty(example)
        RELAX_QC_draw_note(nexttile(layout, 1, [1 3]), "No file with blinks");
    else
        mapRow(layout, 1, {example.topoBaseline, example.topoCleaned, example.topoRemoved}, example.chanlocs, ...
               exampleFile, style);
    end

    if isempty(eligible)
        RELAX_QC_draw_note(nexttile(layout, 5, [1 3]), "No file has enough blinks for a group average");
        RELAX_QC_draw_note(nexttile(layout, 4, [2 1]), "");
    else
        reference = eligible(1).chanlocs;
        mapRow(layout, 2, {groupMean(eligible, 'topoBaseline', reference), groupMean(eligible, 'topoCleaned', reference), ...
                           groupMean(eligible, 'topoRemoved', reference)}, reference, "mean of " + numel(eligible) + " files", style);

        ax = nexttile(layout, 4, [2 1]); hold(ax, 'on');
        ms = 1000 * eligible(1).blinkTimes;
        timeCourse(ax, ms, vertcat(eligible.blinkErpBaseline), style.before, 'Before wICA', style);
        timeCourse(ax, ms, vertcat(eligible.blinkErpCleaned), style.cleaned, 'Cleaned', style);
        xlabel(ax, 'Time from blink peak (ms)'); ylabel(ax, 'Blink electrodes (\muV)');
        yRange = ylim(ax); ylim(ax, [yRange(1), yRange(2) + 0.3 * diff(yRange)]);
        legend(ax, 'Location', 'northeast', 'Box', 'off'); ax.YGrid = 'on'; ax.FontSize = style.fontSize;
        title(ax, 'Median and IQR', 'FontWeight', 'normal');
    end
    RELAX_QC_save_figure(fig, outDir, 'fig5_blink_topography');

    caption = string(sprintf(['## Figure 5 — Blink-locked topography\n\n' ...
        'Scalp maps of the average over ±50 ms of the pipeline''s EyeBlinkMax markers, each row on one colour scale. ' ...
        'Blink averages and maps are baseline-corrected to the mean of −500 to −300 ms before the blink peak. ' ...
        'Top: %s. Bottom: mean of %d files, each file weighted equally; %d file%s with fewer than %d blinks left out. ' ...
        'Right: blink-electrode average, median and interquartile range across those files. ' ...
        'Good cleaning: the frontal blink field in "Before wICA" moves into "Removed", and "Cleaned" is flat with ' ...
        'no residual frontal peak.\n'], exampleDescription(example, exampleFile), numel(eligible), tooFew, ...
        plural(tooFew), cache.meta.minBlinksForGroup));
end

function s = plural(n)
    s = 's';
    if n == 1, s = ''; end
end

function [example, file] = blinkExample(cache, ok)
    example = []; file = "";
    preferred = [];
    if ~isempty(cache.picks)
        preferred = find(contains([cache.picks.role], "worst blink"), 1);
        if isempty(preferred), preferred = 1; end
    end
    if ~isempty(preferred)
        i = find(cache.tables.File == cache.picks(preferred).file, 1);
        s = cache.summaries(i);
        if s.status == "ok" && s.blinkN > 0
            example = s; file = cache.tables.File(i); return
        end
    end
    withBlinks = find([ok.blinkN] > 0, 1);
    if ~isempty(withBlinks)
        example = ok(withBlinks);
        okFiles = cache.tables.File([cache.summaries.status] == "ok");
        file = okFiles(withBlinks);
    end
end

function d = exampleDescription(example, file)
    if isempty(example), d = "no file with blinks"; else, d = file + " (" + example.blinkN + " blinks)"; end
end

function m = groupMean(summaries, field, reference)
    labels = {reference.labels};
    stacked = nan(numel(labels), numel(summaries));
    for k = 1:numel(summaries)
        [present, where] = ismember(labels, {summaries(k).chanlocs.labels});
        stacked(present, k) = summaries(k).(field)(where(present));
    end
    m = mean(stacked, 2, 'omitnan');
end

function mapRow(layout, row, maps, chanlocs, rowLabel, style)
    titles = ["Before wICA", "Cleaned", "Removed by wICA + GEDAI"];
    limit = max(abs([maps{:}]), [], 'all');
    if ~(limit > 0), limit = 1; end
    fig = ancestor(layout, 'figure');
    slots = arrayfun(@(j) nexttile(layout, (row - 1) * 4 + j), 1:3);
    drawnow;   % tiled axes only have their final position after a draw; topoplot cannot reposition a tile
    for j = 1:3
        slot = slots(j);
        position = slot.InnerPosition .* [1 1 1 0.85];
        position(3) = 0.78 * position(3);
        ax = axes('Parent', fig, 'Units', slot.Units, 'Position', position);
        axis(slot, 'off');
        topoplot(maps{j}, chanlocs, 'maplimits', [-limit limit], 'electrodes', 'off', 'colormap', style.diverging);
        title(ax, {char(titles(j)), char(rowLabel)}, 'Interpreter', 'none', 'FontSize', style.fontSize, 'FontWeight', 'normal');
    end
    barPosition = [position(1) + 1.04 * position(3), position(2) + 0.2 * position(4), 0.1 * position(3), 0.6 * position(4)];
    scale = colorbar(ax, 'Position', barPosition); scale.Title.String = '\muV';
    colormap(ax, style.diverging);
    clim(ax, [-limit limit]);
end

function timeCourse(ax, ms, erps, colour, name, style)
    q = quantile(erps, [0.25 0.5 0.75], 1);
    if size(erps, 1) == 1, q = [erps; erps; erps]; end
    fill(ax, [ms fliplr(ms)], [q(1, :) fliplr(q(3, :))], colour, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax, ms, q(2, :), 'Color', colour, 'LineWidth', style.bold, 'DisplayName', name);
end
