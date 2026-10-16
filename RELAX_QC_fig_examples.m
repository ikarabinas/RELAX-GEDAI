% RELAX_QC_fig_examples() - Figure 1: before/after traces at the largest blink, the
%                           strongest muscle burst and a repaired pop, one figure per
%                           example file.

function caption = RELAX_QC_fig_examples(cache, outDir)

    style = RELAX_QC_style();
    titles = struct('blink', "Largest blink", 'muscle', "Strongest muscle burst", 'pop', "Repaired pop");
    beforeNames = struct('blink', "Before wICA", 'muscle', "Before wICA", 'pop', "No pop repair (reconstruction)");
    afterNames  = struct('blink', "Cleaned", 'muscle', "Cleaned", 'pop', "Pop-repaired");

    caption = string(sprintf(['## Figure 1 — Example traces\n\n' ...
        'One figure per example file (%d): the median file and the worst files by blink amplitude, muscle power and ' ...
        'pop count, unless hand-picked. Each 6 s window shows the blink electrodes and the two channels with the ' ...
        'most 20–80 Hz power. Blink and muscle columns: grey is the sample-matched data before wICA, blue is the ' ...
        'cleaned file, and the strip below is what wICA removed. Pop column: grey is the input file band-passed and ' ...
        'resampled without any pop repair, a display-only reconstruction rather than saved pipeline data, beside ' ...
        'the pop-repaired data; the strip below shows the repaired channel only, because the other channels differ only by the processing choices of the reconstruction; the repaired interval is shaded, widened to at least 50 ms so it can be seen. ' ...
        'The reconstruction also band-stops the line harmonics rather than applying ZaplinePlus, so line-noise differences there are not pop repair. ' ...
        'Traces and the pop removed strip are clipped to fixed limits, the pop removed strip at 1.2 times its 99th-percentile absolute value (marked "clipped"). Dashed black lines are boundaries. Good cleaning: the ' ...
        'blink and burst are gone from blue while the background rhythm is unchanged, and the removed strip ' ...
        'carries only the artifact.\n'], numel(cache.examples)));

    if isempty(cache.examples)
        fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 18 6], 'Color', 'w');
        RELAX_QC_draw_note(axes(fig), "No example files: no aligned file had blink, muscle or pop scores");
        RELAX_QC_save_figure(fig, outDir, 'fig1_examples_none');
    end

    for k = 1:numel(cache.examples)
        example = cache.examples(k);
        fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 28 15], 'Color', 'w');
        layout = tiledlayout(fig, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(layout, example.file + " (" + example.role + ")", 'Interpreter', 'none', 'FontSize', style.fontSize + 1);
        for j = 1:numel(example.windows)
            w = example.windows(j);
            kind = char(w.kind);
            traceAx = nexttile(layout, j, [3 1]);
            removedAx = nexttile(layout, 9 + j);
            if ~w.available
                RELAX_QC_draw_note(traceAx, w.note);
                axis(removedAx, 'off');
                title(traceAx, titles.(kind), 'Visible', 'on');
                continue
            end
            drawWindow(traceAx, removedAx, w, beforeNames.(kind), afterNames.(kind), style);
            title(traceAx, titles.(kind));
        end
        RELAX_QC_save_figure(fig, outDir, sprintf('fig1_examples_%d_%s', k, example.file));
    end
end

function seconds = minimumShadeSeconds(), seconds = 0.05; end

function label = removedLabel(w, isPop)
    label = 'Removed (\muV)';
    if isPop, label = [w.channels{1} ' removed (\muV)']; end
end

function drawWindow(traceAx, removedAx, w, beforeName, afterName, style)
    before = w.before; after = w.after;
    centre = median(before, 2);
    spacing = max(8 * median(mad(before, 1, 2)), eps);
    offsets = (size(before, 1) - 1:-1:0)' * spacing;
    yRange = [-1.5 * spacing, max(offsets) + 3 * spacing];

    hold(traceAx, 'on');
    shade(traceAx, w.shade, yRange, style);
    hBefore = plot(traceAx, w.t, before - centre + offsets, 'Color', style.before, 'LineWidth', style.thin);
    hAfter  = plot(traceAx, w.t, after - centre + offsets, 'Color', style.cleaned, 'LineWidth', style.thin);
    for b = w.boundaries, xline(traceAx, b, '--', 'Color', style.ink, 'LineWidth', 0.8, 'HandleVisibility', 'off'); end
    [ticks, order] = sort(offsets);
    yticks(traceAx, ticks); yticklabels(traceAx, w.channels(order));
    ylim(traceAx, yRange); xlim(traceAx, [w.t(1) w.t(end)]);
    traceAx.XTickLabel = []; traceAx.FontSize = style.fontSize; traceAx.TickLabelInterpreter = 'none';
    legend(traceAx, [hBefore(1) hAfter(1)], {beforeName, afterName}, 'Location', 'northoutside', ...
           'Orientation', 'horizontal', 'Box', 'off', 'FontSize', style.fontSize - 1);
    traceStacked = [before - centre + offsets; after - centre + offsets];
    markClipped(traceAx, any(traceStacked(:) < yRange(1) | traceStacked(:) > yRange(2)), style);

    hold(removedAx, 'on');
    isPop = w.kind == "pop";
    [removedTraces, removedRange, removedClipped] = RELAX_QC_removed_strip(w.removed, isPop);
    shade(removedAx, w.shade, removedRange, style);
    plot(removedAx, w.t, removedTraces', 'Color', [style.removed 0.6], 'LineWidth', style.thin);
    for b = w.boundaries, xline(removedAx, b, '--', 'Color', style.ink, 'LineWidth', 0.8); end
    xlim(removedAx, [w.t(1) w.t(end)]); ylim(removedAx, removedRange);
    markClipped(removedAx, removedClipped, style);
    xlabel(removedAx, 'Time (s)'); ylabel(removedAx, removedLabel(w, isPop)); removedAx.FontSize = style.fontSize;
end

function markClipped(ax, isClipped, style)
    if ~isClipped, return, end
    text(ax, 0.99, 0.02, 'clipped', 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'bottom', 'Color', style.muted, 'FontSize', style.fontSize - 1, 'FontAngle', 'italic');
end

function shade(ax, interval, yRange, style)
    if isempty(interval), return, end
    centre = mean(interval);
    interval = centre + [-1 1] * max(diff(interval), minimumShadeSeconds()) / 2;
    patch(ax, interval([1 2 2 1]), yRange([1 1 2 2]), style.removed, 'FaceAlpha', 0.25, 'EdgeColor', style.removed, 'LineWidth', 0.5, ...
          'HandleVisibility', 'off');
end
