% RELAX_QC_fig_examples() - Figure 1: before/after traces at the largest blink, the
%                           strongest muscle burst, the largest pop and the largest pop
%                           outside later-excised periods, one figure per example file.

function caption = RELAX_QC_fig_examples(cache, outDir)

    style = RELAX_QC_style();
    titles = struct('blink', "Largest blink", 'muscle', "Strongest muscle burst", 'pop', "Largest pop", 'popKept', "Largest pop kept");
    beforeNames = struct('blink', "Before wICA", 'muscle', "Before wICA", 'pop', "No pop repair", 'popKept', "No pop repair");
    afterNames  = struct('blink', "Cleaned", 'muscle', "Cleaned", 'pop', "Pop-repaired", 'popKept', "Pop-repaired");

    caption = string(sprintf(['## Figure 1 — Example traces\n\n' ...
        'One figure per example file (%d): the median file and the worst files by blink amplitude, muscle power and ' ...
        'pop count, unless hand-picked. Each 6 s window shows the blink electrodes and the two channels with the ' ...
        'most 20–80 Hz power; the two pop columns add the repaired channel at the top. Blink and muscle columns: grey is the sample-matched data before wICA, blue is the ' ...
        'cleaned file, and the strip below is what wICA removed. ' ...
        'Pop columns: the left shows the file''s own pop step with the largest amplitude, the right the largest one whose transition does not overlap any period the pipeline later deletes as extreme ' ...
        '(the next-largest step when that is the same pop as on the left; a note when none qualifies). ' ...
        'Both are drawn on unfiltered data at the input sample rate, one offset per channel (its median over the first second of the window) subtracted so the stacked traces are readable. ' ...
        'Grey is the raw input; blue is the same input with the pipeline''s pop repair replayed on it (RELAX_repair_steps with the saved steps and the participant''s saved options). ' ...
        'The replay is a display reconstruction, not saved pipeline data, and uses the input restricted to the channels of Extremes_Rejected as the montage. ' ...
        '%s' ...
        'The strip below is grey minus blue on the repaired channel only: the step itself and the slow drift that the detrend removes across the whole channel, not just at the pop. ' ...
        'The repaired interval is shaded in orange, widened to at least 50 ms so it can be seen. ' ...
        'The pop columns are on the pre-excision timeline, so periods the pipeline later deletes as extreme are shaded light grey and labelled "later excised"; wild traces there are not repair damage. ' ...
        'Traces and the pop removed strip are clipped to fixed limits, the pop removed strip at 1.2 times its 99th-percentile absolute value (marked "clipped"). Dashed black lines are boundaries. Good cleaning: the ' ...
        'blink and burst are gone from blue while the background rhythm is unchanged, and the removed strip ' ...
        'carries only the artifact.\n'], numel(cache.examples), referencePopSentence(cache)));

    if isempty(cache.examples)
        fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 18 6], 'Color', 'w');
        RELAX_QC_draw_note(axes(fig), "No example files: no aligned file had blink, muscle or pop scores");
        RELAX_QC_save_figure(fig, outDir, 'fig1_examples_none');
    end

    for k = 1:numel(cache.examples)
        example = cache.examples(k);
        fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 36 15], 'Color', 'w');
        layout = tiledlayout(fig, 4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(layout, example.file + " (" + example.role + ")", 'Interpreter', 'none', 'FontSize', style.fontSize + 1);
        for j = 1:numel(example.windows)
            w = example.windows(j);
            kind = char(w.kind);
            traceAx = nexttile(layout, j, [3 1]);
            removedAx = nexttile(layout, 12 + j);
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

function sentence = referencePopSentence(cache)
    replayed = arrayfun(@(e) any(contains([e.windows.note], 'reference pops subtracted')), cache.examples);
    sentence = "No reference pops needed replaying in these files (reference-pop repair off, or none detected). ";
    if any(replayed)
        sentence = "Reference pops were subtracted first in the replay for " + sum(replayed) + " of " + numel(replayed) + " of these files. ";
    end
end

function seconds = minimumShadeSeconds(), seconds = 0.05; end

function is = isPopKind(w), is = any(w.kind == ["pop", "popKept"]); end

function label = removedLabel(w, isPop)
    label = 'Removed (\muV)';
    if isPop, label = [w.channels{1} ' removed (\muV)']; end
end

function drawWindow(traceAx, removedAx, w, beforeName, afterName, style)
    before = w.before; after = w.after;
    centre = median(before, 2) * ~isPopKind(w);
    spacing = max(8 * median(mad(before, 1, 2)), eps);
    offsets = (size(before, 1) - 1:-1:0)' * spacing;
    yRange = [-1.5 * spacing, max(offsets) + 3 * spacing];

    hold(traceAx, 'on');
    shadeExcised(traceAx, w.excised, yRange, true, style);
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
    isPop = isPopKind(w);
    [removedTraces, removedRange, removedClipped] = RELAX_QC_removed_strip(w.removed, isPop);
    shadeExcised(removedAx, w.excised, removedRange, false, style);
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

function shadeExcised(ax, intervals, yRange, labelled, style)
    for k = 1:size(intervals, 1)
        patch(ax, intervals(k, [1 2 2 1]), yRange([1 1 2 2]), [0.85 0.85 0.85], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        if labelled
            text(ax, mean(intervals(k, :)), yRange(2) - 0.02 * diff(yRange), 'later excised', ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', style.muted, 'FontSize', style.fontSize - 1);
        end
    end
end
