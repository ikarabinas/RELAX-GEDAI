% RELAX_QC_fig_examples() - Figure 1: before/after traces at the largest blink, the
%                           strongest muscle burst, the largest pop and the largest pop
%                           outside later-excised periods, one figure per example file.

function caption = RELAX_QC_fig_examples(cache, outDir)

    style = RELAX_QC_style();
    titles = struct('blink', "Largest blink", 'muscle', "Strongest muscle burst", 'pop', "Largest pop", 'popKept', "Largest pop kept");
    beforeNames = struct('blink', "Before wICA", 'muscle', "Before wICA", 'pop', "Raw input", 'popKept', "Raw input");
    afterNames  = struct('blink', "Cleaned", 'muscle', "Cleaned", 'pop', "Cleaned", 'popKept', "Cleaned");

    caption = string(sprintf(['## Figure 1 — Example traces\n\n' ...
        'One figure per example file (%d): the median file and the worst files by blink amplitude, muscle power and ' ...
        'pop count, unless hand-picked. Each 6 s window shows the blink electrodes and the two channels with the ' ...
        'most 20–80 Hz power; the two pop columns add the repaired channel at the top. Blink and muscle columns: grey is the sample-matched data before wICA, blue is the ' ...
        'cleaned file, and the strip below is what wICA and GEDAI removed. ' ...
        'Pop columns: the left shows the file''s own pop step with the largest amplitude, the right the largest one whose transition does not overlap any period the pipeline later deletes as extreme ' ...
        '(the next-largest step when that is the same pop as on the left; a note when none qualifies). ' ...
        'Grey is the raw input, unfiltered, at the input sample rate and original reference, with its median over the first second of the window subtracted so the stacked traces are readable; blue is the final cleaned file ' ...
        '(filtered, resampled, average-referenced, and cleaned by every later stage), laid back onto the pre-excision timeline. ' ...
        'The strip below is grey minus blue on the repaired channel only: the pop step and slow drift, but also everything else the pipeline removed or changed (filtering, re-referencing, MWF, wICA, GEDAI), so it is not the pop repair alone. ' ...
        'The repaired interval is shaded in orange, widened to at least 50 ms so it can be seen. ' ...
        'The pop columns are on the pre-excision timeline, so periods the pipeline later deletes as extreme are shaded light grey and labelled "later excised"; the cleaned trace has a gap there. ' ...
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

function is = isPopKind(w), is = any(w.kind == ["pop", "popKept"]); end

function label = removedLabel(w, isPop)
    label = 'Removed (\muV)';
    if isPop, label = [w.channels{1} ' removed (\muV)']; end
end

function drawWindow(traceAx, removedAx, w, beforeName, afterName, style)
    RELAX_QC_draw_traces(traceAx, w, beforeName, afterName, style);

    hold(removedAx, 'on');
    isPop = isPopKind(w);
    [removedTraces, removedRange, removedClipped] = RELAX_QC_removed_strip(w.removed, isPop);
    RELAX_QC_shade_window(removedAx, w, removedRange, false, style);
    plot(removedAx, w.t, removedTraces', 'Color', [style.removed 0.6], 'LineWidth', style.thin);
    for b = w.boundaries, xline(removedAx, b, '--', 'Color', style.ink, 'LineWidth', 0.8); end
    xlim(removedAx, [w.t(1) w.t(end)]); ylim(removedAx, removedRange);
    RELAX_QC_mark_clipped(removedAx, removedClipped, style);
    xlabel(removedAx, 'Time (s)'); ylabel(removedAx, removedLabel(w, isPop)); removedAx.FontSize = style.fontSize;
end
