% RELAX_screen_draw_sheet() - One file's review sheet: why it passed or failed the residual-noise
%                             screen, and what its loudest remaining muscle bursts and blinks
%                             look like.
%
% Usage:
%   >> RELAX_screen_draw_sheet(sheet, row, screen, meta, reviewDir);
%
% Inputs:
%   sheet     - RELAX_screen_sheet_data output for this file
%   row       - this file's row of screen
%   screen    - the whole table from RELAX_screen_score_table (for the strip plots)
%   meta      - its second output
%   reviewDir - folder that receives <row.File>.png and .pdf
%
% From top to bottom: the file, its outcome and the reason (RELAX_screen_reason); a table with a
% row for each of meta.metricNames giving the value, dataset median, robust z, cutoff and flag
% (the hole depth row gives its limit as the cutoff and no median or z, the interpolated-channel
% row its channel limit), with a strip plot of every cleaned file beside the muscle, blink and hole
% depth rows (robust z on a shared axis, hole depth on its own); the three loudest remaining
% muscle bursts and the three loudest remaining blinks as full-width trace windows; and the
% blink-locked topographies (before wICA, cleaned, removed, one colour scale, interpolated channels
% circled, hole electrodes filled in amber) with the blink-electrode time course; the interpolated and
% hole electrodes are listed by label above the maps, wrapped over as many lines as they need, whether
% or not the file has blinks.
% Windows without a before-wICA trace, maps without a baseline and files without blinks say so in
% place of the missing panel.

function RELAX_screen_draw_sheet(sheet, row, screen, meta, reviewDir)

    style = RELAX_QC_style();
    electrodeLines = electrodeListLines(sheet, row, meta);
    place = sheetLayout(numel(electrodeLines));
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 place.widthCm place.heightCm], 'Color', 'w');
    page = axes(fig, 'Position', [0 0 1 1], 'XLim', [0 1], 'YLim', [0 1], 'Visible', 'off');

    drawHeader(page, row, RELAX_screen_reason(row, meta, sheet.baselineNote), place, style);
    drawMetricsTable(fig, page, row, screen, meta, place, style);
    drawTraceRows(fig, page, "Muscle burst", sheet.muscle, sheet.muscleRelativePower, "%.1f× file median 20–80 Hz power", ...
                  sheet.interpolatedLabels, place.muscleLabelFromTop, place.muscleTopsFromTop, place, style);
    drawTraceRows(fig, page, "Blink", sheet.blink, sheet.blinkResidual, "%.0f µV left at the blink peak", ...
                  sheet.interpolatedLabels, place.blinkLabelFromTop, place.blinkTopsFromTop, place, style);
    drawTopographyRow(fig, page, sheet, row, electrodeLines, place, style);

    RELAX_QC_save_figure(fig, reviewDir, char(row.File));
end

function place = sheetLayout(nElectrodeLines)
    place.electrodeLineCm = 0.55;
    electrodeListCm = nElectrodeLines * place.electrodeLineCm;
    place.widthCm = 30;
    place.heightCm = 74.4 + electrodeListCm;
    place.left = 0.07;
    place.right = 0.97;
    place.titleFromTop = 0.5;
    place.reasonFromTop = 1.8;
    place.tableTopFromTop = 4.8;
    place.tableBottomFromTop = 12.0;
    place.stripAxisGapCm = 1.2;
    place.traceHeightCm = 5.8;
    place.muscleLabelFromTop = 13.6;
    place.muscleTopsFromTop = 14.1 + (0:2) * 7.7;
    place.blinkLabelFromTop = 38.5;
    place.blinkTopsFromTop = 39.0 + (0:2) * 7.7;
    place.topoLabelFromTop = 63.7;
    place.topoTopFromTop = place.topoLabelFromTop + electrodeListCm + 0.6;
    place.topoHeightCm = 7;
    place.fromTop = @(cm) 1 - cm / place.heightCm;
    place.bottomFromTop = @(topCm, heightCm) 1 - (topCm + heightCm) / place.heightCm;
end

function drawHeader(page, row, reason, place, style)
    text(page, place.left, place.fromTop(place.titleFromTop), row.File + " — " + upper(row.ScreenOutcome), 'FontSize', style.fontSize + 6, ...
         'FontWeight', 'bold', 'Interpreter', 'none', 'VerticalAlignment', 'top', 'Color', style.ink);
    reasonLines = strrep(strrep(reason, "; ", newline), " | ", newline);
    text(page, place.left, place.fromTop(place.reasonFromTop), reasonLines, 'FontSize', style.fontSize + 2, 'Interpreter', 'none', ...
         'VerticalAlignment', 'top', 'Color', style.ink);
end

function drawMetricsTable(fig, page, row, screen, meta, place, style)
    columns = struct('metric', place.left, 'value', 0.21, 'median', 0.285, 'z', 0.41, 'cutoff', 0.50, 'flag', 0.56);
    nLines = numel(meta.metricNames) + 1;
    lineHeightCm = (place.tableBottomFromTop - place.tableTopFromTop - 2 * place.stripAxisGapCm) / nLines;
    firstHoleLine = numel(meta.zMetricNames) + 2;
    lineCentre = @(k) place.fromTop(place.tableTopFromTop + (k - 0.5) * lineHeightCm + place.stripAxisGapCm * ((k >= firstHoleLine) + (k > firstHoleLine)));
    stripHeight = 0.6 * lineHeightCm / place.heightCm;
    stripLeft = 0.63;
    stripRight = 0.89;
    fontSize = style.fontSize + 2;

    headers = struct('metric', "Metric", 'value', "This file", 'median', "Dataset median", 'z', "Robust z", ...
                     'cutoff', "Cutoff", 'flag', "Flag");
    for field = string(fieldnames(columns))'
        text(page, columns.(field), lineCentre(1), headers.(field), 'FontWeight', 'bold', 'FontSize', fontSize, 'Color', style.muted);
    end

    cleanedRows = screen.PipelineOutcome == "cleaned";
    zLimits = commonZLimits(screen, meta, cleanedRows);
    holeLimits = [0, max(2.5 * meta.holeDepthLimit, meta.holeDepthLimit + 1)];
    for m = 1:numel(meta.metricNames)
        name = meta.metricNames(m);
        flagged = row.(name + "_flag");
        colour = style.ink; weight = 'normal';
        if flagged, colour = style.channelScreens(5, :); weight = 'bold'; end
        cells = metricCells(row, name, meta);
        cells.metric = name;
        cells.flag = flagText(flagged);
        for field = string(fieldnames(columns))'
            text(page, columns.(field), lineCentre(m + 1), cells.(field), 'FontSize', fontSize, 'Color', colour, ...
                 'FontWeight', weight, 'Interpreter', 'none');
        end
        stripAxes = axes(fig, 'Position', [stripLeft, lineCentre(m + 1) - stripHeight / 2, stripRight - stripLeft, stripHeight]);
        if ismember(name, meta.zMetricNames)
            strip = struct('xLimits', zLimits, 'reference', meta.cutoff, 'noSpread', ismember(name, meta.zeroScaleMetrics), 'label', "");
            if name == meta.zMetricNames(end), strip.label = "Robust z of cleaned files (> = beyond the edge)"; end
            drawStrip(page, stripAxes, screen.(name + "_z")(cleanedRows), row.(name + "_z"), strip, style);
        elseif name == "HoleDepth"
            strip = struct('xLimits', holeLimits, 'reference', row.HoleDepth_limit, 'noSpread', false, ...
                           'label', "Hole depth (× electrode spacing; > = beyond the edge)");
            drawStrip(page, stripAxes, screen.HoleDepth(cleanedRows), row.HoleDepth, strip, style);
        else
            delete(stripAxes);
        end
    end
end

function cells = metricCells(row, name, meta)
    switch name
        case "HoleDepth"
            cells = struct('value', valueText(row.HoleDepth, "%.2f"), 'median', "", 'z', "", ...
                           'cutoff', valueText(row.HoleDepth_limit, "%.1f"));
        case "NInterpolated"
            cells = struct('value', valueText(row.NInterpolated, "%d"), 'median', "", 'z', "", ...
                           'cutoff', valueText(row.NInterpolated_limit, "%d"));
        otherwise
            cells = struct('value', valueText(row.(name), "%.3g"), 'median', valueText(row.(name + "_median"), "%.3g"), ...
                           'z', zText(row.(name + "_z"), name, meta), 'cutoff', sprintf("%.1f", meta.cutoff));
    end
end

function limits = commonZLimits(screen, meta, cleanedRows)
% One range for every z strip, capped so a single wild file cannot flatten the rest.
    allZ = cellfun(@(name) screen.(name + "_z")(cleanedRows), cellstr(meta.zMetricNames), 'UniformOutput', false);
    allZ = vertcat(allZ{:});
    limits = [min([allZ; -3], [], 'omitnan') - 0.5, max(3 * meta.cutoff, meta.cutoff + 1)];
end

function drawStrip(page, ax, cleaned, this, strip, style)
% Every cleaned file as a dot, this file as the large dot, the dashed line at strip.reference. A value
% beyond strip.xLimits(2) is drawn at the edge as '>' and counted beside the strip.
    hold(ax, 'on');
    xlim(ax, strip.xLimits); ylim(ax, [0 1]);
    ax.YTick = []; ax.FontSize = style.fontSize; ax.Box = 'off';
    if strlength(strip.label) == 0, ax.XTickLabel = []; else, xlabel(ax, strip.label); end
    xline(ax, strip.reference, '--', 'Color', style.ink, 'LineWidth', 1);
    if strip.noSpread
        text(ax, mean(strip.xLimits), 0.5, 'no spread across files: no z-scores', 'HorizontalAlignment', 'center', ...
             'Color', style.muted, 'FontSize', style.fontSize);
        return
    end
    present = cleaned(~isnan(cleaned));
    offScale = present > strip.xLimits(2);
    jitter = deterministicJitter(numel(present));
    scatter(ax, present(~offScale), jitter(~offScale), 14, style.before, 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(ax, repmat(strip.xLimits(2), 1, nnz(offScale)), jitter(offScale), 20, style.before, '>', 'filled');
    if ~isnan(this)
        marker = 'o'; if this > strip.xLimits(2), marker = '>'; end
        scatter(ax, min(this, strip.xLimits(2)), 0.5, 70, style.channelScreens(5, :), marker, 'filled', 'MarkerEdgeColor', style.ink, 'LineWidth', 1);
    end
    if any(offScale)
        position = ax.Position;
        text(page, position(1) + position(3) + 0.01, position(2) + position(4) / 2, sprintf('%d off-scale', nnz(offScale)), ...
             'FontSize', style.fontSize + 1, 'Color', style.muted, 'VerticalAlignment', 'middle');
    end
end

function y = deterministicJitter(n)
    goldenRatio = (sqrt(5) - 1) / 2;
    y = 0.15 + 0.7 * mod((1:n)' * goldenRatio, 1);
end

function s = valueText(x, format)
    if isnan(x), s = "missing"; else, s = sprintf(format, x); end
end

function s = zText(z, name, meta)
    if ismember(name, meta.zeroScaleMetrics), s = "no spread";
    elseif isnan(z), s = "missing";
    else, s = sprintf("%.1f", z);
    end
end

function s = flagText(flagged)
    if flagged, s = "FLAGGED"; else, s = "no"; end
end

function drawTraceRows(fig, page, label, windows, strengths, strengthFormat, interpolatedLabels, labelFromTop, topsFromTop, place, style)
    text(page, place.left, place.fromTop(labelFromTop), label + "s, loudest first", 'FontWeight', 'bold', ...
         'FontSize', style.fontSize + 3, 'Color', style.ink, 'VerticalAlignment', 'bottom');
    drawTraceKey(page, any(arrayfun(@(w) ~isempty(w.before), windows)), place.fromTop(labelFromTop) + 0.003, style);
    for k = 1:numel(windows)
        w = windows(k);
        ax = axes(fig, 'Position', [place.left, place.bottomFromTop(topsFromTop(k), place.traceHeightCm), place.right - place.left, place.traceHeightCm / place.heightCm]);
        heading = sprintf("%s %d", label, k);
        if ~w.available
            RELAX_QC_draw_note(ax, heading + ": " + w.note);
            continue
        end
        handles = RELAX_QC_draw_traces(ax, w, "Before wICA", "Cleaned", style);
        legend(ax, 'off');
        trimEmptyHeadroom(ax, handles, numel(w.channels), style);
        markInterpolatedTicks(ax, interpolatedLabels);
        drawScaleBar(ax, handles.spacing, style);
        ax.XTickLabelMode = 'auto';
        xlabel(ax, heading + sprintf(", starts at %.1f s: ", w.t(1)) + sprintf(strengthFormat, strengths(k)), ...
               'FontSize', style.fontSize + 2, 'FontWeight', 'bold');
    end
end

function drawTraceKey(page, hasBefore, y, style)
    entries = {"Cleaned", style.cleaned};
    if hasBefore, entries = [{"Before wICA", style.before}; entries]; end
    x = 0.36;
    for k = 1:size(entries, 1)
        line(page, [x x + 0.03], [y y], 'Color', entries{k, 2}, 'LineWidth', 2);
        text(page, x + 0.038, y, entries{k, 1}, 'FontSize', style.fontSize + 2, 'Color', style.ink, 'VerticalAlignment', 'middle');
        x = x + 0.17;
    end
end

function trimEmptyHeadroom(ax, handles, nChannels, style)
% RELAX_QC_draw_traces leaves room for a legend and a clipping mark; a sheet gives that room to the traces.
    spacing = handles.spacing;
    delete(findobj(ax, 'Type', 'text', 'String', 'clipped'));
    ylim(ax, [-1.2 * spacing, (nChannels - 1) * spacing + 2 * spacing]);
    yRange = ylim(ax);
    for label = findobj(ax, 'Type', 'text', 'String', 'later excised')'
        label.Position(2) = yRange(2) - 0.02 * diff(yRange);
    end
    plotted = [handles.before; handles.after];
    clipped = arrayfun(@(h) any(h.YData < yRange(1) | h.YData > yRange(2)), plotted);
    RELAX_QC_mark_clipped(ax, any(clipped), style);
end

function markInterpolatedTicks(ax, interpolatedLabels)
    labels = ax.YTickLabel;
    isInterpolated = ismember(labels, interpolatedLabels);
    labels(isInterpolated) = cellfun(@(label) [label ' (interp.)'], labels(isInterpolated), 'UniformOutput', false);
    ax.YTickLabel = labels;
end

function drawScaleBar(ax, spacing, style)
% A vertical bar at the bottom left, as long as the 1-2-5 value just under one trace spacing. The bottom
% right belongs to RELAX_QC_mark_clipped's note.
    barMicrovolts = oneTwoFiveBelow(spacing);
    xRange = xlim(ax);
    x = xRange(1) + 0.01 * diff(xRange);
    bottom = -1.1 * spacing;
    line(ax, [x x], bottom + [0 barMicrovolts], 'Color', style.ink, 'LineWidth', 2, 'HandleVisibility', 'off');
    text(ax, x + 0.006 * diff(xRange), bottom + barMicrovolts / 2, sprintf('%g µV', barMicrovolts), 'FontSize', style.fontSize, ...
         'Color', style.ink, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end

function value = oneTwoFiveBelow(x)
    magnitude = 10 ^ floor(log10(x));
    steps = [1 2 5 10] * magnitude;
    value = steps(find(steps <= x * (1 + 1e-9), 1, 'last'));
end

function drawTopographyRow(fig, page, sheet, row, electrodeLines, place, style)
    height = place.topoHeightCm / place.heightCm;
    bottom = place.bottomFromTop(place.topoTopFromTop, place.topoHeightCm);
    text(page, place.left, place.fromTop(place.topoLabelFromTop), blinkAverageHeading(sheet), 'FontWeight', 'bold', ...
         'FontSize', style.fontSize + 3, 'Color', style.ink, 'VerticalAlignment', 'bottom', 'Interpreter', 'none');
    if ~isempty(electrodeLines)
        text(page, place.left, place.fromTop(place.topoLabelFromTop + 0.05), strjoin(electrodeLines, newline), ...
             'FontSize', style.fontSize + 2, 'Color', style.ink, 'VerticalAlignment', 'top', 'Interpreter', 'none');
    end
    mapWidth = 0.19;
    mapLefts = place.left + (0:2) * (mapWidth + 0.015);
    mapAxes = arrayfun(@(left) axes(fig, 'Position', [left, bottom, mapWidth, height]), mapLefts);
    timeAxes = axes(fig, 'Position', [0.79, bottom + 0.02, place.right - 0.79, height - 0.05]);

    if isempty(sheet.blinkAverage)
        arrayfun(@(ax) RELAX_QC_draw_note(ax, "No blinks in this file"), [mapAxes timeAxes]);
        return
    end
    avg = sheet.blinkAverage;
    titles = ["Before wICA", "Cleaned", "Removed by wICA + GEDAI"];
    maps = {avg.topoBaseline, avg.topoCleaned, avg.topoRemoved};
    hasBaseline = ~any(isnan(avg.topoBaseline));
    limit = max(abs([maps{:}]), [], 'all', 'omitnan');
    if ~(limit > 0), limit = 1; end
    interpolatedIdx = find(ismember({sheet.chanlocs.labels}, sheet.interpolatedLabels));
    electrodes = electrodeOptions(interpolatedIdx);
    holeLabels = holeElectrodeLabels(row);
    drawn = false(1, 3);
    for j = 1:3
        ax = mapAxes(j);
        if j ~= 2 && ~hasBaseline
            RELAX_QC_draw_note(ax, titles(j) + ":" + newline + "no before-wICA baseline" + newline + "(reason in the header)");
            continue
        end
        axes(ax); %#ok<LAXES> topoplot draws into the current axes
        topoplot(maps{j}, sheet.chanlocs, 'maplimits', [-limit limit], 'colormap', style.diverging, ...
                 electrodes{:});
        set(findobj(ax, 'Type', 'line', 'Marker', 'o'), 'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 0.75);
        fillHoleElectrodes(ax, sheet.chanlocs, interpolatedIdx, holeLabels, style);
        drawn(j) = true;
    end
    scale = colorbar(mapAxes(3), 'Position', [mapLefts(3) + mapWidth + 0.005, bottom + 0.3 * height, 0.01, 0.4 * height]);
    scale.Title.String = '\muV';
    colormap(mapAxes(3), style.diverging); clim(mapAxes(3), [-limit limit]);
    for j = find(drawn)
        title(mapAxes(j), titles(j), 'FontSize', style.fontSize + 2, 'FontWeight', 'normal');
    end
    drawBlinkTimeCourse(timeAxes, avg, hasBaseline, style);
end

function options = electrodeOptions(interpolatedIdx)
    options = {'electrodes', 'on'};
    if ~isempty(interpolatedIdx)
        options = [options, {'emarker2', {interpolatedIdx, 'o', 'k', 4, 1}}];
    end
end

function fillHoleElectrodes(ax, chanlocs, interpolatedIdx, holeLabels, style)
% topoplot draws one set of marker-2 electrodes, in ascending channel order: the open rings.
    rings = findobj(ax, 'Type', 'line', 'Marker', 'o');
    isHole = ismember(string({chanlocs(interpolatedIdx).labels}), holeLabels);
    if ~any(isHole), return, end
    if numel(rings) ~= 1 || numel(rings.XData) ~= numel(interpolatedIdx)
        warning('RELAX_screen:holeMarkers', 'Could not match the interpolated rings to electrodes; hole electrodes not filled');
        return
    end
    line(ax, rings.XData(isHole), rings.YData(isHole), rings.ZData(isHole), 'LineStyle', 'none', 'Marker', 'o', ...
         'MarkerSize', 8, 'LineWidth', 0.75, 'MarkerFaceColor', style.channelScreens(1, :), 'MarkerEdgeColor', 'k');
end

function labels = holeElectrodeLabels(row)
    labels = strings(1, 0);
    if strlength(row.HoleElectrodes) > 0, labels = split(row.HoleElectrodes)'; end
end

function heading = blinkAverageHeading(sheet)
    if sheet.blinkN == 0
        heading = "Blink-locked average: no blinks in this file";
    else
        heading = "Blink-locked average over all " + sheet.blinkN + " blinks (baseline −500 to −300 ms, ±50 ms at the peak)";
    end
end

function lines = electrodeListLines(sheet, row, meta)
    holePrefix = sprintf("Hole electrodes (no real electrode within %.1f× spacing):", meta.holeDepthLimit);
    lines = [wrapLabels("Circled electrodes are interpolated:", string(sheet.interpolatedLabels(:))'); ...
             wrapLabels(holePrefix, holeElectrodeLabels(row))];
end

function lines = wrapLabels(prefix, labels)
    maxCharacters = 110;
    lines = strings(0, 1);
    if isempty(labels), return, end
    labels(1:end-1) = labels(1:end-1) + ",";
    current = prefix;
    for label = labels
        if strlength(current) + 1 + strlength(label) > maxCharacters
            lines(end + 1, 1) = current; %#ok<AGROW>
            current = label;
        else
            current = current + " " + label;
        end
    end
    lines(end + 1, 1) = current;
end

function drawBlinkTimeCourse(ax, avg, hasBaseline, style)
    hold(ax, 'on');
    ms = 1000 * avg.blinkTimes;
    if hasBaseline
        plot(ax, ms, avg.erpBaseline, 'Color', style.before, 'LineWidth', style.bold, 'DisplayName', 'Before wICA');
    end
    plot(ax, ms, avg.erpCleaned, 'Color', style.cleaned, 'LineWidth', style.bold, 'DisplayName', 'Cleaned');
    xlabel(ax, 'Time from blink peak (ms)'); ylabel(ax, 'Blink electrodes (\muV)');
    yRange = ylim(ax); ylim(ax, [yRange(1), yRange(2) + 0.3 * diff(yRange)]);
    legend(ax, 'Location', 'northeast', 'Box', 'off'); ax.YGrid = 'on'; ax.FontSize = style.fontSize + 1;
end
