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
% From top to bottom: the file, its outcome and the reason (RELAX_screen_reason); a table of
% each metric's value, dataset median, robust z, cutoff and flag, with a strip plot of every
% cleaned file's z beside it; the three loudest remaining muscle bursts and the three loudest
% remaining blinks as full-width trace windows; and the blink-locked topographies (before wICA,
% cleaned, removed, one colour scale, interpolated channels circled) with the blink-electrode
% time course. Windows without a before-wICA trace, maps without a baseline and files
% without blinks say so in place of the missing panel.

function RELAX_screen_draw_sheet(sheet, row, screen, meta, reviewDir)

    style = RELAX_QC_style();
    place = sheetLayout();
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 place.widthCm place.heightCm], 'Color', 'w');
    page = axes(fig, 'Position', [0 0 1 1], 'XLim', [0 1], 'YLim', [0 1], 'Visible', 'off');

    drawHeader(page, row, RELAX_screen_reason(row, meta, sheet.baselineNote), place, style);
    drawMetricsTable(fig, page, row, screen, meta, place, style);
    drawTraceRows(fig, page, "Muscle burst", sheet.muscle, sheet.muscleRelativePower, "%.1f× file median 20–80 Hz power", ...
                  place.muscleLabelFromTop, place.muscleTopsFromTop, place, style);
    drawTraceRows(fig, page, "Blink", sheet.blink, sheet.blinkResidual, "%.0f µV left at the blink peak", ...
                  place.blinkLabelFromTop, place.blinkTopsFromTop, place, style);
    drawTopographyRow(fig, page, sheet, place, style);

    RELAX_QC_save_figure(fig, reviewDir, char(row.File));
end

function place = sheetLayout()
    place.widthCm = 30;
    place.heightCm = 66;
    place.left = 0.07;
    place.right = 0.97;
    place.titleFromTop = 0.5;
    place.reasonFromTop = 1.8;
    place.tableTopFromTop = 5.0;
    place.tableBottomFromTop = 10.5;
    place.traceHeightCm = 3.8;
    place.muscleLabelFromTop = 13.2;
    place.muscleTopsFromTop = 14.8 + (0:2) * 6.4;
    place.blinkLabelFromTop = 34.8;
    place.blinkTopsFromTop = 36.6 + (0:2) * 6.4;
    place.topoLabelFromTop = 57.0;
    place.topoTopFromTop = 57.5;
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
    lineHeightCm = (place.tableBottomFromTop - place.tableTopFromTop) / nLines;
    lineCentre = @(k) place.fromTop(place.tableTopFromTop + (k - 0.5) * lineHeightCm);
    stripHeight = 0.6 * lineHeightCm / place.heightCm;
    fontSize = style.fontSize + 2;

    headers = struct('metric', "Metric", 'value', "This file", 'median', "Dataset median", 'z', "Robust z", ...
                     'cutoff', "Cutoff", 'flag', "Flag");
    for field = string(fieldnames(columns))'
        text(page, columns.(field), lineCentre(1), headers.(field), 'FontWeight', 'bold', 'FontSize', fontSize, 'Color', style.muted);
    end

    cleanedRows = screen.PipelineOutcome == "cleaned";
    zLimits = commonZLimits(screen, meta, cleanedRows);
    for m = 1:numel(meta.metricNames)
        name = meta.metricNames(m);
        flagged = row.(name + "_flag");
        colour = style.ink; weight = 'normal';
        if flagged, colour = style.channelScreens(5, :); weight = 'bold'; end
        cells = struct('metric', name, 'value', numberText(row.(name)), 'median', numberText(row.(name + "_median")), ...
                       'z', zText(row.(name + "_z"), name, meta), 'cutoff', sprintf("%.1f", meta.cutoff), ...
                       'flag', flagText(flagged));
        for field = string(fieldnames(columns))'
            text(page, columns.(field), lineCentre(m + 1), cells.(field), 'FontSize', fontSize, 'Color', colour, ...
                 'FontWeight', weight, 'Interpreter', 'none');
        end
        stripAxes = axes(fig, 'Position', [0.63, lineCentre(m + 1) - stripHeight / 2, place.right - 0.63, stripHeight]);
        drawStrip(stripAxes, screen.(name + "_z")(cleanedRows), row.(name + "_z"), name, meta, zLimits, m == numel(meta.metricNames), style);
    end
end

function limits = commonZLimits(screen, meta, cleanedRows)
    allZ = cellfun(@(name) screen.(name + "_z")(cleanedRows), cellstr(meta.metricNames), 'UniformOutput', false);
    allZ = vertcat(allZ{:});
    limits = [min([allZ; 0], [], 'omitnan') - 0.5, max([allZ; meta.cutoff], [], 'omitnan') + 1];
end

function drawStrip(ax, cleanedZ, thisZ, name, meta, zLimits, isLast, style)
    hold(ax, 'on');
    xlim(ax, zLimits); ylim(ax, [0 1]);
    ax.YTick = []; ax.FontSize = style.fontSize; ax.Box = 'off';
    if ~isLast, ax.XTickLabel = []; end
    if isLast, xlabel(ax, 'Robust z of every cleaned file'); end
    xline(ax, meta.cutoff, '--', 'Color', style.ink, 'LineWidth', 1);
    if ismember(name, meta.zeroScaleMetrics)
        text(ax, mean(zLimits), 0.5, 'no spread across files: no z-scores', 'HorizontalAlignment', 'center', ...
             'Color', style.muted, 'FontSize', style.fontSize);
        return
    end
    present = cleanedZ(~isnan(cleanedZ));
    scatter(ax, present, deterministicJitter(numel(present)), 14, style.before, 'filled', 'MarkerFaceAlpha', 0.6);
    if ~isnan(thisZ)
        scatter(ax, thisZ, 0.5, 60, style.channelScreens(5, :), 'filled', 'MarkerEdgeColor', style.ink, 'LineWidth', 1);
    end
end

function y = deterministicJitter(n)
    goldenRatio = (sqrt(5) - 1) / 2;
    y = 0.15 + 0.7 * mod((1:n)' * goldenRatio, 1);
end

function s = numberText(x)
    if isnan(x), s = "missing"; else, s = sprintf("%.3g", x); end
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

function drawTraceRows(fig, page, label, windows, strengths, strengthFormat, labelFromTop, topsFromTop, place, style)
    text(page, place.left, place.fromTop(labelFromTop), label + "s, loudest first", 'FontWeight', 'bold', ...
         'FontSize', style.fontSize + 3, 'Color', style.ink, 'VerticalAlignment', 'bottom');
    for k = 1:numel(windows)
        w = windows(k);
        ax = axes(fig, 'Position', [place.left, place.bottomFromTop(topsFromTop(k), place.traceHeightCm), place.right - place.left, place.traceHeightCm / place.heightCm]);
        heading = sprintf("%s %d", label, k);
        if ~w.available
            RELAX_QC_draw_note(ax, heading + ": " + w.note);
            continue
        end
        RELAX_QC_draw_traces(ax, w, "Before wICA", "Cleaned", style);
        ax.XTickLabelMode = 'auto';
        xlabel(ax, heading + sprintf(", starts at %.1f s: ", w.t(1)) + sprintf(strengthFormat, strengths(k)), ...
               'FontSize', style.fontSize + 2, 'FontWeight', 'bold');
    end
end

function drawTopographyRow(fig, page, sheet, place, style)
    height = place.topoHeightCm / place.heightCm;
    bottom = place.bottomFromTop(place.topoTopFromTop, place.topoHeightCm);
    text(page, place.left, place.fromTop(place.topoLabelFromTop), blinkAverageHeading(sheet), 'FontWeight', 'bold', ...
         'FontSize', style.fontSize + 3, 'Color', style.ink, 'VerticalAlignment', 'bottom', 'Interpreter', 'none');
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
    electrodes = electrodeOptions(find(ismember({sheet.chanlocs.labels}, sheet.interpolatedLabels)));
    for j = 1:3
        ax = mapAxes(j);
        if j ~= 2 && ~hasBaseline
            RELAX_QC_draw_note(ax, titles(j) + ":" + newline + "no before-wICA baseline" + newline + "(reason in the header)");
            continue
        end
        axes(ax); %#ok<LAXES> topoplot draws into the current axes
        topoplot(maps{j}, sheet.chanlocs, 'maplimits', [-limit limit], 'colormap', style.diverging, ...
                 electrodes{:});
        title(ax, titles(j), 'FontSize', style.fontSize + 2, 'FontWeight', 'normal');
    end
    scale = colorbar(mapAxes(3), 'Position', [mapLefts(3) + mapWidth + 0.005, bottom + 0.3 * height, 0.01, 0.4 * height]);
    scale.Title.String = '\muV';
    colormap(mapAxes(3), style.diverging); clim(mapAxes(3), [-limit limit]);
    drawBlinkTimeCourse(timeAxes, avg, hasBaseline, style);
end

function options = electrodeOptions(interpolatedIdx)
    options = {'electrodes', 'on'};
    if ~isempty(interpolatedIdx)
        options = [options, {'emarker2', {interpolatedIdx, 'o', 'k', 4, 1}}];
    end
end

function heading = blinkAverageHeading(sheet)
    heading = "Blink-locked average over all " + sheet.blinkN + " blinks (baseline −500 to −300 ms, ±50 ms at the peak)";
    if ~isempty(sheet.interpolatedLabels)
        heading = heading + newline + "Circled electrodes are interpolated: " + strjoin(string(sheet.interpolatedLabels), ", ");
    end
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
