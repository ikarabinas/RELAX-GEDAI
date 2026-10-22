% RELAX_QC_draw_traces() - Stacked before/after traces of one window on one axes.
%
% Usage:
%   >> handles = RELAX_QC_draw_traces(ax, w, beforeName, afterName, style);
%
% Draws the shading and excised patches, the traces (one row per channel, first channel at the
% top), boundaries, channel tick labels, legend and "clipped" mark. Pop windows keep each
% channel's own level; other windows are centred on the before median. When w.before is empty
% the offsets and spacing come from w.after and the legend lists only afterName.
%
% Outputs:
%   handles.before, handles.after - line handles, one per channel ([] for a missing before)

function handles = RELAX_QC_draw_traces(ax, w, beforeName, afterName, style)
    hasBefore = ~isempty(w.before);
    reference = w.after;
    if hasBefore, reference = w.before; end
    isPop = any(w.kind == ["pop", "popKept"]);
    centre = median(reference, 2) * ~isPop;
    spacing = max(8 * median(mad(reference, 1, 2)), eps);
    offsets = (size(reference, 1) - 1:-1:0)' * spacing;
    yRange = [-1.5 * spacing, max(offsets) + 3 * spacing];

    hold(ax, 'on');
    RELAX_QC_shade_window(ax, w, yRange, true, style);
    handles.before = [];
    if hasBefore
        handles.before = plot(ax, w.t, w.before - centre + offsets, 'Color', style.before, 'LineWidth', style.thin);
    end
    handles.after = plot(ax, w.t, w.after - centre + offsets, 'Color', style.cleaned, 'LineWidth', style.thin);
    for b = w.boundaries, xline(ax, b, '--', 'Color', style.ink, 'LineWidth', 0.8, 'HandleVisibility', 'off'); end
    [ticks, order] = sort(offsets);
    yticks(ax, ticks); yticklabels(ax, w.channels(order));
    ylim(ax, yRange); xlim(ax, [w.t(1) w.t(end)]);
    ax.XTickLabel = []; ax.FontSize = style.fontSize; ax.TickLabelInterpreter = 'none';
    shown = handles.after(1); names = string(afterName);
    if hasBefore, shown = [handles.before(1); shown]; names = [string(beforeName), names]; end
    legend(ax, shown, cellstr(names), 'Location', 'northoutside', ...
           'Orientation', 'horizontal', 'Box', 'off', 'FontSize', style.fontSize - 1);
    traceStacked = w.after - centre + offsets;
    if hasBefore, traceStacked = [w.before - centre + offsets; traceStacked]; end
    RELAX_QC_mark_clipped(ax, any(traceStacked(:) < yRange(1) | traceStacked(:) > yRange(2)), style);
end
