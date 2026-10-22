% RELAX_QC_shade_window() - Shade a window's later-excised periods and its highlighted interval.
%
% Usage:
%   >> RELAX_QC_shade_window(ax, w, yRange, labelExcised, style);
%
% Excised periods (w.excised) are light grey, labelled "later excised" when labelExcised; the
% highlighted interval (w.shade) is orange, widened to at least 50 ms so it can be seen.

function RELAX_QC_shade_window(ax, w, yRange, labelExcised, style)
    shadeExcised(ax, w.excised, yRange, labelExcised, style);
    shade(ax, w.shade, yRange, style);
end

function seconds = minimumShadeSeconds(), seconds = 0.05; end

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
