% RELAX_QC_mark_clipped() - Label an axes "clipped" when its limits hide some samples.

function RELAX_QC_mark_clipped(ax, isClipped, style)
    if ~isClipped, return, end
    text(ax, 0.99, 0.02, 'clipped', 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'bottom', 'Color', style.muted, 'FontSize', style.fontSize - 1, 'FontAngle', 'italic');
end
