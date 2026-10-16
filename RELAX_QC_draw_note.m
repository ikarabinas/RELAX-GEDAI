% RELAX_QC_draw_note() - Put a centred note where a panel has nothing to plot.

function RELAX_QC_draw_note(ax, note)
    style = RELAX_QC_style();
    axis(ax, 'off');
    text(ax, 0.5, 0.5, note, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
         'Color', style.muted, 'FontSize', style.fontSize, 'Interpreter', 'none');
end
