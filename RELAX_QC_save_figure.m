% RELAX_QC_save_figure() - Save a QC figure as 300 dpi PNG and vector PDF, then close it.

function RELAX_QC_save_figure(fig, outDir, name)
    if ~isfolder(outDir), mkdir(outDir); end
    style = RELAX_QC_style();
    set(findall(fig, '-property', 'FontName'), 'FontName', style.fontName);
    exportgraphics(fig, fullfile(outDir, [name '.png']), 'Resolution', 300);
    exportgraphics(fig, fullfile(outDir, [name '.pdf']), 'ContentType', 'vector');
    close(fig);
end
