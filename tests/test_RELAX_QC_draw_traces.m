function test_RELAX_QC_draw_traces()
% Traces are drawn with or without a before-wICA baseline; without one only the after legend entry exists.

    style = RELAX_QC_style();
    w = RELAX_QC_empty_window("muscle");
    w.available = true; w.channels = {'E1', 'E2', 'E3'}; w.t = (0:299) / 100;
    w.before = 5 + randn(3, 300); w.after = 5 + 0.5 * randn(3, 300); w.removed = w.before - w.after;
    w.boundaries = 1.5;
    fig = figure('Visible', 'off'); cleanup = onCleanup(@() close(fig));

    ax = axes(fig);
    handles = RELAX_QC_draw_traces(ax, w, "Before wICA", "Cleaned", style);
    assert(numel(handles.before) == 3 && numel(handles.after) == 3);
    assert(isequal(ax.YTickLabel(:)', {'E3', 'E2', 'E1'}) && ax.XLim(2) == w.t(end));

    ax = axes(fig);
    w.before = []; w.removed = [];
    handles = RELAX_QC_draw_traces(ax, w, "Before wICA", "Cleaned", style);
    assert(isempty(handles.before) && numel(handles.after) == 3);
    legendEntries = findobj(fig, 'Type', 'Legend');
    assert(isequal(legendEntries(1).String, {'Cleaned'}), 'only the after entry is listed without a baseline');

    fprintf('test_RELAX_QC_draw_traces: PASS\n');
end
