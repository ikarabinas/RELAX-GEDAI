% RELAX_QC_fill_window() - Fill a RELAX_QC_empty_window with the samples around one centre.
%
% Usage:
%   >> w = RELAX_QC_fill_window(w, centre, before, after, labels, windowSeconds);
%
% Inputs:
%   centre        - window centre in seconds on the timeline of after
%   before, after - EEG structs on the same timeline; before = [] for a window with no
%                   baseline, which is filled from after alone (w.before and w.removed stay [])
%   labels        - channels to show, kept in this order when present in both structs
%   windowSeconds - window length; the window slides to stay inside the file

function w = RELAX_QC_fill_window(w, centre, before, after, labels, windowSeconds)
    hasBefore = ~isempty(before);
    present = labels(ismember(labels, {after.chanlocs.labels}));
    if hasBefore, present = present(ismember(present, {before.chanlocs.labels})); end
    [~, afterRows]  = ismember(present, {after.chanlocs.labels});
    half = round(windowSeconds / 2 * after.srate);
    middle = round(centre * after.srate) + 1;
    first = max(1, min(middle - half, after.pnts - 2 * half));
    span = first : min(after.pnts, first + 2 * half);

    w.channels = present;
    w.t = (span - 1) / after.srate;
    w.after = double(after.data(afterRows, span));
    if hasBefore
        [~, beforeRows] = ismember(present, {before.chanlocs.labels});
        w.before = double(before.data(beforeRows, span));
        w.removed = w.before - w.after;
    end
    segments = RELAX_get_segments(after);
    starts = segments(2:end, 1);
    w.boundaries = (starts(starts > span(1) & starts <= span(end))' - 1.5) / after.srate;
    w.available = true;
end
