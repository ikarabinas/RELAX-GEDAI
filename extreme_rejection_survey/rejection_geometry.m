function g = rejection_geometry(mask, segments, windows, srate)
% Shape of one rule's deletion. mask is the RELAX sample mask: 0 = evaluated and kept,
% 1 = rejected, NaN = never evaluated. Never-evaluated samples are not counted as
% rejected stretches. A stretch "splits" a piece of data when kept data lies on both
% sides of it within the same inter-boundary segment, i.e. deleting it adds a new
% discontinuity rather than trimming one that already exists.

    rejected = mask == 1;
    kept = mask == 0;
    lengths = [];
    g.nSplitting = 0;
    for s = 1:size(segments, 1)
        span = segments(s,1):segments(s,2);
        edges = diff([false, rejected(span), false]);
        starts = find(edges == 1);
        ends = find(edges == -1) - 1;
        segKept = kept(span);
        for i = 1:numel(starts)
            keptBefore = starts(i) > 1 && segKept(starts(i) - 1);
            keptAfter = ends(i) < numel(span) && segKept(ends(i) + 1);
            g.nSplitting = g.nSplitting + (keptBefore && keptAfter);
        end
        lengths = [lengths, ends - starts + 1]; %#ok<AGROW>
    end
    g.nStretches = numel(lengths);
    g.medianStretchSeconds = median(lengths) / srate;
    g.nPostHit = count_windows_touching(rejected, windows.post);
    g.nPreHit = count_windows_touching(rejected, windows.pre);
end
