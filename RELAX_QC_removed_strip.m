% RELAX_QC_removed_strip() - Which removed traces Figure 1 plots in a window's strip and the y-limits.
%
% Blink and muscle strips show every channel over the full data range. The pop strip shows only
% the repaired channel (the first), because the other channels differ between the reconstruction
% and the pipeline only by processing choices; its limits are robust so the pop's own residual
% is not flattened by the step it replaced, and clipped says whether any sample falls outside.

function [traces, yRange, clipped] = RELAX_QC_removed_strip(removed, isPop)
    minimumLimit = 1;
    if isPop
        traces = removed(1, :);
        limit = max(1.2 * prctile(abs(traces), 99), minimumLimit);
        yRange = [-limit limit];
        clipped = any(abs(traces) > limit);
    else
        traces = removed;
        yRange = [min(traces, [], 'all'), max(traces, [], 'all')];
        padding = 0.05 * max(diff(yRange), minimumLimit);
        yRange = yRange + [-padding padding];
        clipped = false;
    end
end
