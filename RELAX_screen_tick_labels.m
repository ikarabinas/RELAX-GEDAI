% RELAX_screen_tick_labels() - Trace-axis tick labels with each interpolated channel's removal reason.
%
% Usage:
%   >> labels = RELAX_screen_tick_labels(tickLabels, interpolatedReasons);
%
% Inputs:
%   tickLabels          - cell array of channel labels
%   interpolatedReasons - RELAX_screen_channel_reasons output (.labels, .reasons)
%
% Outputs:
%   labels - cell array the size of tickLabels; an interpolated channel reads "E37 (interp., muscle)"

function labels = RELAX_screen_tick_labels(tickLabels, interpolatedReasons)

    ticks = string(tickLabels);
    reasons = interpolatedReasons.reasons(:);
    [isInterpolated, reasonIdx] = ismember(ticks, interpolatedReasons.labels);
    ticks(isInterpolated) = ticks(isInterpolated) + " (interp., " + reasons(reasonIdx(isInterpolated)) + ")";
    labels = reshape(cellstr(ticks), size(tickLabels));
end
