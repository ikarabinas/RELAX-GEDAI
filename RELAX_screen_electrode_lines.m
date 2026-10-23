% RELAX_screen_electrode_lines() - The electrode lists above a review sheet's maps, wrapped to fit the page.
%
% Usage:
%   >> lines = RELAX_screen_electrode_lines(interpolatedReasons, holeLabels, holeDepthLimit);
%
% Inputs:
%   interpolatedReasons - RELAX_screen_channel_reasons output; its .groups give one list per reason
%   holeLabels          - [1 x h] string, the electrodes in the hole
%   holeDepthLimit      - the hole depth limit, in multiples of electrode spacing
%
% Outputs:
%   lines - [n x 1] string: "Interpolated for <reason> (<count>):" and its labels for each reason in
%           reason order, then the hole electrodes; empty when there are neither

function lines = RELAX_screen_electrode_lines(interpolatedReasons, holeLabels, holeDepthLimit)

    holePrefix = sprintf("Hole electrodes (no real electrode within %.1f× spacing):", holeDepthLimit);
    lines = strings(0, 1);
    for group = interpolatedReasons.groups(:)'
        prefix = sprintf("Interpolated for %s (%d):", group.reason, numel(group.labels));
        lines = [lines; wrapLabels(prefix, group.labels(:)')]; %#ok<AGROW>
    end
    lines = [lines; wrapLabels(holePrefix, holeLabels)];
end

function lines = wrapLabels(prefix, labels)
    maxCharacters = 110;
    lines = strings(0, 1);
    if isempty(labels), return, end
    labels(1:end-1) = labels(1:end-1) + ",";
    current = prefix;
    for label = labels
        if strlength(current) + 1 + strlength(label) > maxCharacters
            lines(end + 1, 1) = current; %#ok<AGROW>
            current = label;
        else
            current = current + " " + label;
        end
    end
    lines(end + 1, 1) = current;
end
