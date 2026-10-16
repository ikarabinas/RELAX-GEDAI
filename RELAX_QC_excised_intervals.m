% RELAX_QC_excised_intervals() - Extremely-bad periods (samples) that overlap a time window,
%                                as [k x 2] seconds clipped to that window.
%
% Usage:
%   >> intervals = RELAX_QC_excised_intervals(periods, srate, windowSeconds);

function intervals = RELAX_QC_excised_intervals(periods, srate, windowSeconds)
    intervals = zeros(0, 2);
    if isempty(periods), return, end
    seconds = (periods - 1) / srate;
    overlapping = seconds(:, 2) >= windowSeconds(1) & seconds(:, 1) <= windowSeconds(2);
    intervals = [max(seconds(overlapping, 1), windowSeconds(1)), min(seconds(overlapping, 2), windowSeconds(2))];
end
