% RELAX_QC_cleaned_on_extremes_timeline() - Cleaned-file samples laid back onto the pre-excision timeline.
%
% Usage:
%   >> data = RELAX_QC_cleaned_on_extremes_timeline(cleanedEEG, periods, extremesPnts, labels, seconds);
%
% Inputs:
%   cleanedEEG   - cleaned EEG, at the sample rate of Extremes_Rejected
%   periods      - RELAX.ExtremelyBadPeriodsForDeletion of Extremes_Rejected, [start end] samples
%   extremesPnts - number of samples of Extremes_Rejected
%   labels       - channels to return, in this order; each must be in cleanedEEG
%   seconds      - times to evaluate, on the Extremes_Rejected timeline
%
% Output is [numel(labels) x numel(seconds)], linearly interpolated, NaN where the pipeline deleted
% the data (and within one sample of such a gap).

function data = RELAX_QC_cleaned_on_extremes_timeline(cleanedEEG, periods, extremesPnts, labels, seconds)
    keep = true(1, extremesPnts);
    for k = 1:size(periods, 1), keep(periods(k, 1):periods(k, 2)) = false; end
    if sum(keep) ~= cleanedEEG.pnts
        error('RELAX_QC:timelineMismatch', 'Extremes_Rejected minus its deleted periods has %d samples, cleaned file %d', sum(keep), cleanedEEG.pnts);
    end
    cleanedIndex = cumsum(keep);

    srate = cleanedEEG.srate;
    grid = max(1, floor(seconds(1) * srate) + 1) : min(extremesPnts, ceil(seconds(end) * srate) + 1);
    [found, rows] = ismember(labels, {cleanedEEG.chanlocs.labels});
    if ~all(found), error('RELAX_QC:channelMissing', 'Not in the cleaned file: %s', strjoin(labels(~found), ', ')); end

    onGrid = NaN(numel(labels), numel(grid));
    present = keep(grid);
    onGrid(:, present) = double(cleanedEEG.data(rows, cleanedIndex(grid(present))));
    data = interp1((grid - 1) / srate, onGrid', seconds, 'linear')';
end
