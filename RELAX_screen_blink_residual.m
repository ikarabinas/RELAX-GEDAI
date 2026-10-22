% RELAX_screen_blink_residual() - How much of each blink is left in the cleaned data.
%
% Usage:
%   >> scores = RELAX_screen_blink_residual(cleanedEEG, peaks, blinkLabels);
%
% Inputs:
%   cleanedEEG  - cleaned continuous EEG
%   peaks       - blink peak samples (RELAX_QC_blink_peaks); each needs 0.5 s of data inside
%                 its segment before it, so find them with a half window of at least 0.5 s
%   blinkLabels - blink electrodes; those absent from cleanedEEG are ignored
%
% Outputs:
%   scores - [1 x numel(peaks)] microvolts: for each peak, the mean over the blink-electrode
%            rows of |x(peak) - mean(x(peak - 0.5 s : peak - 0.3 s))|. NaN when no blink
%            electrode is present

function scores = RELAX_screen_blink_residual(cleanedEEG, peaks, blinkLabels)
    srate = cleanedEEG.srate;
    blinkRows = find(ismember({cleanedEEG.chanlocs.labels}, blinkLabels));
    scores = nan(1, numel(peaks));
    if isempty(blinkRows), return, end

    preBlinkOffsets = -round(0.5 * srate) : -round(0.3 * srate);
    for k = 1:numel(peaks)
        atPeak = double(cleanedEEG.data(blinkRows, peaks(k)));
        beforePeak = mean(double(cleanedEEG.data(blinkRows, peaks(k) + preBlinkOffsets)), 2);
        scores(k) = mean(abs(atPeak - beforePeak));
    end
end
