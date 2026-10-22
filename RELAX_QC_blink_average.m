% RELAX_QC_blink_average() - Blink-locked average of the baseline, cleaned and removed data.
%
% Usage:
%   >> avg = RELAX_QC_blink_average(baseline, cleanedEEG, peaks, blinkLabels, halfSamples, topoHalfSeconds);
%
% Inputs:
%   baseline        - EEG struct before cleaning, sample-aligned with cleanedEEG; [] when it
%                     cannot be rebuilt (baseline and removed fields come back NaN)
%   peaks           - blink peak samples (RELAX_QC_blink_peaks)
%   blinkLabels     - electrodes averaged into the erp fields
%   halfSamples     - samples either side of each peak
%   topoHalfSeconds - the topographies average this far either side of the peak
%
% Outputs (removed = baseline - cleaned; each average has its mean over times <= -0.3 s
% subtracted):
%   avg.blinkTimes                         - [1 x 2*halfSamples+1] seconds from the peak
%   avg.erpBaseline, erpCleaned, erpRemoved - [1 x 2*halfSamples+1] blink-electrode mean
%   avg.topoBaseline, topoCleaned, topoRemoved - [nChannels x 1] mean around the peak

function avg = RELAX_QC_blink_average(baseline, cleanedEEG, peaks, blinkLabels, halfSamples, topoHalfSeconds)
    srate = cleanedEEG.srate;
    nSamples = 2 * halfSamples + 1;
    avg.blinkTimes = (-halfSamples:halfSamples) / srate;
    blinkRows = ismember({cleanedEEG.chanlocs.labels}, blinkLabels);
    preBlink = avg.blinkTimes <= -0.3;
    atPeak = abs(avg.blinkTimes) <= topoHalfSeconds;

    cleaned = zeros(cleanedEEG.nbchan, nSamples);
    for p = peaks
        cleaned = cleaned + double(cleanedEEG.data(:, p - halfSamples : p + halfSamples));
    end
    cleaned = cleaned / numel(peaks); cleaned = cleaned - mean(cleaned(:, preBlink), 2);
    avg.erpCleaned = mean(cleaned(blinkRows, :), 1);
    avg.topoCleaned = mean(cleaned(:, atPeak), 2);

    if isempty(baseline)
        avg.erpBaseline = nan(1, nSamples); avg.erpRemoved = nan(1, nSamples);
        avg.topoBaseline = nan(cleanedEEG.nbchan, 1); avg.topoRemoved = nan(cleanedEEG.nbchan, 1);
        return
    end

    [base, removed] = deal(zeros(cleanedEEG.nbchan, nSamples));
    for p = peaks
        span = p - halfSamples : p + halfSamples;
        base = base + double(baseline.data(:, span));
        removed = removed + (double(baseline.data(:, span)) - double(cleanedEEG.data(:, span)));
    end
    base = base / numel(peaks); base = base - mean(base(:, preBlink), 2);
    removed = removed / numel(peaks); removed = removed - mean(removed(:, preBlink), 2);
    avg.erpBaseline = mean(base(blinkRows, :), 1);
    avg.erpRemoved = mean(removed(blinkRows, :), 1);
    avg.topoBaseline = mean(base(:, atPeak), 2);
    avg.topoRemoved = mean(removed(:, atPeak), 2);
end
