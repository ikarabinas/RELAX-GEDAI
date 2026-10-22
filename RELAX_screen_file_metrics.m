% RELAX_screen_file_metrics() - The four residual-noise screening metrics of one cleaned file.
%
% Usage:
%   >> [m, blinkRatioNote] = RELAX_screen_file_metrics(row, cfg, allchanLabels, keptLabels);
%
% Inputs:
%   row           - one row of T from RELAX_QC_load_tables
%   cfg           - RELAX_cfg of that file's participant
%   allchanLabels - labels of EEG.allchan, the order of the per-channel blink vector
%   keptLabels    - EEG.RELAX.ListOfChannelsAfterRejections of the cleaned file: the channels that
%                   survived the rejection screens. Not the cleaned EEG.chanlocs, which is the
%                   full montage again once the rejected channels are interpolated
%
% Outputs:
%   m - struct of scalar doubles:
%         MuscleProportion, MuscleStrength - cleaned muscle metrics from the table
%         BlinkRatio    - omitnan mean cleaned blink ratio over the blink electrodes used
%                         (RELAX_screen_blink_labels); NaN if none, or if the vector does
%                         not match allchanLabels
%         NInterpolated - numel(allchanLabels) - numel(keptLabels): the channels rejected and
%                         restored by interpolation
%   blinkRatioNote - "" when BlinkRatio is finite, otherwise why it is NaN

function [m, blinkRatioNote] = RELAX_screen_file_metrics(row, cfg, allchanLabels, keptLabels)
    m.MuscleProportion = row.TableCleanedMuscleProportion;
    m.MuscleStrength   = row.TableCleanedMuscleStrength;
    [m.BlinkRatio, blinkRatioNote] = meanBlinkRatio(row.TableCleanedBlinkRatio{1}, ...
        RELAX_screen_blink_labels(cfg, row.BlinkElectrodeFallback), allchanLabels);
    m.NInterpolated    = numel(allchanLabels) - numel(keptLabels);
end

function [ratio, note] = meanBlinkRatio(ratioPerChannel, blinkLabels, allchanLabels)
    ratio = NaN;
    note = "";
    if all(isnan(ratioPerChannel(:))) && numel(ratioPerChannel) ~= numel(allchanLabels)
        note = "No cleaned blink ratio in the tables";
    elseif numel(ratioPerChannel) ~= numel(allchanLabels)
        note = sprintf("Blink-ratio vector has %d entries, allchan has %d", numel(ratioPerChannel), numel(allchanLabels));
    elseif isempty(blinkLabels)
        note = "The blink method used no blink electrodes in this file";
    else
        isBlinkChannel = ismember(allchanLabels(:), blinkLabels);
        ratio = mean(ratioPerChannel(isBlinkChannel), 'omitnan');
        if isnan(ratio), note = "No blink electrode has a cleaned blink ratio"; end
    end
end
