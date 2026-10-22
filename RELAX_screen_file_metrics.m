% RELAX_screen_file_metrics() - The four residual-noise screening metrics of one cleaned file.
%
% Usage:
%   >> m = RELAX_screen_file_metrics(row, cfg, allchanLabels, keptLabels);
%
% Inputs:
%   row           - one row of T from RELAX_QC_load_tables
%   cfg           - RELAX_cfg of that file's participant
%   allchanLabels - labels of EEG.allchan, the order of the per-channel blink vector
%   keptLabels    - labels of the channels left in the cleaned file's EEG.chanlocs
%
% Outputs:
%   m - struct of scalar doubles:
%         MuscleProportion, MuscleStrength - cleaned muscle metrics from the table
%         BlinkRatio    - omitnan mean cleaned blink ratio over the blink electrodes used
%                         (RELAX_screen_blink_labels); NaN if none, or if the vector does
%                         not match allchanLabels
%         NInterpolated - channels in allchanLabels that are not in keptLabels

function m = RELAX_screen_file_metrics(row, cfg, allchanLabels, keptLabels)
    m.MuscleProportion = row.TableCleanedMuscleProportion;
    m.MuscleStrength   = row.TableCleanedMuscleStrength;
    m.BlinkRatio       = meanBlinkRatio(row.TableCleanedBlinkRatio{1}, ...
        RELAX_screen_blink_labels(cfg, row.BlinkElectrodeFallback), allchanLabels);
    m.NInterpolated    = numel(allchanLabels) - numel(keptLabels);
end

function ratio = meanBlinkRatio(ratioPerChannel, blinkLabels, allchanLabels)
    ratio = NaN;
    if numel(ratioPerChannel) ~= numel(allchanLabels), return, end
    isBlinkChannel = ismember(allchanLabels(:), blinkLabels);
    ratio = mean(ratioPerChannel(isBlinkChannel), 'omitnan');
end
