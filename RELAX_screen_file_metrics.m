% RELAX_screen_file_metrics() - The residual-noise screening metrics of one cleaned file.
%
% Usage:
%   >> [m, note] = RELAX_screen_file_metrics(row, cfg, allchan, keptLabels, holeDepthLimit);
%
% Inputs:
%   row            - one row of T from RELAX_QC_load_tables
%   cfg            - RELAX_cfg of that file's participant
%   allchan        - EEG.allchan: struct array with labels, X, Y and Z, in the order of the
%                    per-channel blink vector
%   keptLabels     - EEG.RELAX.ListOfChannelsAfterRejections of the cleaned file: the channels that
%                    survived the rejection screens. Not the cleaned EEG.chanlocs, which is the
%                    full montage again once the rejected channels are interpolated
%   holeDepthLimit - in spacings, an interpolated electrode further than this from its nearest real
%                    electrode is listed in HoleElectrodes (RELAX_screen_hole_depth)
%
% Outputs:
%   m - struct:
%         MuscleProportion, MuscleStrength - cleaned muscle metrics from the table
%         BlinkRatio    - omitnan mean cleaned blink ratio over the blink electrodes used
%                         (RELAX_screen_blink_labels); NaN if none, or if the vector does
%                         not match allchan
%         NInterpolated - number of allchan labels not in keptLabels: the channels rejected and
%                         restored by interpolation, the set the hole is measured on
%         HoleDepth, HoleDeepestElectrode, HoleElectrodes, Spacing - RELAX_screen_hole_depth of that
%                         set (HoleElectrodes joined by spaces; HoleDepth NaN when it cannot be measured)
%         NInterpolatedLimit - floor(p * numel(allchan)), p = cfg.MaxProportionOfElectrodesThatCanBeDeleted
%                         or 0.10 when the field is absent
%         NChannels, MaxProportionDeleted - numel(allchan) and p, for the reason line
%   note - "" when BlinkRatio and HoleDepth are finite, otherwise why either is NaN, joined by "; "

function [m, note] = RELAX_screen_file_metrics(row, cfg, allchan, keptLabels, holeDepthLimit)
    allchanLabels = {allchan.labels};
    interpolatedLabels = string(allchanLabels(~ismember(string(allchanLabels), string(keptLabels))));
    proportionDeleted = 0.10;
    if isfield(cfg, 'MaxProportionOfElectrodesThatCanBeDeleted')
        proportionDeleted = cfg.MaxProportionOfElectrodesThatCanBeDeleted;
    end

    m.MuscleProportion = row.TableCleanedMuscleProportion;
    m.MuscleStrength   = row.TableCleanedMuscleStrength;
    [m.BlinkRatio, blinkNote] = meanBlinkRatio(row.TableCleanedBlinkRatio{1}, ...
        RELAX_screen_blink_labels(cfg, row.BlinkElectrodeFallback), allchanLabels);
    m.NInterpolated    = numel(interpolatedLabels);
    hole = RELAX_screen_hole_depth(allchan, interpolatedLabels, holeDepthLimit);
    m.HoleDepth            = hole.depth;
    m.HoleDeepestElectrode = hole.deepest;
    m.HoleElectrodes       = strjoin(hole.electrodes, " ");
    m.Spacing              = hole.spacing;
    m.NInterpolatedLimit   = floor(proportionDeleted * numel(allchan));
    m.NChannels            = numel(allchan);
    m.MaxProportionDeleted = proportionDeleted;
    notes = [blinkNote, hole.note];
    note = strjoin(notes(strlength(notes) > 0), "; ");
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
