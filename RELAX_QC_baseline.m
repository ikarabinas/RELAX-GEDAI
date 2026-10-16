% RELAX_QC_baseline() - Rebuild what targeted wICA received, from the saved
%                       Extremes_Rejected file, aligned sample for sample with the
%                       cleaned file.
%
% Usage:
%   >> [baseline, info] = RELAX_QC_baseline(extremesEEG, cleanedEEG);
%
% Applies RELAX_Wrapper_slurm.m's own steps in its own order: excise
% RELAX.ExtremelyBadPeriodsForDeletion (eeg_eegrej), RELAX_average_rereference, then
% spherical interpolation to EEG.allchan as the wrapper does before saving. Only valid
% for configs RELAX_QC_check_config accepts.
%
% Outputs:
%   baseline - EEG struct, channels in cleanedEEG's order
%   info     - .aligned true when sample count, channel set and boundary positions all
%              match cleanedEEG; otherwise false with .reason saying which did not.

function [baseline, info] = RELAX_QC_baseline(extremesEEG, cleanedEEG)

    info = struct('aligned', false, 'reason', "");

    baseline = extremesEEG;
    periods = extremesEEG.RELAX.ExtremelyBadPeriodsForDeletion;
    if ~isempty(periods)
        baseline = eeg_eegrej(baseline, periods);
    end
    baseline = RELAX_average_rereference(baseline);
    baseline = pop_interp(baseline, baseline.allchan, 'spherical');

    if baseline.pnts ~= cleanedEEG.pnts
        info.reason = string(sprintf('pnts differ: baseline %d, cleaned %d', baseline.pnts, cleanedEEG.pnts));
        return
    end

    baselineLabels = {baseline.chanlocs.labels};
    cleanedLabels  = {cleanedEEG.chanlocs.labels};
    unmatched = setxor(baselineLabels, cleanedLabels);
    if ~isempty(unmatched)
        info.reason = "channel sets differ: " + strjoin(string(unmatched), ', ');
        return
    end
    [~, order] = ismember(cleanedLabels, baselineLabels);
    baseline.data = baseline.data(order, :);
    baseline.chanlocs = baseline.chanlocs(order);

    if ~isequal(segmentStarts(baseline), segmentStarts(cleanedEEG))
        info.reason = "boundary positions differ";
        return
    end
    info.aligned = true;
end

function starts = segmentStarts(EEG)
    segments = RELAX_get_segments(EEG);
    starts = segments(:, 1)';
end
