% RELAX_QC_pick_pops() - The two pops Figure 1 shows for one file.
%
% Usage:
%   >> picks = RELAX_QC_pick_pops(steps, montageLabels, extremePeriods, inputSrate, extremesSrate);
%
% Inputs:
%   steps          - saved EEG.RELAX.PopScreen.steps (struct or table), [] when none were saved;
%                    .sample .onset .offset in samples at inputSrate
%   montageLabels  - channels present in Extremes_Rejected
%   extremePeriods - [k x 2] EEG.RELAX.ExtremelyBadPeriodsForDeletion, in samples at
%                    extremesSrate, on the same pre-excision timeline as steps
%
% Candidates are a file's own steps: not common-mode, not eye movement, on a channel in the
% montage. picks.largest is the candidate with the largest |amplitude|. picks.kept is the
% largest candidate whose [onset, offset] overlaps no extreme period, other than
% picks.largest itself so the two Figure 1 columns always differ. Each pick is [] when no
% step qualifies, otherwise a struct with .label .amplitude .centreSeconds .shadeSeconds.

function picks = RELAX_QC_pick_pops(steps, montageLabels, extremePeriods, inputSrate, extremesSrate)

    picks = struct('largest', [], 'kept', []);
    if isempty(steps), return, end
    if istable(steps), steps = table2struct(steps, 'ToScalar', true); end

    candidate = ~steps.commonMode(:) & ismember(steps.label(:), montageLabels);
    if isfield(steps, 'eyeMovement'), candidate = candidate & ~steps.eyeMovement(:); end
    amplitude = abs(steps.amplitude(:));

    largest = bestOf(candidate, amplitude);
    if isempty(largest), return, end
    picks.largest = describe(steps, largest, inputSrate);

    transitionSeconds = ([steps.onset(:), steps.offset(:)] - 1) / inputSrate;
    clear = ~overlapsAny(transitionSeconds, (extremePeriods - 1) / extremesSrate);
    notLargest = (1:numel(amplitude))' ~= largest;
    kept = bestOf(candidate & clear & notLargest, amplitude);
    if ~isempty(kept), picks.kept = describe(steps, kept, inputSrate); end
end

function i = bestOf(allowed, amplitude)
    i = [];
    if ~any(allowed), return, end
    amplitude(~allowed) = -Inf;
    [~, i] = max(amplitude);
end

function overlapping = overlapsAny(intervals, periods)
    overlapping = false(size(intervals, 1), 1);
    for k = 1:size(periods, 1)
        overlapping = overlapping | (intervals(:, 2) >= periods(k, 1) & intervals(:, 1) <= periods(k, 2));
    end
end

function pick = describe(steps, i, inputSrate)
    pick = struct('label', string(steps.label{i}), 'amplitude', abs(steps.amplitude(i)), ...
                  'centreSeconds', (steps.sample(i) - 1) / inputSrate, ...
                  'shadeSeconds', ([steps.onset(i), steps.offset(i)] - 1) / inputSrate);
end
