% RELAX_QC_example_windows() - The four Figure 1 windows (blink, muscle, largest pop,
%                              largest pop kept by the pipeline) for one example file.
%
% Usage:
%   >> windows = RELAX_QC_example_windows(pick, summary, paths, cfg, opts);
%
% Blink and muscle windows compare the pre-wICA baseline with the cleaned file on the
% cleaned timeline. The pop windows (RELAX_QC_pick_pops) compare the raw input, unfiltered at
% the input rate, with the final cleaned file laid back onto the pre-excision timeline
% (RELAX_QC_cleaned_on_extremes_timeline), NaN where the pipeline deleted data. A centre in pick.centres overrides the automatic one in summary
% (the pop centre applies to the largest pop only). Outputs [1 x 4] RELAX_QC_empty_window structs.

function windows = RELAX_QC_example_windows(pick, summary, paths, cfg, opts)

    extremesEEG = pop_loadset('filename', paths.extremes);
    cleanedEEG  = pop_loadset('filename', paths.cleaned);
    baseline = RELAX_QC_baseline(extremesEEG, cleanedEEG);
    shownLabels = [reshape(cfg.BlinkElectrodes(1:min(4, end)), 1, []), summary.muscleChannels];

    windows = [ ...
        comparisonWindow("blink", centreOf(pick.centres.blink, summary.centres.blink), baseline, cleanedEEG, shownLabels, opts), ...
        comparisonWindow("muscle", centreOf(pick.centres.muscle, summary.centres.muscle), baseline, cleanedEEG, shownLabels, opts), ...
        popWindows(pick, summary, extremesEEG, cleanedEEG, paths, shownLabels, opts)];
end

function w = comparisonWindow(kind, centre, before, after, labels, opts)
    w = RELAX_QC_empty_window(kind);
    if isnan(centre)
        w.note = "No " + kind + " found in this file";
        return
    end
    w = RELAX_QC_fill_window(w, centre, before, after, labels, opts.exampleSeconds);
end

function windows = popWindows(pick, summary, extremesEEG, cleanedEEG, paths, shownLabels, opts)
    windows = [RELAX_QC_empty_window("pop"), RELAX_QC_empty_window("popKept")];
    handPickedCentre = pick.centres.pop;
    shownSeconds = zeros(0, 2);
    if ~isnan(handPickedCentre), shownSeconds = handPickedCentre + [-1 1] * opts.exampleSeconds / 2; end
    picks = RELAX_QC_pick_pops(RELAX_QC_saved_pop_steps(extremesEEG), {extremesEEG.chanlocs.labels}, ...
                               extremesEEG.RELAX.ExtremelyBadPeriodsForDeletion, summary.inputSrate, extremesEEG.srate, shownSeconds);
    windows(1).note = "No repaired pop in this file";
    windows(2).note = "No repaired pop outside later-excised periods";
    chosen = {picks.largest, picks.kept};
    present = ~cellfun(@isempty, chosen);
    if ~any(present), return, end

    repairedLabels = cellfun(@(p) char(p.label), chosen(present), 'UniformOutput', false);
    cleanedLabels = {cleanedEEG.chanlocs.labels};
    inputEEG = pop_loadset('filename', paths.input);
    inputLabels = {inputEEG.chanlocs.labels};
    replayLabels = unique([repairedLabels, shownLabels], 'stable');
    replayLabels = replayLabels(ismember(replayLabels, cleanedLabels) & ismember(replayLabels, inputLabels));
    raw = pop_select(inputEEG, 'channel', replayLabels);
    clear inputEEG

    centres = [centreOf(pick.centres.pop, NaN), NaN];
    for k = find(present)
        labels = [{char(chosen{k}.label)}, shownLabels(~strcmp(shownLabels, chosen{k}.label))];
        centre = centreOf(centres(k), chosen{k}.centreSeconds);
        windows(k) = popWindow(windows(k), chosen{k}, centre, raw, cleanedEEG, extremesEEG, labels, opts);
        windows(k).note = "Grey: raw input, unfiltered. Blue: final cleaned data";
    end
end

function c = centreOf(override, automatic)
    c = automatic;
    if ~isnan(override), c = override; end
end

function w = popWindow(w, pop, centre, raw, cleanedEEG, extremesEEG, labels, opts)
    w = RELAX_QC_fill_window(w, centre, raw, raw, labels, opts.exampleSeconds);
    firstSecond = w.t < w.t(1) + 1;
    w.before = w.before - median(w.before(:, firstSecond), 2);
    w.after = RELAX_QC_cleaned_on_extremes_timeline(cleanedEEG, extremesEEG.RELAX.ExtremelyBadPeriodsForDeletion, ...
                                                    extremesEEG.pnts, w.channels, w.t);
    w.removed = w.before - w.after;
    if isequal(centre, pop.centreSeconds) || (pop.shadeSeconds(1) >= w.t(1) && pop.shadeSeconds(2) <= w.t(end))
        w.shade = pop.shadeSeconds;
    end
    w.excised = RELAX_QC_excised_intervals(extremesEEG.RELAX.ExtremelyBadPeriodsForDeletion, extremesEEG.srate, w.t([1 end]));
end
