% RELAX_QC_example_windows() - The four Figure 1 windows (blink, muscle, largest pop,
%                              largest pop kept by the pipeline) for one example file.
%
% Usage:
%   >> windows = RELAX_QC_example_windows(pick, summary, paths, cfg, opts);
%
% Blink and muscle windows compare the pre-wICA baseline with the cleaned file on the
% cleaned timeline. The pop windows (RELAX_QC_pick_pops) compare the raw input with the
% pipeline's pop repair replayed on it (RELAX_QC_replay_pop_repair), both unfiltered and at
% the input rate, on the pre-excision timeline; the replay runs once per file. Each pop
% window carries .fidelity: the relative RMS difference, over the window outside later-
% excised periods, between the replay put through RELAX_QC_display_filter and
% Extremes_Rejected. A centre in pick.centres overrides the automatic one in summary
% (the pop centre applies to the largest pop only). Outputs [1 x 4] RELAX_QC_empty_window structs.

function windows = RELAX_QC_example_windows(pick, summary, paths, cfg, opts)

    extremesEEG = pop_loadset('filename', paths.extremes);
    cleanedEEG  = pop_loadset('filename', paths.cleaned);
    baseline = RELAX_QC_baseline(extremesEEG, cleanedEEG);
    shownLabels = [reshape(cfg.BlinkElectrodes(1:min(4, end)), 1, []), summary.muscleChannels];

    windows = [ ...
        comparisonWindow("blink", centreOf(pick.centres.blink, summary.centres.blink), baseline, cleanedEEG, shownLabels, opts), ...
        comparisonWindow("muscle", centreOf(pick.centres.muscle, summary.centres.muscle), baseline, cleanedEEG, shownLabels, opts), ...
        popWindows(pick, summary, extremesEEG, paths, cfg, shownLabels, opts)];
end

function w = comparisonWindow(kind, centre, before, after, labels, opts)
    w = RELAX_QC_empty_window(kind);
    if isnan(centre)
        w.note = "No " + kind + " found in this file";
        return
    end
    w = fillWindow(w, centre, before, after, labels, opts);
end

function windows = popWindows(pick, summary, extremesEEG, paths, cfg, shownLabels, opts)
    windows = [RELAX_QC_empty_window("pop"), RELAX_QC_empty_window("popKept")];
    picks = RELAX_QC_pick_pops(RELAX_QC_saved_pop_steps(extremesEEG), {extremesEEG.chanlocs.labels}, ...
                               extremesEEG.RELAX.ExtremelyBadPeriodsForDeletion, summary.inputSrate, extremesEEG.srate);
    windows(1).note = "No repaired pop in this file";
    windows(2).note = "No repaired pop outside later-excised periods";
    chosen = {picks.largest, picks.kept};
    present = ~cellfun(@isempty, chosen);
    if ~any(present), return, end

    repairedLabels = cellfun(@(p) char(p.label), chosen(present), 'UniformOutput', false);
    replayLabels = unique([repairedLabels, shownLabels], 'stable');
    inputEEG = pop_loadset('filename', paths.input);
    [unrepaired, repaired, replayInfo] = RELAX_QC_replay_pop_repair(inputEEG, extremesEEG, cfg, replayLabels);
    clear inputEEG
    filtered = RELAX_QC_display_filter(repaired, cfg);
    if filtered.pnts ~= extremesEEG.pnts
        [windows(present).note] = deal(string(sprintf('Replay has %d samples after filtering, Extremes_Rejected %d', filtered.pnts, extremesEEG.pnts)));
        return
    end

    centres = [centreOf(pick.centres.pop, NaN), NaN];
    for k = find(present)
        labels = [{char(chosen{k}.label)}, shownLabels(~strcmp(shownLabels, chosen{k}.label))];
        centre = centreOf(centres(k), chosen{k}.centreSeconds);
        windows(k) = popWindow(windows(k), chosen{k}, centre, unrepaired, repaired, filtered, extremesEEG, labels, opts);
        windows(k).note = referenceNote(replayInfo);
    end
end

function note = referenceNote(replayInfo)
    note = "Grey: raw input, unfiltered. Blue: pop repair replayed on it";
    if replayInfo.referencePopsReplayed
        note = note + " (reference pops subtracted first)";
    end
end

function c = centreOf(override, automatic)
    c = automatic;
    if ~isnan(override), c = override; end
end

function w = popWindow(w, pop, centre, unrepaired, repaired, filtered, extremesEEG, labels, opts)
    w = fillWindow(w, centre, unrepaired, repaired, labels, opts);
    firstSecond = w.t < w.t(1) + 1;
    w.before = w.before - median(w.before(:, firstSecond), 2);
    w.after = w.after - median(w.after(:, firstSecond), 2);
    w.removed = w.before - w.after;
    if isequal(centre, pop.centreSeconds) || (pop.shadeSeconds(1) >= w.t(1) && pop.shadeSeconds(2) <= w.t(end))
        w.shade = pop.shadeSeconds;
    end
    w.excised = RELAX_QC_excised_intervals(extremesEEG.RELAX.ExtremelyBadPeriodsForDeletion, extremesEEG.srate, w.t([1 end]));
    w.fidelity = replayFidelity(w, filtered, extremesEEG);
end

function fidelity = replayFidelity(w, replayed, extremesEEG)
    span = floor(w.t(1) * extremesEEG.srate) + 1 : floor(w.t(end) * extremesEEG.srate) + 1;
    seconds = (span - 1) / extremesEEG.srate;
    keep = true(size(span));
    for k = 1:size(w.excised, 1)
        keep = keep & ~(seconds >= w.excised(k, 1) & seconds <= w.excised(k, 2));
    end
    [~, replayRows] = ismember(w.channels, {replayed.chanlocs.labels});
    [~, savedRows] = ismember(w.channels, {extremesEEG.chanlocs.labels});
    usable = replayRows > 0 & savedRows > 0;
    difference = double(replayed.data(replayRows(usable), span(keep))) - double(extremesEEG.data(savedRows(usable), span(keep)));
    saved = double(extremesEEG.data(savedRows(usable), span(keep)));
    fidelity = struct('channels', {w.channels(usable)}, ...
                      'relativeRms', sqrt(sum(difference.^2, 2) ./ sum(saved.^2, 2))', ...
                      'fractionCompared', mean(keep));
end

function w = fillWindow(w, centre, before, after, labels, opts)
    present = labels(ismember(labels, {before.chanlocs.labels}) & ismember(labels, {after.chanlocs.labels}));
    [~, beforeRows] = ismember(present, {before.chanlocs.labels});
    [~, afterRows]  = ismember(present, {after.chanlocs.labels});
    half = round(opts.exampleSeconds / 2 * after.srate);
    middle = round(centre * after.srate) + 1;
    first = max(1, min(middle - half, after.pnts - 2 * half));
    span = first : min(after.pnts, first + 2 * half);

    w.channels = present;
    w.t = (span - 1) / after.srate;
    w.before = double(before.data(beforeRows, span));
    w.after = double(after.data(afterRows, span));
    w.removed = w.before - w.after;
    segments = RELAX_get_segments(after);
    starts = segments(2:end, 1);
    w.boundaries = (starts(starts > span(1) & starts <= span(end))' - 1.5) / after.srate;
    w.available = true;
end
