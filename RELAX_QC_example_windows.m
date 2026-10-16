% RELAX_QC_example_windows() - The three Figure 1 windows (blink, muscle, pop) for one
%                              example file.
%
% Usage:
%   >> windows = RELAX_QC_example_windows(pick, summary, paths, cfg, opts);
%
% Blink and muscle windows compare the pre-wICA baseline with the cleaned file on the
% cleaned timeline. The pop window compares RELAX_QC_unrepaired_reconstruction with
% Extremes_Rejected on the pre-excision timeline. A centre in pick.centres overrides the
% automatic one in summary. Outputs [1 x 3] RELAX_QC_empty_window structs.

function windows = RELAX_QC_example_windows(pick, summary, paths, cfg, opts)

    extremesEEG = pop_loadset('filename', paths.extremes);
    cleanedEEG  = pop_loadset('filename', paths.cleaned);
    baseline = RELAX_QC_baseline(extremesEEG, cleanedEEG);
    shownLabels = [reshape(cfg.BlinkElectrodes(1:min(4, end)), 1, []), summary.muscleChannels];

    windows = [ ...
        comparisonWindow("blink", centreOf(pick.centres.blink, summary.centres.blink), baseline, cleanedEEG, shownLabels, opts), ...
        comparisonWindow("muscle", centreOf(pick.centres.muscle, summary.centres.muscle), baseline, cleanedEEG, shownLabels, opts), ...
        popWindow(pick, summary, extremesEEG, paths, cfg, shownLabels, opts)];
end

function c = centreOf(override, automatic)
    c = automatic;
    if ~isnan(override), c = override; end
end

function w = comparisonWindow(kind, centre, before, after, labels, opts)
    w = RELAX_QC_empty_window(kind);
    if isnan(centre)
        w.note = "No " + kind + " found in this file";
        return
    end
    w = fillWindow(w, centre, before, after, labels, opts);
end

function w = popWindow(pick, summary, extremesEEG, paths, cfg, labels, opts)
    w = RELAX_QC_empty_window("pop");
    if isempty(summary.pop)
        w.note = "No repaired pop in this file";
        return
    end
    labels = [{char(summary.pop.label)}, labels(~strcmp(labels, summary.pop.label))];
    reconstruction = RELAX_QC_unrepaired_reconstruction(paths.input, labels, cfg);
    if reconstruction.pnts ~= extremesEEG.pnts
        w.note = string(sprintf('Reconstruction has %d samples, Extremes_Rejected %d', reconstruction.pnts, extremesEEG.pnts));
        return
    end
    w = fillWindow(w, centreOf(pick.centres.pop, summary.pop.centreSeconds), reconstruction, extremesEEG, labels, opts);
    if isnan(pick.centres.pop) || (summary.pop.shadeSeconds(1) >= w.t(1) && summary.pop.shadeSeconds(2) <= w.t(end))
        w.shade = summary.pop.shadeSeconds;
    end
    w.note = "Grey: input band-passed without pop repair (reconstruction, not saved pipeline data)";
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
