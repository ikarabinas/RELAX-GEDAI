% RELAX_screen_sheet_data() - The data behind one file's review sheet: its loudest remaining
%                             muscle bursts and blinks as windows, and its blink-locked average.
%
% Usage:
%   >> sheet = RELAX_screen_sheet_data(paths, cfg, opts);
%
% Inputs:
%   paths - .input, .extremes, .cleaned (RELAX_QC_file_paths)
%   cfg   - that participant's saved RELAX_cfg
%   opts  - .windowSeconds        window length, also the least distance between window centres
%           .nWindows             windows of each kind
%           .muscleSmoothSeconds  smoothing of the channel-mean 20-80 Hz power that ranks bursts
%           .blinkHalfSeconds     half-width of the blink-locked average
%           .topoHalfSeconds      topography averages this far either side of the blink peak
%           .blinkElectrodeFallback  (optional) the file's BlinkElectrodeFallback covariate;
%                                 absent or NaN is read as 0 (RELAX_screen_blink_labels)
%
% Windows are ranked on the cleaned data. Each shows up to 4 blink electrodes plus the 2
% non-blink channels with the most 20-80 Hz power inside that window. The before-wICA trace is
% included when the baseline rebuilds and aligns; otherwise before = [] and baselineNote says why.
%
% Outputs (sheet):
%   chanlocs            - the cleaned file's
%   interpolatedLabels  - [1 x n] cell, EEG.allchan labels absent from ListOfChannelsAfterRejections
%   baselineNote        - "" when the before-wICA trace is available, otherwise the reason
%   muscle, blink       - [1 x nWindows] RELAX_QC_empty_window structs, strongest first;
%                         available = false with a note when the file has too few candidates
%   muscleRelativePower - [1 x nWindows] smoothed power at each window's centre over the file
%                         median of that trace; NaN for unavailable windows
%   blinkResidual       - [1 x nWindows] microvolts left at each window's blink
%                         (RELAX_screen_blink_residual); NaN for unavailable windows
%   blinkN              - blinks usable by the blink-locked average
%   blinkAverage        - RELAX_QC_blink_average output; [] when blinkN is 0

function sheet = RELAX_screen_sheet_data(paths, cfg, opts)

    cleanedEEG = pop_loadset('filename', paths.cleaned);
    [baseline, baselineNote] = loadBaseline(paths, cfg, cleanedEEG);
    srate = cleanedEEG.srate;
    cleanedLabels = {cleanedEEG.chanlocs.labels};

    blinkLabels = intersect(RELAX_screen_blink_labels(cfg, requestedFallback(opts)), cleanedLabels, 'stable');
    shownBlinkLabels = blinkLabels(1:min(4, numel(blinkLabels)))';

    [trace, ~, filtered] = RELAX_QC_muscle_trace(cleanedEEG, opts.muscleSmoothSeconds);

    sheet.chanlocs = cleanedEEG.chanlocs;
    sheet.interpolatedLabels = interpolatedLabels(cleanedEEG);
    sheet.baselineNote = baselineNote;
    [sheet.muscle, sheet.muscleRelativePower] = muscleWindows(cleanedEEG, baseline, filtered, trace, shownBlinkLabels, blinkLabels, opts);

    peaks = RELAX_QC_blink_peaks(cleanedEEG, round(opts.blinkHalfSeconds * srate));
    sheet.blinkN = numel(peaks);
    [sheet.blink, sheet.blinkResidual] = blinkWindows(cleanedEEG, baseline, filtered, peaks, shownBlinkLabels, blinkLabels, opts);
    sheet.blinkAverage = [];
    if ~isempty(peaks)
        sheet.blinkAverage = RELAX_QC_blink_average(baseline, cleanedEEG, peaks, blinkLabels, ...
                                                    round(opts.blinkHalfSeconds * srate), opts.topoHalfSeconds);
    end
end

function fallback = requestedFallback(opts)
    fallback = NaN;
    if isfield(opts, 'blinkElectrodeFallback'), fallback = opts.blinkElectrodeFallback; end
end

function labels = interpolatedLabels(cleanedEEG)
    labels = setdiff({cleanedEEG.allchan.labels}, cleanedEEG.RELAX.ListOfChannelsAfterRejections, 'stable');
    labels = labels(:)';
end

function [baseline, note] = loadBaseline(paths, cfg, cleanedEEG)
    baseline = [];
    note = "";
    try
        RELAX_QC_check_config(cfg);
    catch configProblem
        note = string(configProblem.message);
        return
    end
    try
        extremesEEG = pop_loadset('filename', paths.extremes);
    catch loadProblem
        note = "Extremes_Rejected file could not be loaded: " + string(loadProblem.message);
        return
    end
    [rebuilt, info] = RELAX_QC_baseline(extremesEEG, cleanedEEG);
    if info.aligned
        baseline = rebuilt;
    else
        note = info.reason;
    end
end

function [windows, relativePower] = muscleWindows(cleanedEEG, baseline, filtered, trace, shownBlinkLabels, blinkLabels, opts)
    srate = cleanedEEG.srate;
    times = (0:cleanedEEG.pnts - 1) / srate;
    centres = RELAX_screen_top_separated(times, trace, opts.nWindows, opts.windowSeconds);
    windows = repmat(RELAX_QC_empty_window("muscle"), 1, opts.nWindows);
    relativePower = nan(1, opts.nWindows);
    typicalPower = median(trace, 'omitnan');
    for k = 1:numel(centres)
        windows(k) = filledWindow("muscle", centres(k), cleanedEEG, baseline, filtered, shownBlinkLabels, blinkLabels, opts);
        relativePower(k) = trace(round(centres(k) * srate) + 1) / typicalPower;
    end
    for k = numel(centres) + 1 : opts.nWindows
        windows(k).note = shortfallNote("muscle burst", "muscle bursts", numel(centres), nnz(~isnan(trace)), opts.windowSeconds);
    end
end

function [windows, residual] = blinkWindows(cleanedEEG, baseline, filtered, peaks, shownBlinkLabels, blinkLabels, opts)
    srate = cleanedEEG.srate;
    windows = repmat(RELAX_QC_empty_window("blink"), 1, opts.nWindows);
    residual = nan(1, opts.nWindows);
    if isempty(blinkLabels)
        [windows.note] = deal("No blink electrodes in this file");
        return
    end
    peakSeconds = (peaks - 1) / srate;
    scores = RELAX_screen_blink_residual(cleanedEEG, peaks, blinkLabels);
    centres = RELAX_screen_top_separated(peakSeconds, scores, opts.nWindows, opts.windowSeconds);
    for k = 1:numel(centres)
        windows(k) = filledWindow("blink", centres(k), cleanedEEG, baseline, filtered, shownBlinkLabels, blinkLabels, opts);
        residual(k) = scores(peakSeconds == centres(k));
    end
    for k = numel(centres) + 1 : opts.nWindows
        windows(k).note = shortfallNote("blink", "blinks", numel(centres), numel(peaks), opts.windowSeconds);
    end
end

function w = filledWindow(kind, centre, cleanedEEG, baseline, filtered, shownBlinkLabels, blinkLabels, opts)
    w = RELAX_QC_fill_window(RELAX_QC_empty_window(kind), centre, baseline, cleanedEEG, shownBlinkLabels, opts.windowSeconds);
    span = round(w.t * cleanedEEG.srate) + 1;
    loudest = RELAX_screen_loudest_channels(filtered.data, {cleanedEEG.chanlocs.labels}, span, blinkLabels, 2);
    w = RELAX_QC_fill_window(RELAX_QC_empty_window(kind), centre, baseline, cleanedEEG, [shownBlinkLabels loudest], opts.windowSeconds);
end

function note = shortfallNote(singular, plural, found, candidates, separationSeconds)
    if found == 0
        note = "No " + plural + " in this file";
        return
    end
    if found == 1, noun = singular; else, noun = plural; end
    note = "Only " + found + " " + noun + " in this file";
    if candidates > found
        note = note + " at least " + separationSeconds + " s apart";
    end
end
