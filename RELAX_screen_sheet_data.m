% RELAX_screen_sheet_data() - The data behind one file's review sheet: its most muscle-like epochs
%                             and its loudest remaining blinks as windows, and its blink-locked average.
%
% Usage:
%   >> sheet = RELAX_screen_sheet_data(paths, cfg, opts);
%
% Inputs:
%   paths - .extremes and .cleaned (the fields of RELAX_QC_file_paths this reads)
%   cfg   - that participant's saved RELAX_cfg
%   opts  - .windowSeconds        window length, also the least distance between window centres
%           .nWindows             windows of each kind
%           .blinkHalfSeconds     half-width of the blink-locked average
%           .topoHalfSeconds      topography averages this far either side of the blink peak
%           .blinkElectrodeFallback  (optional) the file's BlinkElectrodeFallback covariate;
%                                 absent or NaN is read as 0 (RELAX_screen_blink_labels)
%
% Windows are ranked on the cleaned data. Muscle windows are the epochs with the highest muscle slope
% on any real channel (RELAX_screen_muscle_slopes), shown even when none is above cfg.MuscleSlopeThreshold;
% each shows up to 4 blink electrodes plus the 2 other channels, blink electrodes or not, with the highest
% slope in that epoch, so the channel holding the epoch's highest slope is always drawn.
% Blink windows show up to 4 blink electrodes plus the 2 non-blink channels with the most 20-80 Hz power
% inside the window, and are ranked on the blink electrodes that were not interpolated. A window slides to
% stay inside its segment when the segment is long enough. The before-wICA trace is included when the
% baseline rebuilds and aligns; otherwise before = [] and baselineNote says why. When the slopes cannot be
% computed every muscle window is unavailable and says why.
%
% Outputs (sheet):
%   chanlocs            - the cleaned file's
%   interpolatedLabels  - [1 x n] cell, EEG.allchan labels absent from ListOfChannelsAfterRejections
%   interpolatedReasons - RELAX_screen_channel_reasons output for those labels
%   baselineNote        - "" when the before-wICA trace is available, otherwise the reason
%   muscle, blink       - [1 x nWindows] RELAX_QC_empty_window structs, strongest first;
%                         available = false with a note when the file has too few candidates
%   muscleSlope         - [1 x nWindows] the highest slope over channels in each window's epoch;
%                         NaN for unavailable windows
%   muscleEpochStart    - [1 x nWindows] start in seconds of each window's epoch; NaN for unavailable windows
%   muscleSlopeChannel  - [1 x nWindows] string, the channel holding that slope; "" for unavailable windows
%   muscleSlopeThreshold - cfg.MuscleSlopeThreshold
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

    filtered = RELAX_filtbutter(cleanedEEG, 20, 80, 4, 'bandpass', 'acausal', true);
    slopes = RELAX_screen_muscle_slopes(cleanedEEG, cfg);

    sheet.chanlocs = cleanedEEG.chanlocs;
    sheet.interpolatedLabels = interpolatedLabels(cleanedEEG);
    sheet.interpolatedReasons = RELAX_screen_channel_reasons(cleanedEEG);
    sheet.baselineNote = baselineNote;
    sheet.muscleSlopeThreshold = slopes.threshold;
    [sheet.muscle, sheet.muscleSlope, sheet.muscleSlopeChannel, sheet.muscleEpochStart] = muscleWindows(cleanedEEG, baseline, slopes, shownBlinkLabels, opts);

    peaks = unique(RELAX_QC_blink_peaks(cleanedEEG, round(max(opts.blinkHalfSeconds, 0.5) * srate)));
    sheet.blinkN = numel(peaks);
    rankedBlinkLabels = setdiff(blinkLabels, sheet.interpolatedLabels, 'stable');
    [sheet.blink, sheet.blinkResidual] = blinkWindows(cleanedEEG, baseline, filtered, peaks, shownBlinkLabels, blinkLabels, rankedBlinkLabels, opts);
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

function [windows, slope, slopeChannel, epochStart] = muscleWindows(cleanedEEG, baseline, slopes, shownBlinkLabels, opts)
    srate = cleanedEEG.srate;
    windows = repmat(RELAX_QC_empty_window("muscle"), 1, opts.nWindows);
    slope = nan(1, opts.nWindows);
    slopeChannel = strings(1, opts.nWindows);
    epochStart = nan(1, opts.nWindows);
    if isempty(slopes.values)
        [windows.note] = deal(slopes.note);
        return
    end
    segments = RELAX_get_segments(cleanedEEG);
    sampleWindowCentres = windowCentresInsideSegments(segments, cleanedEEG.pnts, srate, opts.windowSeconds);
    epochCentreSamples = round((slopes.startSeconds + 0.5) * srate) + 1;
    epochWindowCentres = sampleWindowCentres(epochCentreSamples);
    epochSlopes = max(slopes.values, [], 1, 'omitnan');
    centres = RELAX_screen_top_separated(epochWindowCentres, epochSlopes, opts.nWindows, opts.windowSeconds);
    for k = 1:numel(centres)
        sharingCentre = find(epochWindowCentres == centres(k));
        [slope(k), best] = max(epochSlopes(sharingCentre));
        epoch = sharingCentre(best);
        [~, loudestRow] = max(slopes.values(:, epoch));
        slopeChannel(k) = slopes.labels(loudestRow);
        epochStart(k) = slopes.startSeconds(epoch);
        windows(k) = filledMuscleWindow(centres(k), cleanedEEG, baseline, slopes, epoch, shownBlinkLabels, opts);
    end
    for k = numel(centres) + 1 : opts.nWindows
        windows(k).note = shortfallNote("muscle epoch", "muscle epochs", numel(centres), nnz(~isnan(epochSlopes)), opts.windowSeconds);
    end
end

function w = filledMuscleWindow(centre, cleanedEEG, baseline, slopes, epoch, shownBlinkLabels, opts)
    unshownRows = find(~ismember(slopes.labels, string(shownBlinkLabels)));
    [~, order] = sort(slopes.values(unshownRows, epoch), 'descend', 'MissingPlacement', 'last');
    steepest = cellstr(slopes.labels(unshownRows(order(1:min(2, numel(order))))));
    w = RELAX_QC_fill_window(RELAX_QC_empty_window("muscle"), centre, baseline, cleanedEEG, [shownBlinkLabels steepest(:)'], opts.windowSeconds);
end

function centres = windowCentresInsideSegments(segments, nSamples, srate, windowSeconds)
% Where RELAX_QC_fill_window must be centred for each sample's burst to get a window that stays inside
% the segment (2 * half + 1 samples). Separating these, not the burst times, keeps shown windows apart.
    half = round(windowSeconds / 2 * srate);
    centres = ((1:nSamples) - 1) / srate;
    for s = 1:size(segments, 1)
        first = segments(s, 1); last = segments(s, 2);
        if last - first < 2 * half, continue, end
        span = first:last;
        centres(span) = (min(max(span, first + half), last - half) - 1) / srate;
    end
end

function [windows, residual] = blinkWindows(cleanedEEG, baseline, filtered, peaks, shownBlinkLabels, blinkLabels, rankedBlinkLabels, opts)
    srate = cleanedEEG.srate;
    windows = repmat(RELAX_QC_empty_window("blink"), 1, opts.nWindows);
    residual = nan(1, opts.nWindows);
    if isempty(blinkLabels)
        [windows.note] = deal("No blink electrodes in this file");
        return
    end
    if isempty(rankedBlinkLabels)
        [windows.note] = deal("Every blink electrode in this file is interpolated");
        return
    end
    peakSeconds = (peaks - 1) / srate;
    scores = RELAX_screen_blink_residual(cleanedEEG, peaks, rankedBlinkLabels);
    centres = RELAX_screen_top_separated(peakSeconds, scores, opts.nWindows, opts.windowSeconds);
    for k = 1:numel(centres)
        windows(k) = filledWindow("blink", centres(k), cleanedEEG, baseline, filtered, shownBlinkLabels, blinkLabels, opts);
        residual(k) = scores(find(peakSeconds == centres(k), 1));
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
