% RELAX_QC_file_summary() - Stage-1 QC summary of one cleaned file.
%
% Usage:
%   >> summary = RELAX_QC_file_summary(paths, cfg, opts);
%
% Inputs:
%   paths - .input, .extremes, .cleaned: full paths of the input, Extremes_Rejected and
%           cleaned .set files
%   cfg   - that participant's saved RELAX_cfg
%   opts  - .windowSeconds, .maxHz, .blinkHalfSeconds, .topoHalfSeconds,
%           .muscleWindowSeconds (see RELAX_QC_figures)
%
% A file whose rebuilt baseline does not align with its cleaned file comes back with
% status "skipped" and the reason; everything aligned-only is left empty.

function summary = RELAX_QC_file_summary(paths, cfg, opts)

    summary = RELAX_QC_empty_summary();
    extremesEEG = pop_loadset('filename', paths.extremes);
    cleanedEEG  = pop_loadset('filename', paths.cleaned);
    [baseline, info] = RELAX_QC_baseline(extremesEEG, cleanedEEG);
    if ~info.aligned
        summary.status = "skipped";
        summary.skipReason = info.reason;
        return
    end

    srate = cleanedEEG.srate;
    labels = {cleanedEEG.chanlocs.labels};
    segments = RELAX_get_segments(cleanedEEG);
    summary.durationSeconds = cleanedEEG.pnts / srate;
    summary.nChannels = cleanedEEG.nbchan;
    summary.chanlocs = cleanedEEG.chanlocs;

    removedData = double(baseline.data) - double(cleanedEEG.data);
    [summary.freqs, p] = RELAX_QC_segment_psd(baseline.data, srate, segments, opts.windowSeconds, opts.maxHz);
    summary.psdBaseline = mean(p, 1);
    [~, p] = RELAX_QC_segment_psd(cleanedEEG.data, srate, segments, opts.windowSeconds, opts.maxHz);
    summary.psdCleaned = mean(p, 1);
    [~, p] = RELAX_QC_segment_psd(removedData, srate, segments, opts.windowSeconds, opts.maxHz);
    summary.psdRemoved = mean(p, 1);
    summary.varianceRemovedFraction = sum(var(removedData, 0, 2)) / sum(var(double(baseline.data), 0, 2));

    inputEEG = pop_loadset('filename', paths.input);
    summary.inputSrate = inputEEG.srate;
    inputEEG = pop_select(inputEEG, 'channel', intersect({inputEEG.chanlocs.labels}, labels, 'stable'));
    [inputFreqs, p] = RELAX_QC_segment_psd(inputEEG.data, inputEEG.srate, RELAX_get_segments(inputEEG), ...
                                           opts.windowSeconds, opts.maxHz);
    summary.psdInput = interp1(inputFreqs, mean(p, 1), summary.freqs);
    clear inputEEG

    summary = addBlinks(summary, baseline, cleanedEEG, cfg.BlinkElectrodes, opts);
    clear removedData

    [trace, channelPower] = RELAX_QC_muscle_trace(baseline, opts.muscleWindowSeconds);
    [summary.scores.musclePower, peak] = max(trace);
    summary.centres.muscle = (peak - 1) / srate;
    channelPower(ismember(labels, cfg.BlinkElectrodes)) = -Inf;
    [~, loudest] = sort(channelPower, 'descend');
    summary.muscleChannels = labels(loudest(1:2));

    summary = addPop(summary, extremesEEG);

    summary.metricsBaseline = RELAX_QC_artifact_metrics(baseline, cfg, cfg.BlinkElectrodes);
    summary.metricsCleaned  = RELAX_QC_artifact_metrics(cleanedEEG, cfg, cfg.BlinkElectrodes);
    summary.status = "ok";
end

function summary = addBlinks(summary, baseline, cleanedEEG, blinkLabels, opts)
    srate = cleanedEEG.srate;
    half = round(opts.blinkHalfSeconds * srate);
    peaks = RELAX_QC_blink_peaks(cleanedEEG, half);
    summary.blinkN = numel(peaks);
    summary.blinkTimes = (-half:half) / srate;
    if isempty(peaks), return, end

    avg = RELAX_QC_blink_average(baseline, cleanedEEG, peaks, blinkLabels, half, opts.topoHalfSeconds);
    summary.blinkErpBaseline = avg.erpBaseline;
    summary.blinkErpCleaned  = avg.erpCleaned;
    summary.blinkErpRemoved  = avg.erpRemoved;
    summary.topoBaseline = avg.topoBaseline;
    summary.topoCleaned  = avg.topoCleaned;
    summary.topoRemoved  = avg.topoRemoved;

    blinkRows = ismember({cleanedEEG.chanlocs.labels}, blinkLabels);
    [summary.scores.blinkAmplitude, biggest] = max(mean(abs(double(baseline.data(blinkRows, peaks))), 1));
    summary.centres.blink = (peaks(biggest) - 1) / srate;
end

function summary = addPop(summary, extremesEEG)
    steps = RELAX_QC_saved_pop_steps(extremesEEG);
    if isempty(steps), return, end
    summary.scores.popSteps = numel(steps.sample);
    % Pop steps are detected before resampling, at the input file's rate:
    picks = RELAX_QC_pick_pops(steps, {extremesEEG.chanlocs.labels}, zeros(0, 2), summary.inputSrate, extremesEEG.srate);
    summary.pop = picks.largest;
end
