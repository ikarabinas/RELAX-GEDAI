function S = GEDAI_alpha_session_summary(name, opts)
% GEDAI_ALPHA_SESSION_SUMMARY  Run broadband GEDAI on one production session and summarise it.
%
%   S = GEDAI_alpha_session_summary(name, opts)
%
%   name is the base name, e.g. 'm544_dlpfc_day2_stitched'. The pre-GEDAI input is the
%   cleaned file reduced to EEG.RELAX.ListOfChannelsAfterRejections. status is 'ok', or
%   'skipped' (reason says which leadfield or label failed); errors are left to the caller.
%
%   Fields: name, participant, status, reason, threshold, enova, labels, chanlocs, B (see
%   GEDAI_alpha_band_summary), M (counts only: nExcluded, nTinySegments, nSamplesUsed,
%   nRest, restIsWholeFile), gofInput, gofRemoved (single-dipole fits of the leading 8-13 Hz
%   topographies), tmpl and resid.
%
%   tmpl: blink and muscle topographies and spectra from the pre-wICA baseline, and their
%   correlation with removed alpha power (corrBlink, corrMuscle) and with input alpha power
%   (corrInputAlpha for blink, corrInputAlphaMuscle for muscle). Spectra use one-second
%   windows (tmpl.specFreq) because a blink window is one second long. The muscle mask is the
%   top muscleTopFraction of the 20-80 Hz trace of RELAX_QC_muscle_trace, which stands in for
%   RELAX_muscle's epoch marks. resid: mean removed power per sample inside / outside the
%   blink and muscle samples, valid rest samples only. All are NaN when the baseline does
%   not align with the cleaned file (tmpl.reason says why).

    p = inputParser;
    addRequired(p, 'name', @ischar);
    addRequired(p, 'opts', @isstruct);
    parse(p, name, opts);

    S = struct('name', name, 'participant', strtok(name, '_'), 'status', 'ok', 'reason', '');
    cleanedFull = pop_loadset('filename', [name '_RELAX.set'], 'filepath', opts.cleanedDir);
    EEG = pop_select(cleanedFull, 'channel', cleanedFull.RELAX.ListOfChannelsAfterRejections);
    labels = {EEG.chanlocs.labels};

    try
        [refCOV, Lref] = GEDAI_alpha_refcov(labels, opts.leadfieldDir, S.participant);
    catch err
        if startsWith(err.identifier, 'GEDAI_alpha:leadfield') || startsWith(err.identifier, 'RELAX_load_ppt_leadfield')
            S.status = 'skipped';
            S.reason = err.message;
            return
        end
        rethrow(err);
    end

    M = GEDAI_alpha_masks(EEG, opts);
    A = GEDAI_alpha_run(EEG, refCOV, opts);
    B = GEDAI_alpha_band_summary(A, M, EEG.srate, opts);

    S.threshold = A.threshold;
    S.enova = A.enova;
    S.labels = labels;
    S.chanlocs = EEG.chanlocs;
    S.B = B;
    S.M = countsOnly(M, EEG.srate, opts);
    S.gofInput = GEDAI_alpha_dipfit(B.vecInput, Lref);
    S.gofRemoved = GEDAI_alpha_dipfit(B.vecRemoved, Lref);

    extremes = pop_loadset('filename', [name '_Extremes_Rejected.set'], 'filepath', opts.extremesDir);
    [baseline, info] = RELAX_QC_baseline(extremes, cleanedFull);
    [S.tmpl, S.resid] = templatesAndResidual(baseline, info, labels, A, M, B, opts);
end

function counts = countsOnly(M, srate, opts)
    windowSamples = round(opts.winSec * srate);
    paddedMask = [false, M.rest(:)', false];
    runs = [find(diff(paddedMask) == 1)', find(diff(paddedMask) == -1)' - 1];
    runLengths = runs(:, 2) - runs(:, 1) + 1;
    counts.nExcluded = M.nExcluded;
    counts.nTinySegments = M.nTinySegments;
    counts.nRest = sum(M.rest);
    counts.nSamplesUsed = sum(runLengths(runLengths >= windowSamples));
    counts.restIsWholeFile = M.restIsWholeFile;
end

function [tmpl, resid] = templatesAndResidual(baseline, info, labels, A, M, B, opts)
    nch = numel(labels);
    tmplWinSec = 1;
    nBins = floor(tmplWinSec * baseline.srate / 2) + 1;
    specFreq = (0:nBins - 1) / tmplWinSec;
    specFreq = specFreq(specFreq <= opts.maxHz);

    tmpl = struct('aligned', info.aligned, 'reason', char(info.reason), 'specFreq', specFreq, ...
        'blinkTopo', nan(nch, 1), 'muscleTopo', nan(nch, 1), ...
        'blinkSpec', nan(1, numel(specFreq)), 'muscleSpec', nan(1, numel(specFreq)), ...
        'corrBlink', NaN, 'corrMuscle', NaN, 'corrInputAlpha', NaN, 'corrInputAlphaMuscle', NaN);
    resid = struct('blinkInOut', NaN, 'muscleInOut', NaN);
    if ~info.aligned
        return
    end

    [~, order] = ismember(labels, {baseline.chanlocs.labels});
    baseline.data = baseline.data(order, :);
    baseline.chanlocs = baseline.chanlocs(order);
    baseline.nbchan = nch;
    srate = baseline.srate;

    halfWindow = round(opts.blinkHalfSec * srate);
    peaks = RELAX_QC_blink_peaks(baseline, halfWindow);
    blinkMask = false(1, baseline.pnts);
    if ~isempty(peaks)
        tmpl.blinkTopo = blinkTopography(baseline.data, peaks, srate);
        for peak = peaks
            blinkMask(peak - halfWindow : peak + halfWindow) = true;
        end
        [~, blinkPower] = GEDAI_alpha_mask_psd(baseline.data, srate, blinkMask, tmplWinSec, opts.maxHz);
        blinkRows = ismember(labels, opts.blinkElectrodes);
        tmpl.blinkSpec = mean(blinkPower(blinkRows, :), 1);
    end

    [trace, channelPower] = RELAX_QC_muscle_trace(baseline, 2);
    tmpl.muscleTopo = channelPower;
    usable = M.valid & ~isnan(trace);
    muscleMask = usable & trace >= quantile(trace(usable), 1 - opts.muscleTopFraction);
    [~, musclePower] = GEDAI_alpha_mask_psd(baseline.data, srate, muscleMask, tmplWinSec, opts.maxHz);
    [~, strongest] = sort(channelPower, 'descend');
    tmpl.muscleSpec = mean(musclePower(strongest(1:opts.nTopAlpha), :), 1);

    removedAlpha = B.band.removed(:, 3);
    offAlpha = B.band.off(:, 3);
    tmpl.corrBlink = corr(abs(tmpl.blinkTopo), removedAlpha);
    tmpl.corrMuscle = corr(abs(tmpl.muscleTopo), removedAlpha);
    tmpl.corrInputAlpha = corr(abs(tmpl.blinkTopo), offAlpha);
    tmpl.corrInputAlphaMuscle = corr(abs(tmpl.muscleTopo), offAlpha);

    removedPower = sum(A.removed .^ 2, 1);
    resid.blinkInOut = inOutRatio(removedPower, blinkMask, M.rest);
    resid.muscleInOut = inOutRatio(removedPower, muscleMask, M.rest);
end

function ratio = inOutRatio(power, marked, restMask)
    inside = restMask & marked;
    outside = restMask & ~marked;
    if ~any(inside) || ~any(outside)
        ratio = NaN;
        return
    end
    ratio = mean(power(inside)) / mean(power(outside));
end

function topo = blinkTopography(data, peaks, srate)
    baselineOffsets = -round(0.5 * srate) : -round(0.3 * srate);
    topo = zeros(size(data, 1), 1);
    for peak = peaks
        atPeak = double(data(:, peak));
        beforePeak = mean(double(data(:, peak + baselineOffsets)), 2);
        topo = topo + (atPeak - beforePeak);
    end
    topo = topo / numel(peaks);
end
