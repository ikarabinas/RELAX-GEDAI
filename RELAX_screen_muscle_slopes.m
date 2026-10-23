% RELAX_screen_muscle_slopes() - The pipeline's muscle slope of every real channel in every epoch of a cleaned file.
%
% Usage:
%   >> slopes = RELAX_screen_muscle_slopes(cleanedEEG, cfg);
%
% Inputs:
%   cleanedEEG - a cleaned file with allchan and RELAX.ListOfChannelsAfterRejections, whose
%                RELAX.Data_has_been_cleaned is 1
%   cfg        - that participant's saved RELAX_cfg (MuscleSlopeThreshold, LowPassFilter)
%
% Slopes come from RELAX_metrics_muscle on the channels that were not interpolated: the 7-75 Hz log-log
% slope of each 1 s epoch, one epoch every 0.5 s, as in the cleaned file's MuscleProportion metric.
%
% Outputs:
%   slopes.labels       - [1 x c] string, the real channels
%   slopes.values       - [c x e] slopes; empty when they cannot be computed
%   slopes.startSeconds - [1 x e] each epoch's start in the continuous data; empty when values are
%   slopes.threshold    - cfg.MuscleSlopeThreshold
%   slopes.note         - "" or why values is empty (a low-pass under 75 Hz, or the error from RELAX_metrics_muscle)

function slopes = RELAX_screen_muscle_slopes(cleanedEEG, cfg)

    slopes.labels = string(cleanedEEG.RELAX.ListOfChannelsAfterRejections(:)');
    slopes.values = [];
    slopes.startSeconds = [];
    slopes.threshold = cfg.MuscleSlopeThreshold;
    slopes.note = "";

    if cfg.LowPassFilter < 75
        slopes.note = "Muscle slopes need data kept above 75 Hz, but the low-pass filter was " + cfg.LowPassFilter + " Hz";
        return
    end
    try
        [values, startSeconds] = measureSlopes(cleanedEEG, cellstr(slopes.labels), cfg);
    catch problem
        slopes.note = "Muscle slopes could not be computed: " + string(problem.message);
        return
    end
    slopes.values = values;
    slopes.startSeconds = startSeconds;
end

function [values, startSeconds] = measureSlopes(cleanedEEG, labels, cfg)
    realChannelsEEG = pop_select(cleanedEEG, 'channel', labels);
    noAutoPool = withoutPoolAutoCreate();
    [~, epochedEEG] = RELAX_metrics_muscle(realChannelsEEG, [], cfg);
    clear noAutoPool
    values = epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels;
    startSeconds = epochStartSeconds(epochedEEG, cleanedEEG.srate);
end

function restore = withoutPoolAutoCreate()
% RELAX_metrics_muscle's parfor would otherwise start a pool on the client for a few hundred small fits.
    poolSettings = parallel.Settings;
    wasAutoCreating = poolSettings.Pool.AutoCreate;
    poolSettings.Pool.AutoCreate = false;
    restore = onCleanup(@() setAutoCreate(wasAutoCreating));
end

function setAutoCreate(value)
    poolSettings = parallel.Settings;
    poolSettings.Pool.AutoCreate = value;
end

function startSeconds = epochStartSeconds(epochedEEG, srate)
    urLatencies = [epochedEEG.urevent.latency];
    startSamples = arrayfun(@(e) urLatencies(epochedEEG.epoch(e).eventurevent(1)), 1:epochedEEG.trials);
    startSeconds = (startSamples - 1) / srate;
end
