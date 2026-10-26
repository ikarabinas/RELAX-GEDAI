% RELAX_QC_display_filter() - The pipeline's band-pass and resampling applied to raw-rate
%                             data, to compare a pop replay with Extremes_Rejected.
%
% Usage:
%   >> EEG = RELAX_QC_display_filter(EEG, cfg);
%
% Display only. Line noise is band-stopped at each harmonic instead of ZaplinePlus, so
% differences at the line harmonics are not pop repair.

function EEG = RELAX_QC_display_filter(EEG, cfg)
    if isinf(cfg.LowPassFilter)
        EEG = RELAX_filtbutter(EEG, cfg.HighPassFilter, [], 4, 'highpass', ...
                               cfg.causal_or_acausal_filter, cfg.RespectBoundaries);
    else
        EEG = RELAX_filtbutter(EEG, cfg.HighPassFilter, cfg.LowPassFilter, 4, 'bandpass', ...
                               cfg.causal_or_acausal_filter, cfg.RespectBoundaries);
    end
    topHarmonicHz = min(cfg.LowPassFilter, cfg.DownSample_to_X_Hz / 2) - 4;
    for f = cfg.LineNoiseFrequency : cfg.LineNoiseFrequency : topHarmonicHz
        EEG = RELAX_filtbutter(EEG, f - 4, f + 4, 4, 'bandstop', 'acausal', cfg.RespectBoundaries);
    end
    EEG = pop_resample(EEG, cfg.DownSample_to_X_Hz);
end
