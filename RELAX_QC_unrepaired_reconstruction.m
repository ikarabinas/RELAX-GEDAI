% RELAX_QC_unrepaired_reconstruction() - The input file put through the pipeline's
%                                        band-pass and resampling with no pop repair, for
%                                        display beside the repaired Extremes_Rejected data.
%
% Usage:
%   >> EEG = RELAX_QC_unrepaired_reconstruction(inputFile, labels, cfg);
%
% Display only: shows what the zero-phase band-pass does to an unrepaired step. Line
% noise is band-stopped at each harmonic instead of ZaplinePlus, and neither reference
% nor per-channel pop repair is applied.

function EEG = RELAX_QC_unrepaired_reconstruction(inputFile, labels, cfg)
    EEG = pop_loadset('filename', inputFile);
    EEG = pop_select(EEG, 'channel', labels(ismember(labels, {EEG.chanlocs.labels})));
    EEG = RELAX_filtbutter(EEG, cfg.HighPassFilter, cfg.LowPassFilter, 4, 'bandpass', ...
                           cfg.causal_or_acausal_filter, cfg.RespectBoundaries);
    for f = cfg.LineNoiseFrequency : cfg.LineNoiseFrequency : cfg.LowPassFilter - 4
        EEG = RELAX_filtbutter(EEG, f - 4, f + 4, 4, 'bandstop', 'acausal', cfg.RespectBoundaries);
    end
    EEG = pop_resample(EEG, cfg.DownSample_to_X_Hz);
end
