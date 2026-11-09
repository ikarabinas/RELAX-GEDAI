function slopes = RELAX_epoch_drift_slopes(epochedEEG)
% Channels x epochs slope of the log-frequency log-power spectrum (1-75 Hz, Hanning
% multitaper FFT), which the extreme pass compares against ExtremeDriftSlopeThreshold:
% a steep negative slope means the epoch holds little but slow drift.
    temp = eeglab2fieldtrip(epochedEEG, 'preprocessing');
    cfg              = [];
    cfg.output       = 'pow';
    cfg.channel      = {'all'};
    cfg.method       = 'mtmfft';
    cfg.taper        = 'hanning';
    cfg.foi          = 1:75;
    cfg.toi          = 0:0.05:size(temp.time{1,1},2);
    cfg.pad          = 'nextpow2';
    cfg.keeptrials   = 'yes';

    warning('ignore the following warnings about trial definition, they are not relevant to this function');
    FFTPower = ft_freqanalysis(cfg, temp);

    logFrequency = log(FFTPower.cfg.foi(1,:));
    slopes = zeros(size(FFTPower.powspctrm,2), size(FFTPower.powspctrm,1));
    for chan = 1:size(FFTPower.powspctrm,2)
        parfor trial = 1:size(FFTPower.powspctrm,1)
            powspctrm = squeeze(FFTPower.powspctrm(trial,chan,:))';
            p = polyfit(logFrequency, log(powspctrm(1,1:75)), 1);
            slopes(chan,trial) = p(1);
        end
    end
end
