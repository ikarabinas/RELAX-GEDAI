function validate_pop_repair(outFile)
% VALIDATE_POP_REPAIR - Measure what the electrode-pop stage does on real files.
%
% Answers the question deferred by decision 7 of
% docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md: guard-band samples are
% currently left as recorded, so how much ringing survives at a repaired step site after
% the 1 Hz highpass? If it is large, guard-band interpolation gets added.
%
% Run on a compute node, never a login node:
%   srun -p sackler -c 4 --mem=32G matlab -nodisplay -nosplash -batch \
%     "addpath(pwd); validate_pop_repair('pop_validation.mat')"

    if nargin < 1; outFile = 'pop_validation.mat'; end

    eeglab nogui;
    RELAX_add_vendored_noisetools();

    dataDir = '~/TMS_response_variability/store/mdd_dlpfc/preprocessed_Aug26/tms_artifact_removed';
    files = {'m159_dlpfc_day1_treatment_20230227_123711_interTrainIntervals.set', ...
             'm497_dlpfc_day1_treatment_20240728_102919_interTrainIntervals.set', ...
             'm208_dlpfc_day1_treatment_20230717_120446_interTrainIntervals.set', ...
             'm583_dlpfc_day1_treatment_20250804_124646_interTrainIntervals.set', ...
             'm155_dlpfc_day1_treatment_20230313_012312_interTrainIntervals.set'};

    results = struct([]);
    for k = 1:numel(files)
        fprintf('\n===== %s\n', files{k});
        EEG = pop_loadset('filename', files{k}, 'filepath', dataDir);

        [steps, dInfo] = RELAX_detect_steps(EEG);
        fprintf('detected %d steps on %d channels; %d flagged as bad contacts\n', ...
            dInfo.nSteps, dInfo.nChannelsAffected, numel(dInfo.badChannels));

        rawEEG = EEG;
        [EEG, rInfo] = RELAX_repair_steps(EEG, steps);

        % --- the ringing measurement ---
        % Highpass both the raw and the repaired data exactly as the pipeline will, then
        % compare the energy in a 2 s window centred on each repaired step site against a
        % matched step-free control window on the same channel.
        hpRaw  = pop_eegfiltnew(rawEEG, 1, []);
        hpRep  = pop_eegfiltnew(EEG,    1, []);

        halfWin = round(1.0 * EEG.srate);
        ringRaw = []; ringRep = []; ctrlRep = [];
        for s = 1:numel(steps.sample)
            c = steps.channel(s); i = steps.sample(s);
            if i-halfWin < 1 || i+halfWin > EEG.pnts; continue; end
            ringRaw(end+1) = std(hpRaw.data(c, i-halfWin:i+halfWin)); %#ok<AGROW>
            ringRep(end+1) = std(hpRep.data(c, i-halfWin:i+halfWin)); %#ok<AGROW>
            j = mod(i + round(EEG.pnts/2), EEG.pnts - 2*halfWin) + halfWin + 1;
            ctrlRep(end+1) = std(hpRep.data(c, j-halfWin:j+halfWin)); %#ok<AGROW>
        end

        r = struct();
        r.file                  = files{k};
        r.detect                = dInfo;
        r.repair                = rInfo;
        r.ringingRawMedian      = median(ringRaw);
        r.ringingRepairedMedian = median(ringRep);
        r.controlMedian         = median(ctrlRep);
        r.residualRatio         = median(ringRep) / median(ctrlRep);
        results = [results, r]; %#ok<AGROW>

        fprintf(['ringing at step sites after 1 Hz HP: raw %.1f uV -> repaired %.1f uV; ' ...
                 'step-free control %.1f uV; residual ratio %.2fx\n'], ...
            r.ringingRawMedian, r.ringingRepairedMedian, r.controlMedian, r.residualRatio);
        fprintf('ladder: oa=%d o3=%d o2=%d o1=%d mean=%d fallback=%d; %.1f%% samples mean-only\n', ...
            rInfo.treatmentCounts.order3oa, rInfo.treatmentCounts.order3, ...
            rInfo.treatmentCounts.order2, rInfo.treatmentCounts.order1, ...
            rInfo.treatmentCounts.meanOnly, rInfo.treatmentCounts.fallback, ...
            100*rInfo.fractionSamplesMeanOnly);
    end

    save(outFile, 'results');
    fprintf('\nSaved %s\n', outFile);
end
