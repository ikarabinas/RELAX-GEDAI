% RELAX_QC_replay_pop_repair() - Re-run the wrapper's pop repair on the raw input file, for
%                                unfiltered display beside the unrepaired data.
%
% Usage:
%   >> [unrepaired, repaired, info] = RELAX_QC_replay_pop_repair(inputEEG, extremesEEG, cfg, showLabels);
%
% The montage is the input restricted to the channels of Extremes_Rejected plus those the
% saved steps name, an approximation of what the wrapper had at that point. Reference pops
% are subtracted first when cfg asks for it and their steps were saved on the EEG; the
% reference trace is recomputed because the wrapper does not save it. Repair options and
% BasisExcludeLabels follow RELAX_Wrapper_slurm.m. Both outputs hold only showLabels (those
% present), in dataset order, still at the input rate and unfiltered.
%
% info - .referencePopsReplayed (logical) .nSteps .nMontageChannels

function [unrepaired, repaired, info] = RELAX_QC_replay_pop_repair(inputEEG, extremesEEG, cfg, showLabels)

    steps = RELAX_QC_saved_pop_steps(extremesEEG);
    inputLabels = {inputEEG.chanlocs.labels};
    wanted = union({extremesEEG.chanlocs.labels}, reshape(steps.label, 1, []));
    montage = inputLabels(ismember(inputLabels, wanted));
    EEG = pop_select(inputEEG, 'channel', montage);
    showLabels = showLabels(ismember(showLabels, montage));
    info = struct('referencePopsReplayed', false, 'nSteps', numel(steps.sample), 'nMontageChannels', EEG.nbchan);

    unrepaired = pop_select(EEG, 'channel', showLabels);

    if isfield(cfg, 'Perform_reference_pop_repair') && cfg.Perform_reference_pop_repair == 1 ...
            && isfield(extremesEEG.RELAX, 'RefPopSteps') && ~isempty(extremesEEG.RELAX.RefPopSteps.onset)
        referenceSteps = extremesEEG.RELAX.RefPopSteps;
        if istable(referenceSteps), referenceSteps = table2struct(referenceSteps, 'ToScalar', true); end
        reference = RELAX_reference_channel(EEG);
        EEG = RELAX_subtract_reference_pops(EEG, referenceSteps, reference.data, 'WindowMs', cfg.PopDetect_WindowMs);
        info.referencePopsReplayed = true;
    end

    [~, steps.channel] = ismember(steps.label, {EEG.chanlocs.labels});
    onMontage = steps.channel > 0;
    steps = structfun(@(column) column(onMontage), steps, 'UniformOutput', false);

    [EEG, repairInfo] = RELAX_repair_steps(EEG, steps, ...
        'MaxPolyOrder',       cfg.PopRepair_MaxPolyOrder, ...
        'MaxCutoffHz',        cfg.PopRepair_MaxCutoffHz, ...
        'MaxWindowSeconds',   cfg.PopRepair_MaxWindowSeconds, ...
        'DetrendThresh',      cfg.PopRepair_DetrendThresh, ...
        'DetrendIterations',  cfg.PopRepair_DetrendIterations, ...
        'AnchorMs',           cfg.PopRepair_AnchorMs, ...
        'FillMethod',         cfg.PopRepair_FillMethod, ...
        'BasisExcludeLabels', basisExcludeLabels(EEG, extremesEEG, cfg), ...
        'MinBasisChannels',   cfg.PopRepair_MinBasisChannels, ...
        'ExciseCommonMode',   cfg.PopRepair_ExciseCommonMode, ...
        'MaxExciseMs',        cfg.PopRepair_MaxExciseMs, ...
        'MinSegmentSeconds',  cfg.PopRepair_MinSegmentSeconds);
    if ~repairInfo.available
        error('RELAX_QC:replayUnavailable', 'Pop repair replay unavailable: %s', repairInfo.status);
    end
    repaired = pop_select(EEG, 'channel', showLabels);
end

function labels = basisExcludeLabels(EEG, extremesEEG, cfg)
    labels = reshape(extremesEEG.RELAX.PopScreen.degenerateChannels, 1, []);
    rejectsImpedance = isfield(cfg, 'Perform_impedance_bad_channel_rejection') && ...
                       cfg.Perform_impedance_bad_channel_rejection == 1 && isfield(EEG.chanlocs, 'impedance');
    if ~rejectsImpedance, return, end
    limit = 500;
    if isfield(cfg, 'MaxImpedanceKOhm'), limit = cfg.MaxImpedanceKOhm; end
    impedance = NaN(1, EEG.nbchan);
    for c = 1:EEG.nbchan
        value = EEG.chanlocs(c).impedance;
        if isnumeric(value) && isscalar(value), impedance(c) = double(value); end
    end
    labels = union(labels, {EEG.chanlocs(impedance > limit).labels});
end
