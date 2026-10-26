% Per-file channel rejections by screen from a finished RELAX_Wrapper_slurm.m run.
statsDir = '/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Oct26/rest_treat_joined/RELAXProcessed_broadbandGEDAI/PerParticipantStats';
files = dir(fullfile(statsDir, 'RELAXProcessingExtremeRejectionsAllParticipants_*.mat'));
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name));
    T = S.RELAXProcessingExtremeRejectionsAllParticipants;
    for r = 1:height(T)
        if isempty(T.aFileName{r}); continue; end
        m = T.MuscleBasedElectrodesToReject{r}; if iscell(m) && numel(m)==1 && iscell(m{1}); m = m{1}; end
        e = T.ExtremeDataBasedChannelToReject{r}; if iscell(e) && numel(e)==1 && iscell(e{1}); e = e{1}; end
        fprintf('%s muscle=%d recommended=%d extreme=%d retained=%d\n', string(T.aFileName{r}), ...
            numel(m), T.NumberOfMuscleContaminatedChannelsRecomendedToDelete(r), numel(e), T.RetainedChannelCount(r));
    end
end
