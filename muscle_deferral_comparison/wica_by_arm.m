% wica_by_arm() - Write each arm's per-file targeted-wICA covariates to results/wica_by_arm.csv.
% The cleaned .set files don't carry RELAXProcessing_wICA, so compare_muscle_arms can't
% read them; the per-participant ProcessingStatistics_wICA tables do.
function wica_by_arm()
    myPath = '/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Oct26/rest_treat_joined/';
    keep = {'aFileName','ProportionICs_was_Muscle','ProportionVariance_was_MuscleICs', ...
            'ProportionICs_was_Brain','ProportionVariance_was_BrainICs','Proportion_artifactICs_reduced_by_wICA'};
    tables = {};
    for arm = ["Before" "After"]
        statsFiles = dir(fullfile(myPath, "RELAXProcessed_muscle" + arm, 'PerParticipantStats', 'ProcessingStatistics_wICA_*.mat'));
        for k = 1:numel(statsFiles)
            loaded = load(fullfile(statsFiles(k).folder, statsFiles(k).name));
            stats = loaded.(char(fieldnames(loaded)));
            stats = stats(:, intersect(keep, stats.Properties.VariableNames, 'stable'));
            stats.aFileName = string(stats.aFileName);
            stats = stats(stats.aFileName ~= "" & ~ismissing(stats.aFileName), :);
            stats.Arm = repmat(arm, height(stats), 1);
            tables{end+1} = stats; %#ok<AGROW>
        end
    end
    writetable(vertcat(tables{:}), fullfile(fileparts(mfilename('fullpath')), 'results', 'wica_by_arm.csv'));
end
