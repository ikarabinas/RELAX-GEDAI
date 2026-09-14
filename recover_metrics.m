%% Recover RELAX and GEDAI cleaning metrics from individual files
% The RELAX pipeline normally aggregates cleaning metrics automatically.
% This script can be used to recover metrics in the case of an interrupted RELAX run, 
% or to save a set of user-specified RELAX and/or GEDAI metrics.
% Outputs are .mat and .csv files containing the cleaned metrics for the preprocessed files.

eeglab;

datadir = '/athena/grosenicklab/scratch/imk2003/acc_tmseeg/eeg_data/RELAX_GEDAI/RELAX_twICA_GEDAI-dlpfc';
preprocessing_label = 'RELAX_twICA_GEDAI-dlpfc';

% Point to RELAX-created cleaned data folder and define savepath
cleaned_path = fullfile(datadir, 'RELAXProcessed', 'Cleaned_Data');
savepath = fullfile(datadir, 'RELAXProcessed');

% Get list of all cleaned files
all_cleaned_files = dir(fullfile(cleaned_path, '*RELAX*.set'));
fprintf('Total cleaned files found: %d\n', length(all_cleaned_files));

% Initialize struct
CleanedMetrics_recovered = struct();

for f = 1:length(all_cleaned_files)
    fprintf('Loading file %d/%d: %s\n', f, length(all_cleaned_files), all_cleaned_files(f).name);
    EEG = pop_loadset('filename', all_cleaned_files(f).name, 'filepath', cleaned_path);
    
    % Extract filename for tracking
    [~, basename, ~] = fileparts(all_cleaned_files(f).name);
    CleanedMetrics_recovered(f).filename = basename;
    
    % RELAX_Metrics.Cleaned fields
    %CleanedMetrics_recovered(f).All_SER                                                 = EEG.RELAX_Metrics.Cleaned.All_SER;
    %CleanedMetrics_recovered(f).All_ARR                                                 = EEG.RELAX_Metrics.Cleaned.All_ARR;
    CleanedMetrics_recovered(f).MeanMuscleStrengthFromOnlySuperThresholdValues          = EEG.RELAX_Metrics.Cleaned.MeanMuscleStrengthFromOnlySuperThresholdValues;
    CleanedMetrics_recovered(f).ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel = EEG.RELAX_Metrics.Cleaned.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel;
    
    % BlinkAmplitudeRatio may be missing for files where few or no blinks were detected
    if isfield(EEG.RELAX_Metrics.Cleaned, 'BlinkAmplitudeRatio')
        frontal_idx = [18, 19, 25, 26, 31, 32, 33, 37];
        %CleanedMetrics_recovered(f).BAR = EEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio;
        CleanedMetrics_recovered(f).MeanBAR = mean(EEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio);
        CleanedMetrics_recovered(f).fBAR = mean(EEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio(frontal_idx));
    else
        %CleanedMetrics_recovered(f).BlinkAmplitudeRatio = NaN;
        CleanedMetrics_recovered(f).MeanBAR = NaN;
        CleanedMetrics_recovered(f).fBAR = NaN;
        fprintf('    WARNING: BlinkAmplitudeRatio missing in %s\n', basename);
    end
    fprintf('RELAX CleanedMetrics successfully recovered.\n');
    
    % RELAX_issues_to_check fields
    CleanedMetrics_recovered(f).ElectrodeRejectionRecommendationsMetOrExceededThreshold = EEG.RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold;
    CleanedMetrics_recovered(f).HighProportionExcludedAsExtremeOutlier                  = EEG.RELAX_issues_to_check.HighProportionExcludedAsExtremeOutlier;
    CleanedMetrics_recovered(f).NoBlinksDetected                                        = EEG.RELAX_issues_to_check.NoBlinksDetected;
    fprintf('RELAX issues_to_check metrics successfully recovered.\n');
    
    % GEDAI fields
    if isfield(EEG.etc, 'GEDAI')
        CleanedMetrics_recovered(f).SENSAI_score                = EEG.etc.GEDAI.SENSAI_score;
        CleanedMetrics_recovered(f).SENSAI_score_per_band       = EEG.etc.GEDAI.SENSAI_score_per_band;
        %CleanedMetrics_recovered(f).GEDAI_artifact_threshold_per_band = EEG.etc.GEDAI.artifact_threshold_per_band;
        CleanedMetrics_recovered(f).ENOVA                  = EEG.etc.GEDAI.mean_ENOVA;
        CleanedMetrics_recovered(f).ENOVA_per_band              = EEG.etc.GEDAI.ENOVA_per_band;
        %CleanedMetrics_recovered(f).ENOVA_per_epoch             = EEG.etc.GEDAI.ENOVA_per_epoch;
        CleanedMetrics_recovered(f).GEDAI_total_epochs                = EEG.etc.GEDAI.total_epochs;
        CleanedMetrics_recovered(f).GEDAI_epochs_percent_rejected            = EEG.etc.GEDAI.percentage_rejected;

        if isfield(EEG.etc.GEDAI, 'visualization_metrics')
            CleanedMetrics_recovered(f).ssi_before_mean              = EEG.etc.GEDAI.visualization_metrics.ssi_before_mean;
            CleanedMetrics_recovered(f).SSSI                         = EEG.etc.GEDAI.visualization_metrics.ssi_after_mean;
            CleanedMetrics_recovered(f).NSSI                          = EEG.etc.GEDAI.visualization_metrics.ssi_artifacts_mean;
            CleanedMetrics_recovered(f).lda_accuracy                  = EEG.etc.GEDAI.visualization_metrics.lda_accuracy;
            CleanedMetrics_recovered(f).SSI_silhouette                      = EEG.etc.GEDAI.visualization_metrics.signal_silhouette;
            %CleanedMetrics_recovered(f).SENSAI_ideal_power_target_db         = EEG.etc.GEDAI.visualization_metrics.ideal_power_target_db;
        end
        
        % iterate over band_names to extract scores per band sequentially
        % These categories are valid for data sampled at 500Hz
        %band_names = {'broadband', 'highgamma_190Hz', 'midgamma_94Hz', 'lowgamma_47Hz', 'beta_23Hz', 'alphabeta_12Hz', 'thetaalpha_6Hz', 'deltatheta_3Hz', 'delta_1p5Hz', 'slow_0p73Hz'};
        
        % For data sampled at 250Hz - wavelet lower frequencies
        band_names = {'broadband', '62Hz', '31Hz', '16Hz', '7p8Hz', '3p9Hz', '2Hz', '0p98Hz', '0p49Hz'};
        
        % SENSAI scores per band
        for b = 1:length(band_names)
            sensai_field_name = sprintf('GEDAI_SENSAI_score_%s', band_names{b});
            CleanedMetrics_recovered(f).(sensai_field_name) = EEG.etc.GEDAI.SENSAI_score_per_band(b);
        end
        % ENOVA score per band
        for b = 1:length(band_names)
            enova_field_name = sprintf('GEDAI_ENOVA_%s', band_names{b});
            CleanedMetrics_recovered(f).(enova_field_name) = EEG.etc.GEDAI.ENOVA_per_band(b);
        end
        % Artifact threshold per band
        %for b = 1:length(band_names)
            %artifactthresh_field_name = sprintf('GEDAI_artifact_threshold_%s', band_names{b});
            %CleanedMetrics_recovered(f).(artifactthresh_field_name) = EEG.etc.GEDAI.artifact_threshold_per_band(b);
        %end
        
        fprintf('GEDAI metrics successfully recovered.\n');
    end
    
end

% Save metrics to file as .mat and .csv
mat_filename = sprintf('%s_CleanedMetrics.mat', preprocessing_label);
csv_filename = sprintf('%s_CleanedMetrics.csv', preprocessing_label);

%save(fullfile(savepath, mat_filename), 'CleanedMetrics_recovered');

T = struct2table(CleanedMetrics_recovered);
writetable(T, fullfile(savepath, csv_filename));
fprintf('Cleaned metrics saved to file as %s and %s.\n', mat_filename, csv_filename);
