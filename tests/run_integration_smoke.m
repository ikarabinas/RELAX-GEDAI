function run_integration_smoke()
% Run RELAX_Wrapper_slurm end to end on the cropped fixture and assert the
% boundary-related outcomes. Validates plumbing, not cleaning quality.

    fixtureDir = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
    % Same vendored-path setup the SET_PARAMETERS entry point does after eeglab():
    % the vendored zapline-plus must shadow the EEGLAB plugin copy, and nt_detrend
    % (RELAX_repair_steps) is only available from vendor/noisetools.
    RELAX_add_vendored_zapline();
    RELAX_add_vendored_noisetools();
    assert(exist(fullfile(fixtureDir,'boundary_fixture.set'),'file') == 2, ...
        'run make_boundary_fixture first');

    RELAX_cfg = struct();
    RELAX_cfg.caploc = [];
    RELAX_cfg.myPath = fixtureDir;
    RELAX_cfg.all_data_in_1_folder_or_BIDS_format = 'folder';
    RELAX_cfg.filename = [];
    RELAX_cfg.files = {'boundary_fixture.set'};
    RELAX_cfg.dirList = dir(fullfile(fixtureDir,'boundary_fixture.set'));
    RELAX_cfg.FilesToProcess = 1;
    RELAX_cfg.ParticipantTag = 'smoke';

    % The three new flags under test
    RELAX_cfg.RespectBoundaries = 1;
    RELAX_cfg.Perform_PREP_bad_channel_rejection = 0;
    RELAX_cfg.ComputeSER_ARR = 0;
    % Electrode-pop stage at its shipped defaults (on in the SET_PARAMETERS script since
    % 2026-09-17; exercised here for plumbing - eye mask, spatial fill, excision path).
    RELAX_cfg.Perform_pop_artifact_repair = 1;

    % Everything else copied from RELAX_SET_PARAMETERS_AND_RUN_slurm.m
    RELAX_cfg.Perform_targeted_wICA = 1;
    RELAX_cfg.Run_GEDAI = 0;
    RELAX_cfg.Do_MWF_Once = 0; RELAX_cfg.Do_MWF_Twice = 0; RELAX_cfg.Do_MWF_Thrice = 0;
    RELAX_cfg.Perform_wICA_on_ICLabel = 0;
    RELAX_cfg.Perform_ICA_subtract = 0;
    RELAX_cfg.ICA_method = 'picard';
    RELAX_cfg.Report_all_ICA_info = 'yes';
    RELAX_cfg.Clean_other_comps = 'no';
    RELAX_cfg.ICLabel_thresholds = [0.5 0.5 0.5 0.5 0.5 0.5 0.5];
    RELAX_cfg.computerawmetrics = 1;
    RELAX_cfg.computecleanedmetrics = 1;
    RELAX_cfg.MWFRoundToCleanBlinks = 2;
    RELAX_cfg.LowPassFilterAt_6Hz_BeforeDetectingBlinks = 'no';
    RELAX_cfg.ProbabilityDataHasNoBlinks = 0;
    RELAX_cfg.DriftSeverityThreshold = 10;
    RELAX_cfg.ProportionWorstEpochsForDrift = 0.30;
    RELAX_cfg.ExtremeVoltageShiftThreshold = 25;
    RELAX_cfg.ExtremeAbsoluteVoltageThreshold = 1000;
    RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold = 10;
    RELAX_cfg.ExtremeSingleChannelKurtosisThreshold = 10;
    RELAX_cfg.ExtremeAllChannelKurtosisThreshold = 10;
    RELAX_cfg.ExtremeDriftSlopeThreshold = -4;
    RELAX_cfg.ExtremeBlinkShiftThreshold = 10;
    RELAX_cfg.MinimumArtifactDuration = 1200;
    RELAX_cfg.MinimumBlinkArtifactDuration = 800;
    RELAX_cfg.BlinkElectrodes = {'E18'; 'E37'; 'E19'; 'E33'; 'E32'; 'E25'; 'E31'; 'E26'};
    RELAX_cfg.HEOGLeftpattern  = ["E46", "E54", "E252", "E47", "E248", "E241", "E244"];
    RELAX_cfg.HEOGRightpattern = ["E10", "E1", "E226", "E2", "E230", "E238", "E234"];
    RELAX_cfg.BlinkMaskFocus = 150;
    RELAX_cfg.HorizontalEyeMovementType = 2;
    RELAX_cfg.HorizontalEyeMovementThreshold = 2;
    RELAX_cfg.HorizontalEyeMovementThresholdIQR = 1.5;
    RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold = 25;
    RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow = 49;
    RELAX_cfg.HorizontalEyeMovementFocus = 200;
    RELAX_cfg.LowPassFilterBeforeMWF = 'yes';
    RELAX_cfg.DownSample = 'yes';
    RELAX_cfg.DownSample_to_X_Hz = 500;
    RELAX_cfg.FilterType = 'Butterworth';
    RELAX_cfg.causal_or_acausal_filter = 'acausal';
    RELAX_cfg.HighPassFilter = 1;
    RELAX_cfg.LowPassFilter = 200;
    RELAX_cfg.NotchFilterType = 'ZaplinePlus';
    RELAX_cfg.LineNoiseFrequency = 60;
    RELAX_cfg.ElectrodesToDelete = {'E257'};
    RELAX_cfg.electrodes_2_keep_but_not_clean = {''};
    RELAX_cfg.HighPassFilter_aux_elecs = 1;
    RELAX_cfg.LowPassFilter_aux_elecs = 200;
    RELAX_cfg.KeepAllInfo = 0;
    RELAX_cfg.saveextremesrejected = 1;
    RELAX_cfg.saveround1 = 0; RELAX_cfg.saveround2 = 0; RELAX_cfg.saveround3 = 0;
    RELAX_cfg.OnlyIncludeTaskRelatedEpochs = 0;
    RELAX_cfg.MuscleSlopeThreshold = -0.31;
    RELAX_cfg.MaxProportionOfDataCanBeMarkedAsMuscle = 0.50;
    RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel = 0.50;
    RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel = 0.25;
    RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0.25;
    RELAX_cfg.InterpolateRejectedElectrodesAfterCleaning = 'yes';
    RELAX_cfg.MWFDelayPeriod_for_eye_movements = 10;
    RELAX_cfg.MWFDelayPeriod_for_muscle_artifacts = 10;
    RELAX_cfg.MWF_delay_spacing_for_eye_movements = 16;
    RELAX_cfg.MWF_delay_spacing_for_muscle_artifacts = 2;

    RELAX_Wrapper_slurm(RELAX_cfg);

    % --- Assertions on the output --------------------------------------------------
    cleanedDir = fullfile(fixtureDir, 'RELAXProcessed', 'Cleaned_Data');
    cleaned = dir(fullfile(cleanedDir, '*RELAX*.set'));
    assert(~isempty(cleaned), 'no cleaned file was written to %s', cleanedDir);

    assert_smoke_outcome(fullfile(cleanedDir, cleaned(1).name), ...
                         fullfile(fixtureDir, 'boundary_fixture.set'));

    fprintf('run_integration_smoke: PASS\n');
end
