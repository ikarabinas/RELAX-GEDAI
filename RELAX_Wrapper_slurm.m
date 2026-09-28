%% RELAX EEG CLEANING PIPELINE, Copyright (C) (2022) Neil Bailey

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see https://www.gnu.org/licenses/.

%% RELAX_Wrapper:
function [RELAX_cfg, FileNumber, CleanedMetrics, RawMetrics, RELAXProcessingMWFStepOneAllParticipants, RELAXProcessingMWFStepTwoAllParticipants, RELAXProcessing_wICA_AllParticipants,...
        RELAXProcessing_ICA_AllParticipants, RELAXProcessingMWFStepThreeAllParticipants, RELAX_issues_to_check, RELAX_issues_to_check_2nd_run, RELAXProcessingExtremeRejectionsAllParticipants] = RELAX_Wrapper_slurm (RELAX_cfg)

% Load pre-processing statistics file for these participants if it already
% exists (note that this can cause errors if the number of variables
% inserted into the output table differs between participants, which can be
% caused by using different parameters in the preceding section):

tic;

% Turn off auto fig display
set(0, 'DefaultFigureVisible', 'off')

% Detect number of workers actually allocated to this job by SLURM
n_workers = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(n_workers) || n_workers < 1
    n_workers = str2double(getenv('SLURM_CPUS_ON_NODE'));  % set for any SLURM job
end
if isnan(n_workers) || n_workers < 1
    n_workers = feature('numcores');
end
% Give each SLURM array task its own job storage location, so that concurrent
% parpool startups don't collide on the shared default (~/.matlab/local_cluster_jobs)
if isempty(gcp('nocreate'))
    c = parcluster('local');
    scratchRoot = getenv('TMPDIR');
    if isempty(scratchRoot); scratchRoot = tempdir; end
    jobID = getenv('SLURM_JOB_ID');
    if isempty(jobID); jobID = num2str(feature('getpid')); end % fallback for non-SLURM/interactive runs
    taskID = getenv('SLURM_ARRAY_TASK_ID');
    if isempty(taskID); taskID = '0'; end % fallback for non-array jobs
    jsl = fullfile(scratchRoot, sprintf('mjs_%s_%s', jobID, taskID));
    if ~exist(jsl, 'dir'); mkdir(jsl); end
    c.JobStorageLocation = jsl;
    fprintf('Using parpool job storage location: %s\n', jsl);
    parpool(c, n_workers);
end
fprintf('Started parpool with %d workers\n', n_workers);

% Set output path
RELAX_cfg.OutputPath=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep];   % use fileseparators for increased compatability 
if ~exist(RELAX_cfg.OutputPath, 'dir'); mkdir(RELAX_cfg.OutputPath); end % make dir if not present
fprintf('Output path set to: %s\n', RELAX_cfg.OutputPath);

cd(RELAX_cfg.OutputPath);
dirList=dir('*.mat');

% Use a per-participant tag so parallel array tasks never collide on filenames
tag = 'local';
if isfield(RELAX_cfg,'ParticipantTag') && ~isempty(RELAX_cfg.ParticipantTag)
    tag = RELAX_cfg.ParticipantTag;
end

% Each SLURM task processes one participant sequentially
fprintf('Processing participant %s with SLURM task %s\n', tag, getenv('SLURM_ARRAY_TASK_ID'));

statsDir = fullfile(RELAX_cfg.OutputPath, 'PerParticipantStats');
if ~exist(statsDir, 'dir'); mkdir(statsDir); end

ProcessingStatisticsFile = fullfile(statsDir, sprintf('ProcessingStatistics_%s.mat', tag));
RawMetricsFile = fullfile(statsDir, sprintf('RawMetrics_%s.mat', tag));
CleanedMetricsFile = fullfile(statsDir, sprintf('CleanedMetrics_%s.mat', tag));
ProcessingStatistics_wICAFile = fullfile(statsDir, sprintf('ProcessingStatistics_wICA_%s.mat', tag));
RELAX_issues_to_checkFile = fullfile(statsDir, sprintf('RELAX_issues_to_check_%s.mat', tag));
RELAXProcessingExtremeRejectionsAllParticipantsFile = fullfile(statsDir, sprintf('RELAXProcessingExtremeRejectionsAllParticipants_%s.mat', tag));
% Files this run could not clean (see the catch at the end of the file loop): one row
% per file with why, so a missing cleaned file is never silent. Loaded on resume so a
% file that fails twice has one row, and a file that succeeds on a rerun loses its row.
RELAX_unprocessed_filesFile = fullfile(statsDir, sprintf('RELAX_unprocessed_files_%s.mat', tag));

statsFiles = {ProcessingStatisticsFile, RawMetricsFile, CleanedMetricsFile, ProcessingStatistics_wICAFile, RELAX_issues_to_checkFile, RELAXProcessingExtremeRejectionsAllParticipantsFile, RELAX_unprocessed_filesFile};
% Only load if this participant's own files already exist (e.g. resuming a partial run)
for x=1:numel(statsFiles)    
    if exist(statsFiles{x}, 'file')
        load(statsFiles{x});
    end
end

% for x=1:numel(dirList)
%     if  strcmp(dirList(x).name,'ProcessingStatisticsRoundOne.mat')==1
%         load('ProcessingStatisticsRoundOne.mat');
%     end
% end
% for x=1:numel(dirList)
%     if  strcmp(dirList(x).name,'ProcessingStatisticsRoundTwo.mat')==1
%         load('ProcessingStatisticsRoundTwo.mat');
%     end
% end
% for x=1:numel(dirList)
%     if  strcmp(dirList(x).name,'ProcessingStatisticsRoundThree.mat')==1
%         load('ProcessingStatisticsRoundThree.mat');
%     end
% end
% for x=1:numel(dirList)
%     if  strcmp(dirList(x).name,'RawMetrics.mat')==1
%         load('RawMetrics.mat');
%     end
% end
% for x=1:numel(dirList)
%     if  strcmp(dirList(x).name,'CleanedMetrics.mat')==1
%         load('CleanedMetrics.mat');
%     end
% end
% for x=1:numel(dirList)
%     if  strcmp(dirList(x).name,'ProcessingStatistics_wICA.mat')==1
%         load('ProcessingStatistics_wICA.mat');
%     end
% end
% for x=1:numel(dirList)
%     if  strcmp(dirList(x).name,'ProcessingStatistics_ICA.mat')==1
%         load('ProcessingStatistics_ICA.mat');
%     end
% end
% for x=1:numel(dirList)
%     if  strcmp(dirList(x).name,'RELAX_issues_to_check.mat')==1
%         load('RELAX_issues_to_check.mat');
%     end
% end
% for x=1:numel(dirList)
%     if  strcmp(dirList(x).name,'RELAXProcessingExtremeRejectionsAllParticipants.mat')==1
%         load('RELAXProcessingExtremeRejectionsAllParticipants.mat');
%     end
% end

if ~isempty(RELAX_cfg.filename)
    RELAX_cfg.FilesToProcess = 1;
    RELAX_cfg.SingleFile = 1; % 1 for single file
else
    RELAX_cfg.SingleFile = 0; % 0 for multiple files
end

WarningAboutFileNumber=0;
if size(RELAX_cfg.FilesToProcess,2) > size(RELAX_cfg.files,2)
    RELAX_cfg.FilesToProcess=RELAX_cfg.FilesToProcess(1,1):size(RELAX_cfg.files,2);
    WarningAboutFileNumber=1;
end

%% Loop selected files in the directory list:
for FileNumber=RELAX_cfg.FilesToProcess(1,1:size(RELAX_cfg.FilesToProcess,2))
    
    if RELAX_cfg.SingleFile == 0
        RELAX_cfg.filename=RELAX_cfg.files{FileNumber};
        if strcmp(RELAX_cfg.all_data_in_1_folder_or_BIDS_format,'BIDS') % v2.0.1 NWB adjusted to make BIDS compatible 8/8/2025
            RELAX_cfg.foldername=RELAX_cfg.folders{FileNumber};
        else
            RELAX_cfg.foldername=RELAX_cfg.myPath;
        end
    end

    % clearvars -except 'RELAX_cfg' 'FileNumber' 'CleanedMetrics' 'RawMetrics' 'RELAXProcessingMWFStepOneAllParticipants' 'RELAXProcessingMWFStepTwoAllParticipants' 'RELAXProcessing_wICA_AllParticipants'...
    %     'RELAXProcessing_ICA_AllParticipants' 'RELAXProcessingMWFStepThreeAllParticipants' 'Warning' 'RELAX_issues_to_check' 'RELAX_issues_to_check_2nd_run'...
    %     'RELAXProcessingExtremeRejectionsAllParticipants' 'WarningAboutFileNumber';
    clearvars -except 'RELAX_cfg' 'FileNumber' 'CleanedMetrics' 'RawMetrics' ...
            'RELAXProcessingMWFStepOneAllParticipants' 'RELAXProcessingMWFStepTwoAllParticipants' ...
            'RELAXProcessing_wICA_AllParticipants' 'RELAXProcessing_ICA_AllParticipants' ...
            'RELAXProcessingMWFStepThreeAllParticipants' 'Warning' 'RELAX_issues_to_check' ...
            'RELAX_issues_to_check_2nd_run' 'RELAXProcessingExtremeRejectionsAllParticipants' ...
            'WarningAboutFileNumber' ...
            'ProcessingStatisticsFile' 'ProcessingStatisticsRoundTwoFile' 'ProcessingStatisticsRoundThreeFile' ...
            'RawMetricsFile' 'CleanedMetricsFile' 'ProcessingStatistics_wICAFile' 'ProcessingStatistics_ICAFile' ...
            'RELAX_issues_to_checkFile' 'RELAXProcessingExtremeRejectionsAllParticipantsFile' 'statsDir' 'tag' ...
            'RELAX_unprocessed_files' 'RELAX_unprocessed_filesFile';

    % Everything from here to the end of the loop body runs under a try/catch: one file
    % that cannot be cleaned must not take the participant's remaining files down with
    % it (a MATLAB error here otherwise ends the whole SLURM array task). Nothing is
    % checkpointed until the very end of the body, so a failure part-way leaves no
    % half-written aggregate rows behind.
    try
    %% Load data (assuming the data is in EEGLAB .set format):

    %  1.1.4: fix error where PREP seems to be removed from the path after an
    % EEGLAB update:
    PrepFileLocation = which('pop_prepPipeline','-all');
    PrepFolderLocation=extractBefore(PrepFileLocation,'pop_prepPipeline.m');

    cd(RELAX_cfg.foldername); % v2.0.1 NWB adjusted to make BIDS compatible 8/8/2025
    EEG = pop_loadset(RELAX_cfg.filename);

    FileName = extractBefore(RELAX_cfg.filename,".");
    if RELAX_cfg.SingleFile == 1 % RELAX v1.1.3 NWB added to stop RELAX crashing when trying to save due to whole folder being included twice in save file
        last_slash_pos = find(RELAX_cfg.filename == '\', 1, 'last');
        FileName = extractBetween(RELAX_cfg.filename,last_slash_pos+1,".");
        FileName = FileName{1};
    end

    EEG.RELAXProcessing.aFileName=cellstr(FileName);
    EEG.RELAXProcessingExtremeRejections.aFileName=cellstr(FileName);
    
    EEG.RELAX.Data_has_been_averagerereferenced=0;
    EEG.RELAX.Data_has_been_cleaned=0;
    RELAX_cfg.ms_per_sample=(1000/EEG.srate);

    if ~isfield(RELAX_cfg,'RespectBoundaries'); RELAX_cfg.RespectBoundaries = 0; end

    % Fail fast if this file contains discontinuities and a boundary-blind stage is on:
    RELAX_check_boundary_support(EEG, RELAX_cfg);

    if ~exist([RELAX_cfg.foldername filesep 'RELAXProcessed'], 'dir')
        mkdir([RELAX_cfg.foldername filesep 'RELAXProcessed']); % make dir if not present
    end
    if ~exist([RELAX_cfg.foldername filesep 'RELAXProcessed' filesep 'RELAX_cfg'], 'dir')
        mkdir([RELAX_cfg.foldername filesep 'RELAXProcessed' filesep 'RELAX_cfg']); % per-file cfg copies, so parallel array tasks don't clobber a shared file
    end
    % Save this file's config under its own name, so parallel array tasks do not clobber
    % a shared file. fileparts strips both the directory (RELAX_cfg.filename carries a
    % full path in single-file mode) and the .set extension - saving with a .set
    % extension would write a MAT-format file that the BIDS-mode '**/*.set' glob would
    % later pick up as if it were a recording.
    [~, cfg_basename, ~] = fileparts(RELAX_cfg.filename);
    savefileone=[RELAX_cfg.foldername filesep 'RELAXProcessed' filesep 'RELAX_cfg' filesep cfg_basename '.mat'];
    save(savefileone,'RELAX_cfg')

    %% Select channels 
    % v2.0.1 NWB adjusted to prevent from crashing if RELAX_cfg.caploc variable was not provided, 8/8/2025:
    if isfield(RELAX_cfg,'caploc')
        if ~isempty(RELAX_cfg.caploc)
            EEG=pop_chanedit(EEG,  'lookup', RELAX_cfg.caploc);
        end
    end

    % v2.0.1 NWB added option to preserve additional electrodes that are
    % not cleaned, and to apply the same bad period rejections:
    %% extract non-EEG electrodes to preserve for later analysis, and filter them using different settings:
    if ~isempty(RELAX_cfg.electrodes_2_keep_but_not_clean{1})
        Non_cleaned_electrodes=pop_select(EEG,'channel',RELAX_cfg.electrodes_2_keep_but_not_clean);
    end
        
    %% Delete channels that are not relevant if present       
    EEG=pop_select(EEG,'nochannel',RELAX_cfg.ElectrodesToDelete);
    EEG = eeg_checkset( EEG );
    EEG.allchan=EEG.chanlocs; % take list of all included channels before any rejections

    % The shared channel-rejection budget is taken against the montage EEG.allchan just
    % recorded, so it must be created here, before any screen has removed a channel.
    if ~isfield(RELAX_cfg,'MaxProportionOfElectrodesThatCanBeDeleted')
        RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0.10;
    end
    EEG = RELAX_channel_budget(EEG, 'init', RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted);

    %% Reject bad electrodes on recorded impedance:
    % Runs BEFORE the pop stage and BEFORE PREP, and the order is the point. The pop
    % detector's common-mode fraction is a fraction of the channels present and its
    % spatial fill learns from the neighbours; findNoisyChannels builds a robust average
    % reference and channel-to-channel correlations - so a dead electrode skews their
    % verdicts on every OTHER channel, not just on itself. Screening first means neither
    % sees a channel already known to be a bad contact on impedance alone. (Until
    % 2026-09-18 this screen ran after the pop stage and the filters; it reads only
    % metadata, so moving it earlier needed no filtered signal.) The bridge screen runs
    % after the pop stage instead - see its own comment below for why.
    %
    % The impedance is metadata written at import by attach_mff_impedance (in the
    % TMS_response_variability repo) and rides through pop_saveset in .chanlocs, so this
    % reads nothing from the data and is unaffected by boundaries.
    %
    % Absent = off: unlike Perform_PREP_bad_channel_rejection this is new behaviour, and
    % must not switch itself on for a dataset whose files carry no impedances.
    if ~isfield(RELAX_cfg,'Perform_impedance_bad_channel_rejection')
        RELAX_cfg.Perform_impedance_bad_channel_rejection = 0;
    end
    if ~isfield(RELAX_cfg,'MaxImpedanceKOhm')
        RELAX_cfg.MaxImpedanceKOhm = 500; % matches BAD_CHAN_MAX_IMPEDANCE_KOHM upstream
    end
    if ~isfield(RELAX_cfg,'MaxProportionOfElectrodesThatCanBeDeleted')
        % RELAX_excluding_channels_and_epoching defaults this too, but only later, and the
        % cap below needs it now. Same default as there:
        RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0.10;
    end
    EEG.RELAXProcessingExtremeRejections.ImpedanceBasedChannelToReject={};
    if RELAX_cfg.Perform_impedance_bad_channel_rejection == 1
        % Spends from the shared budget created just above, so the impedance and pop
        % screens cannot between them delete more of the montage than
        % MaxProportionOfElectrodesThatCanBeDeleted allows. (Before the budget existed
        % each screen independently capped itself at that proportion, so multiple screens
        % could remove multiple times it.) PREP does NOT currently spend from this budget -
        % if RELAX_cfg.Perform_PREP_bad_channel_rejection is enabled (0 by default here),
        % total deletions across all three screens can still exceed
        % MaxProportionOfElectrodesThatCanBeDeleted. RELAX_excluding_channels_and_epoching
        % keeps its own later, independent cap - see RELAX_channel_budget.m for why.
        popBudget = EEG.RELAX.ChannelBudget;
        [ImpedanceBasedChannelToReject, ImpedanceScreenInfo] = RELAX_impedance_bad_channels( ...
            EEG.chanlocs, RELAX_cfg.MaxImpedanceKOhm, ...
            RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted, ...
            max(0, popBudget.maxRejectable - popBudget.nSpent));
        % Wrapped in an outer cell so struct2table('AsArray',true) below stores this as
        % ONE opaque per-row entry regardless of how many channels were rejected - an
        % unwrapped 1xN/Nx1 cell's own shape otherwise leaks into the table variable's
        % shape and breaks row appends. Every *ChannelToReject list is stored this way.
        EEG.RELAXProcessingExtremeRejections.ImpedanceBasedChannelToReject = {ImpedanceBasedChannelToReject};
        EEG.RELAX.ImpedanceScreen = ImpedanceScreenInfo;
        [~, EEG] = RELAX_channel_budget(EEG, 'spend', 'impedance', ImpedanceBasedChannelToReject);
        if ~ImpedanceScreenInfo.available
            warning('No usable impedances in %s - no channels rejected on impedance.', FileName);
        elseif isempty(ImpedanceBasedChannelToReject)
            fprintf('Impedance screen: no channels over %g kOhm (median %.1f, max %.1f).\n', ...
                RELAX_cfg.MaxImpedanceKOhm, ImpedanceScreenInfo.medianImpedance, ImpedanceScreenInfo.maxImpedance);
        else
            fprintf('Impedance screen: removing %d channel(s) over %g kOhm: %s\n', ...
                numel(ImpedanceBasedChannelToReject), RELAX_cfg.MaxImpedanceKOhm, ...
                strjoin(ImpedanceBasedChannelToReject', ', '));
            EEG=pop_select(EEG,'nochannel',ImpedanceBasedChannelToReject);
            EEG = eeg_checkset( EEG );
        end
    else
        fprintf('Skipping impedance-based channel rejection (RELAX_cfg.Perform_impedance_bad_channel_rejection=0).\n');
    end

    %% Detect and repair electrode pops:
    % A pop is an abrupt change in one channel's DC baseline. This runs BEFORE the notch
    % and band-pass filters, at the native sample rate, because a zero-phase highpass
    % converts a step into a multi-second bipolar excursion extending in BOTH time
    % directions - after which there is no step left to find and the damage has already
    % spread to data recorded before the pop occurred.
    %
    % Absent = off: this is new behaviour and must not switch itself on for a dataset that
    % has never been checked for pops.
    if ~isfield(RELAX_cfg,'Perform_pop_artifact_repair')
        RELAX_cfg.Perform_pop_artifact_repair = 0;
    end
    EEG.RELAXProcessingExtremeRejections.PopBasedChannelToReject={};
    if RELAX_cfg.Perform_pop_artifact_repair == 1
        PopDefaults = struct( ...
            'PopDetect_WindowMs',                            151, ...
            'PopDetect_TransitionMs',                        4, ...
            'PopDetect_ZThreshold',                          20, ...
            'PopDetect_ChannelZThreshold',                   35, ...
            'PopDetect_MinAmplitudeMicrovolts',              15, ...
            'PopDetect_RefractoryMs',                        10, ...
            'PopDetect_RefineTransition',                    true, ...
            'PopDetect_StepsPerMinuteToRejectChannel',       6, ...
            'PopDetect_MaxTransitionFractionToRejectChannel', 0.02, ...
            'PopDetect_CommonModeChannelFraction',           0.15, ...
            'PopDetect_SlopeBaselineMs',                     1000, ...
            'PopDetect_SlopeNoiseFloorMicrovolts',           0.1, ...
            'PopDetect_SlopeSpreadQuantile',                 0.85, ...
            'PopRepair_MaxPolyOrder',                        3, ...
            'PopRepair_MaxCutoffHz',                         1.0, ...
            'PopRepair_MaxWindowSeconds',                    2.0, ...
            'PopRepair_DetrendThresh',                       3, ...
            'PopRepair_DetrendIterations',                   3, ...
            'PopRepair_AnchorMs',                            3, ...
            'PopRepair_FillMethod',                          'spatial', ...
            'PopRepair_MinBasisChannels',                    16, ...
            'PopRepair_ExciseCommonMode',                    false, ...
            'PopRepair_MaxExciseMs',                         500, ...
            'PopRepair_MinSegmentSeconds',                   1.0);
        for PopField = fieldnames(PopDefaults)'
            if ~isfield(RELAX_cfg, PopField{1}); RELAX_cfg.(PopField{1}) = PopDefaults.(PopField{1}); end
        end
        % The derivative-based detector's options have no equivalent in the two-window
        % detector. A config that still sets them would otherwise run silently with the
        % new defaults while its author believes MaxVarianceRatio is doing the tuning.
        PopStaleOptions = intersect(fieldnames(RELAX_cfg), ...
            {'PopDetect_MedianFilterOrder'; 'PopDetect_SustainMs'; 'PopDetect_MaxVarianceRatio'});
        if ~isempty(PopStaleOptions)
            error('RELAX_Wrapper_slurm:stalePopOptions', ...
                ['RELAX_cfg sets %s, which belonged to the previous (derivative-based) pop ' ...
                 'detector and is ignored by RELAX_detect_steps. Remove it and tune the ' ...
                 'PopDetect_* options listed in RELAX_SET_PARAMETERS_AND_RUN_slurm.m instead.'], ...
                strjoin(PopStaleOptions, ', '));
        end

        %% Reference-electrode pops, before the per-channel detector sees them:
        % A pop on the reference contact (E257 / Cz on these files) is the same level change
        % on every channel and would otherwise be a common-mode event for the per-channel
        % detector below - excised with a new boundary, or left in if too long. It is found
        % as a step on the cross-channel median (RELAX_detect_reference_pops), checked to be
        % the same amplitude on every channel, and subtracted from every channel and from the
        % auxiliary electrodes, which share the reference. Absent = off, like the parent flag.
        % Survey and decision: docs/superpowers/reports/2026-09-21-reference-pop-survey.md.
        if ~isfield(RELAX_cfg,'Perform_reference_pop_repair')
            RELAX_cfg.Perform_reference_pop_repair = 0;
        end
        if RELAX_cfg.Perform_reference_pop_repair == 1
            RefPopDefaults = struct( ...
                'RefPop_MinChannels',            8, ...
                'RefPop_MinAmplitudeMicrovolts', 10, ...
                'RefPop_MinSignAgreement',       0.9, ...
                'RefPop_MaxAmplitudeSpread',     0.5);
            for RefPopField = fieldnames(RefPopDefaults)'
                if ~isfield(RELAX_cfg, RefPopField{1}); RELAX_cfg.(RefPopField{1}) = RefPopDefaults.(RefPopField{1}); end
            end
            [RefPopSteps, RefPopInfo] = RELAX_detect_reference_pops(EEG, ...
                'MinChannels',        RELAX_cfg.RefPop_MinChannels, ...
                'MinSignAgreement',   RELAX_cfg.RefPop_MinSignAgreement, ...
                'MaxAmplitudeSpread', RELAX_cfg.RefPop_MaxAmplitudeSpread, ...
                'Detect', {'WindowMs',               RELAX_cfg.PopDetect_WindowMs, ...
                           'TransitionMs',           RELAX_cfg.PopDetect_TransitionMs, ...
                           'ZThreshold',             RELAX_cfg.PopDetect_ZThreshold, ...
                           'ChannelZThreshold',      RELAX_cfg.PopDetect_ChannelZThreshold, ...
                           'MinAmplitudeMicrovolts', RELAX_cfg.RefPop_MinAmplitudeMicrovolts, ...
                           'RefractoryMs',           RELAX_cfg.PopDetect_RefractoryMs, ...
                           'RefineTransition',       RELAX_cfg.PopDetect_RefineTransition, ...
                           'SlopeBaselineMs',        RELAX_cfg.PopDetect_SlopeBaselineMs, ...
                           'SlopeNoiseFloorMicrovolts', RELAX_cfg.PopDetect_SlopeNoiseFloorMicrovolts, ...
                           'SlopeSpreadQuantile',    RELAX_cfg.PopDetect_SlopeSpreadQuantile});
            RefPopRepairInfo = struct('available', false, 'status', 'not run');
            if ~RefPopInfo.available
                warning('Reference-pop screen unavailable for %s (%s) - nothing subtracted.', FileName, RefPopInfo.status);
            else
                fprintf('Reference-pop screen: %d candidate(s) on the cross-channel median, %d rejected by the uniformity gate, %d merged into an overlapping neighbour, %d step(s) to subtract (median of %d channels, robust sd %.2f uV).\n', ...
                    RefPopInfo.nCandidates, RefPopInfo.nRejectedUniformity, RefPopInfo.nMerged, RefPopInfo.nSteps, ...
                    RefPopInfo.channelInfo.nChannelsUsed, RefPopInfo.channelInfo.robustSd);
                if RefPopInfo.nSteps > 0
                    [EEG, RefPopRepairInfo] = RELAX_subtract_reference_pops(EEG, RefPopSteps, ...
                        RefPopInfo.refChannel.data, 'WindowMs', RELAX_cfg.PopDetect_WindowMs);
                    if exist('Non_cleaned_electrodes', 'var')
                        Non_cleaned_electrodes = RELAX_subtract_reference_pops(Non_cleaned_electrodes, RefPopSteps, ...
                            RefPopInfo.refChannel.data, 'WindowMs', RELAX_cfg.PopDetect_WindowMs);
                    end
                    fprintf('Reference-pop repair: %d pop(s) subtracted from %d channel(s) (median |amplitude| %.0f uV, %.2f s of transition).\n', ...
                        RefPopRepairInfo.nSteps, RefPopRepairInfo.nChannelsCorrected, ...
                        RefPopRepairInfo.medianAmplitude, RefPopRepairInfo.transitionSeconds);
                end
            end
            % The virtual channel and the correction trace are pnts-long; keep the .set small.
            EEG.RELAX.RefPopScreen = rmfield(RefPopInfo, 'refChannel');
            EEG.RELAX.RefPopSteps  = RefPopSteps;
            if isfield(RefPopRepairInfo, 'correction'); RefPopRepairInfo = rmfield(RefPopRepairInfo, 'correction'); end
            EEG.RELAX.RefPopRepair = RefPopRepairInfo;
            RefPopCovariates = RELAX_refpop_covariates(RefPopInfo, RefPopRepairInfo);
            for RefPopField = fieldnames(RefPopCovariates)'
                EEG.RELAXProcessingExtremeRejections.(RefPopField{1}) = RefPopCovariates.(RefPopField{1});
            end
        else
            fprintf('Skipping reference-pop repair (RELAX_cfg.Perform_reference_pop_repair=0).\n');
        end

        [PopSteps, PopScreenInfo] = RELAX_detect_steps(EEG,...
            'WindowMs',                             RELAX_cfg.PopDetect_WindowMs, ...
            'TransitionMs',                         RELAX_cfg.PopDetect_TransitionMs, ...
            'ZThreshold',                           RELAX_cfg.PopDetect_ZThreshold, ...
            'ChannelZThreshold',                    RELAX_cfg.PopDetect_ChannelZThreshold, ...
            'MinAmplitudeMicrovolts',               RELAX_cfg.PopDetect_MinAmplitudeMicrovolts, ...
            'RefractoryMs',                         RELAX_cfg.PopDetect_RefractoryMs, ...
            'RefineTransition',                     RELAX_cfg.PopDetect_RefineTransition, ...
            'StepsPerMinuteToRejectChannel',        RELAX_cfg.PopDetect_StepsPerMinuteToRejectChannel, ...
            'MaxTransitionFractionToRejectChannel', RELAX_cfg.PopDetect_MaxTransitionFractionToRejectChannel, ...
            'CommonModeChannelFraction',            RELAX_cfg.PopDetect_CommonModeChannelFraction, ...
            'SlopeBaselineMs',                      RELAX_cfg.PopDetect_SlopeBaselineMs, ...
            'SlopeNoiseFloorMicrovolts',             RELAX_cfg.PopDetect_SlopeNoiseFloorMicrovolts, ...
            'SlopeSpreadQuantile',                  RELAX_cfg.PopDetect_SlopeSpreadQuantile);
        EEG.RELAX.PopScreen = PopScreenInfo;

        if ~PopScreenInfo.available
            warning('Pop screen unavailable for %s (%s) - no pops detected or repaired.', ...
                FileName, PopScreenInfo.status);
        else
            fprintf(['Pop screen: %d step(s) on %d channel(s) in %d event(s); ' ...
                     '%d step(s) belong to %d common-mode event(s); median rise %.1f ms.\n'], ...
                PopScreenInfo.nSteps, PopScreenInfo.nChannelsAffected, numel(PopScreenInfo.events.id), ...
                PopScreenInfo.nCommonModeSteps, PopScreenInfo.nCommonModeEvents, PopScreenInfo.medianRiseMs);

            % Captured BEFORE the Population-A filter below prunes PopSteps, so this
            % covariate shares its population with PopNSteps/PopNChannelsAffected
            % (which come from PopScreenInfo, i.e. pre-rejection). Recording one
            % covariate post-filter and its siblings pre-filter would make the set
            % silently incomparable across files.
            PopMedianStepAmplitudeAll = median(abs(PopSteps.amplitude));

            % Population B look-alikes FIRST: a horizontal saccade is a same-sign level shift
            % on a cluster of peri-ocular electrodes, indistinguishable here from a local
            % cluster of pops - and, at several per minute, it also makes those electrodes
            % look like chronically popping contacts. Repairing a saccade would take EOG
            % structure away from the later ICA, and rejecting the eye electrodes would
            % blind every later blink / HEOG stage, so clusters whose eye electrodes
            % (BlinkElectrodes plus the HEOG pair the mask resolved; two of them unless the
            % event is common-mode) have their mid-crossing inside an eye-movement period
            % (RELAX's own blink detector - and the HEOG detector if PopEye_IncludeHEOG -
            % run on a temporary copy of the still-complete montage), and none of whose
            % members exceeds PopEye_MaxAmplitudeMicrovolts, are flagged before anything
            % is rejected or repaired. Steps on the raw timeline: this runs before excision.
            if ~isfield(RELAX_cfg,'BlinkElectrodes')
                RELAX_cfg.BlinkElectrodes = {'E18'; 'E37'; 'E19'; 'E33'; 'E32'; 'E25'; 'E31'; 'E26'};
            end
            if ~isfield(RELAX_cfg,'PopEye_MaxAmplitudeMicrovolts'); RELAX_cfg.PopEye_MaxAmplitudeMicrovolts = 500; end
            [PopEyeMask, PopEyeMaskInfo] = RELAX_pop_eye_mask(EEG, RELAX_cfg);
            PopEyeLabels = union(reshape(RELAX_cfg.BlinkElectrodes, [], 1), reshape(PopEyeMaskInfo.heogChannels, [], 1));
            [PopSteps, PopEyeInfo] = RELAX_flag_eye_steps(PopSteps, PopEyeMask, PopEyeLabels, ...
                'MaxAmplitudeMicrovolts', RELAX_cfg.PopEye_MaxAmplitudeMicrovolts);
            EEG.RELAX.PopScreen.eyeMask = PopEyeMaskInfo;
            EEG.RELAX.PopScreen.eyeLabels  = PopEyeLabels;
            EEG.RELAX.PopScreen.nEyeSteps  = PopEyeInfo.nEyeSteps;
            EEG.RELAX.PopScreen.nEyeEvents = PopEyeInfo.nEyeEvents;
            EEG.RELAX.PopScreen.nEyeEventsOverAmplitude = PopEyeInfo.nEventsOverAmplitude;
            fprintf(['Pop screen: %d step(s) in %d event(s) flagged as eye movements and left ' ...
                     'alone (blink mask %.1f%%, HEOG mask %.1f%% of samples, HEOG %s; %d event(s) ' ...
                     'kept for repair by the %g uV ceiling).\n'], ...
                PopEyeInfo.nEyeSteps, PopEyeInfo.nEyeEvents, ...
                100*PopEyeMaskInfo.blinkFraction, 100*PopEyeMaskInfo.heogFraction, PopEyeMaskInfo.heogStatus, ...
                PopEyeInfo.nEventsOverAmplitude, RELAX_cfg.PopEye_MaxAmplitudeMicrovolts);

            % Population A: chronically bad contacts, ranked over the channel's own steps
            % (neither common-mode nor eye-movement); channels constant within every
            % segment (dead / railed contacts) come first as 'flat'. Repair is meaningless
            % on a channel with no stable baseline, so it is rejected rather than repaired.
            % badChannels arrives ordered worst-first, which is what the budget's truncation
            % needs. The eye electrodes are exempt: when one of them trips a criterion it is
            % far more likely to be showing eye movement the mask missed than to be dead
            % (m407: HEOG mask 1.1 %, all 14 periocular electrodes over the rate limit), and
            % the blink / HEOG stages below cannot run without them - so its steps are left
            % alone instead (like eye-movement steps) and it is counted, not removed.
            PopScreenInfo = RELAX_rank_pop_channels(PopSteps, PopScreenInfo);
            [PopSteps, PopScreenInfo] = RELAX_protect_eye_channels(PopSteps, PopScreenInfo, PopEyeLabels);
            if ~isempty(PopScreenInfo.protectedEyeChannels)
                fprintf('Pop screen: %d eye electrode(s) over a bad-contact limit kept and left alone: %s\n', ...
                    numel(PopScreenInfo.protectedEyeChannels), strjoin(PopScreenInfo.protectedEyeChannels', ', '));
            end
            % the re-ranked lists supersede the detector's own on the stored info struct
            EEG.RELAX.PopScreen.badChannels       = PopScreenInfo.badChannels;
            EEG.RELAX.PopScreen.badChannelReason  = PopScreenInfo.badChannelReason;
            EEG.RELAX.PopScreen.stepRatePerMinute = PopScreenInfo.stepRatePerMinute;
            EEG.RELAX.PopScreen.transitionFraction = PopScreenInfo.transitionFraction;
            EEG.RELAX.PopScreen.protectedEyeChannels = PopScreenInfo.protectedEyeChannels;
            [PopChannelToReject, EEG] = RELAX_channel_budget(EEG, 'spend', 'pop', PopScreenInfo.badChannels);
            % Wrapped in an outer cell - see the ImpedanceBasedChannelToReject comment above.
            % (Kept here in full for the record: struct2table('AsArray',true) below stores
            % this as ONE opaque per-row entry regardless of how many channels were rejected;
            % an unwrapped 1xN/Nx1 cell's own shape otherwise leaks into the table variable's
            % per-row width, and two files rejecting a different number of channels then
            % crash RELAX_append_table_row's T(RowNumber,:)=NewRow with a table-variable
            % size mismatch. See tests/test_RELAX_append_table_row.m for the convention.
            EEG.RELAXProcessingExtremeRejections.PopBasedChannelToReject = {PopChannelToReject};
            % Which of the two bad-contact criteria bound, over the channels actually
            % removed (the budget may have truncated the detector's list).
            PopRejectReason = PopScreenInfo.badChannelReason(ismember(PopScreenInfo.badChannels, PopChannelToReject));
            PopNRejectedByRate     = sum(contains(PopRejectReason, 'steps/min'));
            PopNRejectedByFraction = sum(contains(PopRejectReason, 'transition'));
            PopNRejectedByFlat     = sum(contains(PopRejectReason, 'flat'));
            if ~isempty(PopChannelToReject)
                fprintf('Pop screen: removing %d chronically popping channel(s): %s\n', ...
                    numel(PopChannelToReject), strjoin(PopChannelToReject, ', '));
                keep = ~ismember(PopSteps.label, PopChannelToReject);
                PopSteps = structfun(@(f) f(keep), PopSteps, 'UniformOutput', false);
                EEG = pop_select(EEG,'nochannel',PopChannelToReject);
                EEG = eeg_checkset( EEG );
                % RELAX_detect_steps returns channel INDICES; re-map them after pop_select
                [~, PopSteps.channel] = ismember(PopSteps.label, {EEG.chanlocs.labels});
            end

            % Channels the spatial fill must not learn from: flat channels, and any contact
            % over the impedance limit that the impedance screen above had to leave in
            % place because the shared budget was already spent.
            PopBasisExclude = reshape(PopScreenInfo.degenerateChannels, 1, []);
            if isfield(RELAX_cfg,'Perform_impedance_bad_channel_rejection') && ...
                    RELAX_cfg.Perform_impedance_bad_channel_rejection == 1 && isfield(EEG.chanlocs,'impedance')
                if ~isfield(RELAX_cfg,'MaxImpedanceKOhm'); RELAX_cfg.MaxImpedanceKOhm = 500; end
                PopImpedance = NaN(1, EEG.nbchan);
                for PopCh = 1:EEG.nbchan
                    if isnumeric(EEG.chanlocs(PopCh).impedance) && isscalar(EEG.chanlocs(PopCh).impedance)
                        PopImpedance(PopCh) = double(EEG.chanlocs(PopCh).impedance);
                    end
                end
                PopBasisExclude = union(PopBasisExclude, {EEG.chanlocs(PopImpedance > RELAX_cfg.MaxImpedanceKOhm).labels});
            end

            [EEG, PopRepairInfo] = RELAX_repair_steps(EEG, PopSteps, ...
                'MaxPolyOrder',       RELAX_cfg.PopRepair_MaxPolyOrder, ...
                'MaxCutoffHz',        RELAX_cfg.PopRepair_MaxCutoffHz, ...
                'MaxWindowSeconds',   RELAX_cfg.PopRepair_MaxWindowSeconds, ...
                'DetrendThresh',      RELAX_cfg.PopRepair_DetrendThresh, ...
                'DetrendIterations',  RELAX_cfg.PopRepair_DetrendIterations, ...
                'AnchorMs',           RELAX_cfg.PopRepair_AnchorMs, ...
                'FillMethod',         RELAX_cfg.PopRepair_FillMethod, ...
                'BasisExcludeLabels', PopBasisExclude, ...
                'MinBasisChannels',   RELAX_cfg.PopRepair_MinBasisChannels, ...
                'ExciseCommonMode',   RELAX_cfg.PopRepair_ExciseCommonMode, ...
                'MaxExciseMs',        RELAX_cfg.PopRepair_MaxExciseMs, ...
                'MinSegmentSeconds',  RELAX_cfg.PopRepair_MinSegmentSeconds);
            EEG.RELAX.PopRepair = PopRepairInfo;
            % Kept for auditing which intervals were touched. Sample indices are on the
            % PRE-excision timeline: once common-mode events are cut out (PopRepairInfo.
            % excised, also pre-excision), a sample after the first cut no longer lands
            % where steps.sample says.
            EEG.RELAX.PopScreen.steps = PopSteps;
            if ~PopRepairInfo.available
                warning('Pop repair unavailable for %s (%s) - steps detected but NOT repaired.', ...
                    FileName, PopRepairInfo.status);
            else
                fprintf(['Pop repair: %d interval(s) filled (%d spatial, %d linear; median %.0f ms; ' ...
                         '%.2f%% of samples, worst channel %.2f%%); %d common-mode event(s) excised ' ...
                         '(%.2f s), %d too long to excise; %d sub-segment(s) detrended, %d fallback(s).\n'], ...
                    PopRepairInfo.nIntervals, PopRepairInfo.nFillsSpatial, PopRepairInfo.nFillsLinear, ...
                    PopRepairInfo.medianTransitionMs, 100*PopRepairInfo.fractionSamplesInterpolated, ...
                    100*PopRepairInfo.maxChannelFractionInterpolated, PopRepairInfo.nExcised, ...
                    PopRepairInfo.excisedSeconds, PopRepairInfo.nCommonModeTooLong, ...
                    PopRepairInfo.nSubSegments, PopRepairInfo.treatmentCounts.fallback);
            end

            % Per-file covariates. RELAX thresholds are computed from within-recording
            % statistics, so a poppy session gets cleaned harder than a quiet one in a way
            % that tracks its electrode quality. A downstream model predicting treatment
            % response from session-to-session EEG change could learn the cleaning instead
            % of the neural signal, so how hard this stage worked has to be recorded per
            % file - see CLAUDE.md, "Saving covariates for downstream ML classification
            % tasks". Scalars and cells only: this struct is struct2table'd below.
            EEG.RELAXProcessingExtremeRejections.PopNSteps                = PopScreenInfo.nSteps;
            EEG.RELAXProcessingExtremeRejections.PopNChannelsAffected     = PopScreenInfo.nChannelsAffected;
            EEG.RELAXProcessingExtremeRejections.PopNEvents               = numel(PopScreenInfo.events.id);
            EEG.RELAXProcessingExtremeRejections.PopMedianStepAmplitude   = PopMedianStepAmplitudeAll;
            % Steps actually handed to the repair, i.e. after Population-A channels
            % were removed. Equal to PopMedianStepAmplitude when nothing was rejected.
            EEG.RELAXProcessingExtremeRejections.PopMedianRepairedStepAmplitude = median(abs(PopSteps.amplitude));
            EEG.RELAXProcessingExtremeRejections.PopMedianRiseMs          = PopScreenInfo.medianRiseMs;
            EEG.RELAXProcessingExtremeRejections.PopNChannelsRejected     = numel(PopChannelToReject);
            EEG.RELAXProcessingExtremeRejections.PopNChannelsRejectedByRate     = PopNRejectedByRate;
            EEG.RELAXProcessingExtremeRejections.PopNChannelsRejectedByFraction = PopNRejectedByFraction;
            EEG.RELAXProcessingExtremeRejections.PopNChannelsRejectedByFlat     = PopNRejectedByFlat;
            EEG.RELAXProcessingExtremeRejections.PopNEyeElectrodesProtected     = numel(PopScreenInfo.protectedEyeChannels);
            EEG.RELAXProcessingExtremeRejections.PopNCommonModeSteps      = PopScreenInfo.nCommonModeSteps;
            EEG.RELAXProcessingExtremeRejections.PopNCommonModeEvents     = PopScreenInfo.nCommonModeEvents;
            EEG.RELAXProcessingExtremeRejections.PopDegenerateChannelCount = numel(PopScreenInfo.degenerateChannels);
            EEG.RELAXProcessingExtremeRejections.PopNEyeMovementStepsSkipped = PopEyeInfo.nEyeSteps;
            EEG.RELAXProcessingExtremeRejections.PopNEyeEventsOverAmplitude  = PopEyeInfo.nEventsOverAmplitude;
            EEG.RELAXProcessingExtremeRejections.PopEyeMaskFraction       = PopEyeMaskInfo.maskFraction;
            if PopRepairInfo.available
                EEG.RELAXProcessingExtremeRejections.PopNIntervalsFilled        = PopRepairInfo.nIntervals;
                EEG.RELAXProcessingExtremeRejections.PopFractionSamplesInterpolated    = PopRepairInfo.fractionSamplesInterpolated;
                EEG.RELAXProcessingExtremeRejections.PopMaxChannelFractionInterpolated = PopRepairInfo.maxChannelFractionInterpolated;
                EEG.RELAXProcessingExtremeRejections.PopMedianTransitionMs      = PopRepairInfo.medianTransitionMs;
                EEG.RELAXProcessingExtremeRejections.PopNFillsSpatial           = PopRepairInfo.nFillsSpatial;
                EEG.RELAXProcessingExtremeRejections.PopNFillsLinear            = PopRepairInfo.nFillsLinear;
                EEG.RELAXProcessingExtremeRejections.PopNCommonModeExcised      = PopRepairInfo.nExcised;
                EEG.RELAXProcessingExtremeRejections.PopExcisedSeconds          = PopRepairInfo.excisedSeconds;
                EEG.RELAXProcessingExtremeRejections.PopNCommonModeTooLong      = PopRepairInfo.nCommonModeTooLong;
                EEG.RELAXProcessingExtremeRejections.PopMedianSubSegmentDuration = PopRepairInfo.medianSubSegmentDuration;
                EEG.RELAXProcessingExtremeRejections.PopFractionSamplesMeanOnly  = PopRepairInfo.fractionSamplesMeanOnly;
                EEG.RELAXProcessingExtremeRejections.PopMedianLowFreqVarianceRemoved = median(PopRepairInfo.lowFreqVarianceRemoved);
                EEG.RELAXProcessingExtremeRejections.PopNOrder3OA  = PopRepairInfo.treatmentCounts.order3oa;
                % Order-specific OA counters only tell the whole story at the shipped
                % MaxPolyOrder=3; this total stays correct if a rederivation pass lowers it.
                EEG.RELAXProcessingExtremeRejections.PopNOrderOA   = PopRepairInfo.treatmentCounts.order3oa ...
                    + PopRepairInfo.treatmentCounts.order2oa + PopRepairInfo.treatmentCounts.order1oa;
                EEG.RELAXProcessingExtremeRejections.PopNOrder3    = PopRepairInfo.treatmentCounts.order3;
                EEG.RELAXProcessingExtremeRejections.PopNOrder2    = PopRepairInfo.treatmentCounts.order2;
                EEG.RELAXProcessingExtremeRejections.PopNOrder1    = PopRepairInfo.treatmentCounts.order1;
                EEG.RELAXProcessingExtremeRejections.PopNMeanOnly  = PopRepairInfo.treatmentCounts.meanOnly;
                EEG.RELAXProcessingExtremeRejections.PopNFallback  = PopRepairInfo.treatmentCounts.fallback;
            end
        end
    else
        fprintf('Skipping electrode-pop repair (RELAX_cfg.Perform_pop_artifact_repair=0).\n');
    end

    %% Reject bridged electrode clusters detected via eBridge:
    % Bridged electrodes (electrolyte spreading between adjacent contacts) report a
    % near-identical, spatially smeared voltage; RELAX_bridged_channels flags them by an
    % electrical-distance criterion computed directly from the data (Alschuler et al.
    % 2014). This is a data-driven measurement, not metadata that could itself be wrong,
    % so nothing here is capped the way the impedance screen caps its rejections - every
    % flagged channel is removed. A bridge cannot be undone (both electrodes measure one
    % shorted voltage), so the removed channels are reconstructed later by spherical
    % interpolation from the surviving channels, the same treatment any other bad channel
    % receives (PREP- and impedance-rejected channels included). Note that interpolation
    % runs before GEDAI (RELAX_cfg.Run_GEDAI), same as for those other screens - see
    % stage 10 in this file's header comment / CLAUDE.md for why GEDAI needs bridged
    % channels gone in the first place.
    %
    % Runs after the pop stage (moved back here 2026-09-18; briefly ran before the pop
    % stage - see the pop-repair design doc addendum for that rationale and why it was
    % reverted): unlike the impedance screen, the bridge screen does not spend from the
    % shared channel-rejection budget the pop stage above also spends from, so its
    % position relative to that budget doesn't matter. It still runs before PREP.
    %
    % Absent = off: eBridge is an external EEGLAB plugin not bundled with this repo, and
    % this is new behaviour, so it must not switch itself on for a user who has not
    % installed it or opted in.
    if ~isfield(RELAX_cfg,'Perform_bridged_channel_rejection')
        RELAX_cfg.Perform_bridged_channel_rejection = 0;
    end
    if ~isfield(RELAX_cfg,'BridgeMaxProportionChannels')
        RELAX_cfg.BridgeMaxProportionChannels = 0.10; % QC flag only, see RELAX_bridged_channels
    end
    if ~isfield(RELAX_cfg,'BridgeMaxClusterSize')
        RELAX_cfg.BridgeMaxClusterSize = 4; % QC flag only, see RELAX_bridged_channels
    end
    EEG.RELAXProcessingExtremeRejections.BridgedChannelToReject={};
    if RELAX_cfg.Perform_bridged_channel_rejection == 1
        [BridgedChannelToReject, BridgeScreenInfo] = RELAX_bridged_channels(EEG, ...
            'RespectBoundaries', RELAX_cfg.RespectBoundaries, ...
            'MaxProportionChannels', RELAX_cfg.BridgeMaxProportionChannels, ...
            'MaxClusterSize', RELAX_cfg.BridgeMaxClusterSize, ...
            'FiltHigh', 20);
        % Wrapped in an outer cell - see the ImpedanceBasedChannelToReject comment above.
        EEG.RELAXProcessingExtremeRejections.BridgedChannelToReject = {BridgedChannelToReject};
        EEG.RELAX.BridgeScreen = BridgeScreenInfo;
        if ~BridgeScreenInfo.available
            warning('Bridge screen unavailable for %s (%s) - no channels rejected on bridging.', FileName, BridgeScreenInfo.status);
        elseif isempty(BridgedChannelToReject)
            fprintf('Bridge screen: no bridged channels detected (%d epochs).\n', BridgeScreenInfo.numEpochs);
        else
            if BridgeScreenInfo.excessiveBridging
                % RELAX_bridged_channels caps nothing itself (see its header) and
                % only flags this. Acting on the flag belongs here: with most of the
                % montage bridged (m396: 177 of 193 channels) there is no recording
                % left to clean, and deleting the channels and continuing produced a
                % file that crashed later at the blink stage anyway. Raised as an
                % exclusion so the file loop's catch records it as such.
                error('RELAX:excludeRecording', ...
                    ['excessive bridging: %d of %d channels (%.0f%%) bridged in %d cluster(s), ' ...
                     'above BridgeMaxProportionChannels (%.0f%%) - recording excluded'], ...
                    numel(BridgedChannelToReject), EEG.nbchan, 100*numel(BridgedChannelToReject)/EEG.nbchan, ...
                    numel(BridgeScreenInfo.clusters), 100*RELAX_cfg.BridgeMaxProportionChannels);
            end
            fprintf('Bridge screen: removing %d bridged channel(s) in %d cluster(s): %s\n', ...
                numel(BridgedChannelToReject), numel(BridgeScreenInfo.clusters), strjoin(BridgedChannelToReject', ', '));
            EEG=pop_select(EEG,'nochannel',BridgedChannelToReject);
            EEG = eeg_checkset( EEG );
        end
    else
        fprintf('Skipping bridged-channel rejection (RELAX_cfg.Perform_bridged_channel_rejection=0).\n');
    end

    %% Notch filter data:
    if strcmp(RELAX_cfg.NotchFilterType,'Butterworth')
        % Loop through list of freqs to reduce via notch filter
        if length(RELAX_cfg.LineNoiseFrequency) > 1
            for freq = RELAX_cfg.LineNoiseFrequency
                % Apply butterworth filter: 
                EEG = RELAX_filtbutter( EEG, freq-4, freq+4, 4, 'bandstop','acausal', RELAX_cfg.RespectBoundaries);
                if ~isempty(RELAX_cfg.electrodes_2_keep_but_not_clean{1})
                    Non_cleaned_electrodes = RELAX_filtbutter( Non_cleaned_electrodes, freq-4, freq+4, 4, 'bandstop','acausal', RELAX_cfg.RespectBoundaries);
                end
            end
        else
            % Apply butterworth filter 
            EEG = RELAX_filtbutter( EEG, RELAX_cfg.LineNoiseFrequency-4, RELAX_cfg.LineNoiseFrequency+4, 4, 'bandstop','acausal', RELAX_cfg.RespectBoundaries);
            if ~isempty(RELAX_cfg.electrodes_2_keep_but_not_clean{1})
                Non_cleaned_electrodes = RELAX_filtbutter( Non_cleaned_electrodes, RELAX_cfg.LineNoiseFrequency-4, RELAX_cfg.LineNoiseFrequency+4, 4, 'bandstop','acausal', RELAX_cfg.RespectBoundaries);
            end
        end
    end

    if strcmp(RELAX_cfg.NotchFilterType,'PMnotch')
        % v2.0.1 NWB adjusted to make PMnotch filtering available as an option, 8/8/2025
        % apply PMnotch filter via ERPLAB function (requires ERPLAB to be installed):
        EEG  = pop_basicfilter( EEG,  1:EEG.nbchan , 'Boundary', 'boundary', 'Cutoff',  RELAX_cfg.LineNoiseFrequency, 'Design', 'notch', 'Filter', 'PMnotch', 'Order',  180, 'RemoveDC','on');
        if ~isempty(RELAX_cfg.electrodes_2_keep_but_not_clean{1})
            Non_cleaned_electrodes  = pop_basicfilter( Non_cleaned_electrodes,  1:Non_cleaned_electrodes.nbchan , 'Boundary', 'boundary', 'Cutoff',  RELAX_cfg.LineNoiseFrequency, 'Design', 'notch', 'Filter', 'PMnotch', 'Order',  180, 'RemoveDC','on');
        end
    end

    %% Band Pass filter data: 
    if strcmp(RELAX_cfg.LowPassFilterBeforeMWF,'no') % updated implementation, avoiding low pass filtering prior to MWF reduces chances of rank deficiencies, increasing potential values for MWF delay period 
        if strcmp(RELAX_cfg.FilterType,'Butterworth')
            EEG = RELAX_filtbutter( EEG, RELAX_cfg.HighPassFilter, [], 4, 'highpass', RELAX_cfg.causal_or_acausal_filter, RELAX_cfg.RespectBoundaries);
            if ~isempty(RELAX_cfg.electrodes_2_keep_but_not_clean{1})
                % However, still low pass filter the auxilary electrodes, because they are not cleaned with MWF:
                if (RELAX_cfg.HighPassFilter_aux_elecs~=RELAX_cfg.HighPassFilter) || (RELAX_cfg.LowPassFilter_aux_elecs~=RELAX_cfg.LowPassFilter)
                    disp('Applying filter settings to auxilary electrodes:')
                end
                Non_cleaned_electrodes = RELAX_filtbutter( Non_cleaned_electrodes, RELAX_cfg.HighPassFilter_aux_elecs, RELAX_cfg.LowPassFilter_aux_elecs, 4, 'bandpass', RELAX_cfg.causal_or_acausal_filter, RELAX_cfg.RespectBoundaries);
            end
        end
        if strcmp(RELAX_cfg.FilterType,'pop_eegfiltnew')
            EEG = pop_eegfiltnew(EEG,RELAX_cfg.HighPassFilter,[]);
            if ~isempty(RELAX_cfg.electrodes_2_keep_but_not_clean{1})
                if (RELAX_cfg.HighPassFilter_aux_elecs~=RELAX_cfg.HighPassFilter) || (RELAX_cfg.LowPassFilter_aux_elecs~=RELAX_cfg.LowPassFilter)
                    disp('Applying filter settings to auxilary electrodes:')
                end
                % However, still low pass filter the auxilary electrodes, because they are not cleaned with MWF:
                Non_cleaned_electrodes = pop_eegfiltnew(Non_cleaned_electrodes,RELAX_cfg.HighPassFilter_aux_elecs,RELAX_cfg.LowPassFilter_aux_elecs);
            end
        end
    end
    
    if strcmp(RELAX_cfg.LowPassFilterBeforeMWF,'yes') % original implementation, not recommended before MWF unless downsampling, as increases chances of rank deficiencies
        if strcmp(RELAX_cfg.FilterType,'Butterworth')
            EEG = RELAX_filtbutter( EEG, RELAX_cfg.HighPassFilter, RELAX_cfg.LowPassFilter, 4, 'bandpass', RELAX_cfg.causal_or_acausal_filter, RELAX_cfg.RespectBoundaries);
            if ~isempty(RELAX_cfg.electrodes_2_keep_but_not_clean{1})
                if (RELAX_cfg.HighPassFilter_aux_elecs~=RELAX_cfg.HighPassFilter) || (RELAX_cfg.LowPassFilter_aux_elecs~=RELAX_cfg.LowPassFilter)
                    disp('Applying filter settings to auxilary electrodes:')
                end
                Non_cleaned_electrodes = RELAX_filtbutter( Non_cleaned_electrodes, RELAX_cfg.HighPassFilter_aux_elecs, RELAX_cfg.LowPassFilter_aux_elecs, 4, 'bandpass', RELAX_cfg.causal_or_acausal_filter, RELAX_cfg.RespectBoundaries);
            end
        end
        if strcmp(RELAX_cfg.FilterType,'pop_eegfiltnew')
            EEG = pop_eegfiltnew(EEG,RELAX_cfg.HighPassFilter,RELAX_cfg.LowPassFilter_aux_elecs);
            if ~isempty(RELAX_cfg.electrodes_2_keep_but_not_clean{1})
                if (RELAX_cfg.HighPassFilter_aux_elecs~=RELAX_cfg.HighPassFilter) || (RELAX_cfg.LowPassFilter_aux_elecs~=RELAX_cfg.LowPassFilter)
                    disp('Applying filter settings to auxilary electrodes:')
                end
                Non_cleaned_electrodes = pop_eegfiltnew(Non_cleaned_electrodes,RELAX_cfg.HighPassFilter_aux_elecs,RELAX_cfg.LowPassFilter_aux_elecs);
            end
        end
    end

    if strcmp(RELAX_cfg.DownSample,'yes')
        EEG = pop_resample(EEG,RELAX_cfg.DownSample_to_X_Hz); % downsample data (if applied, should always be applied after low pass filtering)
        RELAX_cfg.ms_per_sample=(1000/EEG.srate);
        if ~isempty(RELAX_cfg.electrodes_2_keep_but_not_clean{1})
            Non_cleaned_electrodes = pop_resample(Non_cleaned_electrodes,RELAX_cfg.DownSample_to_X_Hz);
        end
    end

    if RELAX_cfg.ms_per_sample<0.7
        warning('The sampling rate for this file is quite high. Depending on your processing power, RELAX may run slowly or even stall, especially if applying MWF cleaning');
        warning('RELAX was validated using 1000Hz sampling rates, which is still a high sample rate for most analyses. You could downsample your data by setting the relevant options in RELAX');
    end

    if strcmp(RELAX_cfg.NotchFilterType,'ZaplinePlus')
        zaplineConfig = struct('plotResults',0, ...
                               'detailedFreqBoundsUpper', [-0.2 0.2], ...
                               'detailedFreqBoundsLower', [-0.8 0.4], ...
                               'searchIndividualNoise', false);
        [EEG ] = clean_data_with_zapline_plus_eeglab_wrapper(EEG, zaplineConfig); % requires the zapline plus plugin. Best applied after downsampling to 250Hz or 500Hz.
        if ~isempty(RELAX_cfg.electrodes_2_keep_but_not_clean{1})
            % zapline plus won't fix the auxilary electrodes, so Butterworth filtering them instead:
            disp('Applying notch filter to auxilary electrodes as zapline will not work on them:')
            Non_cleaned_electrodes = RELAX_filtbutter( Non_cleaned_electrodes, RELAX_cfg.LineNoiseFrequency-3, RELAX_cfg.LineNoiseFrequency+3, 4, 'bandstop','acausal', RELAX_cfg.RespectBoundaries);
        end
    end

    %% Clean flat channels and bad channels showing improbable data:
    % PREP pipeline: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4471356/
    % NOTE: findNoisyChannels states its own precondition "No segments of the EEG data
    % have been removed" and tiles its correlation windows flatly from sample 1, so on
    % data with boundary events a large fraction of windows straddle a discontinuity and
    % channel rejection becomes close to arbitrary. It is therefore skipped by default
    % here; the epoching stage immediately afterwards (RELAX_excluding_channels_and_epoching,
    % or its _boundary variant when RELAX_cfg.RespectBoundaries=1 - see CLAUDE.md's
    % "Boundary/discontinuity handling") performs its own channel rejection on epochs that
    % never straddle a discontinuity. Set Perform_PREP_bad_channel_rejection=1 to restore
    % the original behaviour on data with no discontinuities.
    if ~isfield(RELAX_cfg,'Perform_PREP_bad_channel_rejection')
        RELAX_cfg.Perform_PREP_bad_channel_rejection = 1; % absent = original behaviour
    end
    EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject={};
    noisyOut = struct('noisyChannels', struct('all', []));
    if RELAX_cfg.Perform_PREP_bad_channel_rejection == 1
        addpath(genpath(PrepFolderLocation{1,1})); %  1.1.4: fix error where PREP seems to be removed from the path after an EEGLAB update
        noisyOut = findNoisyChannels(EEG);
        for x=1:size(noisyOut.noisyChannels.all,2) % loop through output of PREP's findNoisyChannels and take a record of noisy electrodes for deletion:
            PREPBasedChannelToReject{x}=EEG.chanlocs(noisyOut.noisyChannels.all(x)).labels;
        end
        if size(noisyOut.noisyChannels.all,2) > 0
            % Wrapped in an outer cell - see the PopBasedChannelToReject comment above.
            EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject = {PREPBasedChannelToReject'};
        end
    else
        fprintf('Skipping PREP findNoisyChannels (RELAX_cfg.Perform_PREP_bad_channel_rejection=0).\n');
    end

    % Print channel names to be removed by PREP Noisy Channels
    if ~isempty(noisyOut.noisyChannels.all)
        noisy_labels = {EEG.chanlocs(noisyOut.noisyChannels.all).labels};
        fprintf('Removing noisy channels: %s\n', strjoin(noisy_labels, ', '));
    else
        fprintf('No noisy channels detected by PREP\n');
    end
    
    EEG=pop_select(EEG,'nochannel',noisyOut.noisyChannels.all); % delete noisy electrodes detected by PREP

    continuousEEG=EEG;

    if RELAX_cfg.RespectBoundaries==1
        % Boundary-aware variant: epochs are tiled within each inter-boundary segment
        % (via RELAX_get_segments) instead of on a fixed grid across the whole
        % recording, so data adjacent to a 'boundary' event is still evaluated rather
        % than left un-epoched (NaN) and later swept into ExtremelyBadPeriodsForDeletion
        % regardless of whether it's actually extreme.
        [continuousEEG, epochedEEG] = RELAX_excluding_channels_and_epoching_boundary(continuousEEG, RELAX_cfg); % Epoch data, detect extremely bad data, delete channels if over the set threshold for proportion of data affected by extreme outlier for each electrode
    else
        [continuousEEG, epochedEEG] = RELAX_excluding_channels_and_epoching(continuousEEG, RELAX_cfg); % Epoch data, detect extremely bad data, delete channels if over the set threshold for proportion of data affected by extreme outlier for each electrode
    end
    [continuousEEG, epochedEEG] = RELAX_excluding_extreme_values(continuousEEG, epochedEEG, RELAX_cfg); % Mark extreme periods for exclusion from MWF cleaning, and deletion before wICA cleaning

    % Use the continuous data to detect eye blinks and mark
    % these in the EEG.event as well as in the mask. The output is
    % continuous data but includes all the previous extreme period 
    % markings from the epoched data.
    if RELAX_cfg.ProbabilityDataHasNoBlinks<2
        [continuousEEG, epochedEEG] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, RELAX_cfg); % use an IQR threshold method to detect and mark blinks
        if continuousEEG.RELAX.IQRmethodDetectedBlinks(1,1)==0 % If a participants doesn't show any blinks, make a note
            NoBlinksDetected{FileNumber,1}=FileName; 
            warning('No blinks were detected - if blinks are expected then you should visually inspect the file');
        end
        if RELAX_cfg.computerawmetrics==1
        [continuousEEG, epochedEEG] = RELAX_metrics_blinks(continuousEEG, epochedEEG); % record blink amplitude ratio from raw data for comparison.
        end
    end

    % Record the channel count retained through this point, so the effect of skipping
    % PREP bad-channel rejection (RELAX_cfg.Perform_PREP_bad_channel_rejection=0, see
    % above) is visible in the aggregated per-participant table, not just in stdout.
    epochedEEG.RELAXProcessingExtremeRejections.RetainedChannelCount = epochedEEG.nbchan;

    % Record extreme artifact rejection details for all participants in single table:
    ExtremeRejectionsRow = struct2table(epochedEEG.RELAXProcessingExtremeRejections,'AsArray',true);
    % Via RELAX_append_table_row, because a table checkpointed by an earlier version of
    % this wrapper can be missing columns this version writes
    % (ImpedanceBasedChannelToReject was added after some runs had already been
    % checkpointed) and plain row assignment errors on a differing variable set, which
    % would break resuming a partially processed participant:
    if ~exist('RELAXProcessingExtremeRejectionsAllParticipants','var')
        RELAXProcessingExtremeRejectionsAllParticipants = [];
    end
    RELAXProcessingExtremeRejectionsAllParticipants = RELAX_append_table_row( ...
        RELAXProcessingExtremeRejectionsAllParticipants, FileNumber, ExtremeRejectionsRow);

    rawEEG=continuousEEG; % Take a copy of the not yet cleaned data for calculation of all cleaning SER and ARR at the end
    
    %% Mark artifacts for calculating SER and ARR, regardless of whether MWF is performed (RELAX v1.1.3 update):
    % NOTE: this block exists only to build NoiseMaskFullLengthR1, whose sole consumer
    % is RELAX_metrics_final_SER_and_ARR - and that call is commented out further down
    % this file. It is off by default because RELAX_muscle, RELAX_horizontaleye,
    % RELAX_drift and RELAX_pad_brief_mask_periods are all boundary-blind and would need
    % boundary handling before this block can be re-enabled on discontinuous data.
    if ~isfield(RELAX_cfg,'ComputeSER_ARR'); RELAX_cfg.ComputeSER_ARR = 0; end
    if RELAX_cfg.ComputeSER_ARR==1 && RELAX_cfg.computecleanedmetrics==1 && (RELAX_cfg.Do_MWF_Once==0 || RELAX_cfg.Do_MWF_Twice==0 || RELAX_cfg.Do_MWF_Thrice==0)
        [Marking_artifacts_for_SER_ARR, ~] = RELAX_muscle(continuousEEG, epochedEEG, RELAX_cfg); 
        Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(Marking_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength==1)=1; 
        [Marking_artifacts_for_SER_ARR] = RELAX_horizontaleye(continuousEEG, RELAX_cfg); 
        Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(Marking_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength==1)=1; 
        [Marking_artifacts_for_SER_ARR, ~] = RELAX_drift(continuousEEG, epochedEEG, RELAX_cfg); % Use epoched data to add periods showing excessive drift to the mask 
        Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(Marking_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength==1)=1; 
        Marking_all_artifacts_for_SER_ARR.RELAX.NaNsForExtremeOutlierPeriods=continuousEEG.RELAX.NaNsForExtremeOutlierPeriods; 
        [Marking_all_artifacts_for_SER_ARR] = RELAX_pad_brief_mask_periods (Marking_all_artifacts_for_SER_ARR, RELAX_cfg, 'notblinks'); % If period has been marked as shorter than RELAX_cfg.MinimumArtifactDuration, then pad it out. 
        if isfield(continuousEEG.RELAX,'eyeblinkmask')
                Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(continuousEEG.RELAX.eyeblinkmask==1)=1; 
                [Marking_all_artifacts_for_SER_ARR] = RELAX_pad_brief_mask_periods (Marking_all_artifacts_for_SER_ARR, RELAX_cfg, 'blinks'); 
        end
        continuousEEG.RELAX.NoiseMaskFullLengthR1=Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength; 
        rawEEG.RELAX.NoiseMaskFullLengthR1=Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength; 
    end
    
    if RELAX_cfg.saveextremesrejected==1
        if ~exist([RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep 'Extremes_Rejected'], 'dir')
            mkdir([RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep 'Extremes_Rejected'])
        end
        SaveSetExtremes_Rejected =[RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep 'Extremes_Rejected', filesep FileName '_Extremes_Rejected.set'];    
        EEG = pop_saveset( rawEEG, SaveSetExtremes_Rejected ); % If desired, save data here with bad channels deleted, filtering applied, extreme outlying data periods marked
    end

    %% THIS SECTION CONTAINS FUNCTIONS WHICH MARK AND CLEAN MUSCLE ARTIFACTS
    % Any one of these functions can be commented out to ignore those artifacts
    % when creating the mask    
    if RELAX_cfg.Do_MWF_Once==1

        % Use epoched data and FFT to detect slope of log frequency log
        % power, add periods exceeding muscle threshold to mask:
        [continuousEEG, epochedEEG] = RELAX_muscle(continuousEEG, epochedEEG, RELAX_cfg);  
        if RELAX_cfg.computerawmetrics==1
            [continuousEEG, epochedEEG] = RELAX_metrics_muscle(continuousEEG, epochedEEG, RELAX_cfg); % record muscle contamination metrics from raw data for comparison.
        end

        EEG=continuousEEG; % Return continuousEEG to the "EEG" variable for MWF processing

        % If including eye blink cleaning in first round MWF, then insert
        % eye blink mask into noise mask:
        if RELAX_cfg.MWFRoundToCleanBlinks==1
            EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
            EEG.RELAX.eyeblinkmask(isnan(EEG.RELAXProcessing.Details.NaNsForNonEvents))=NaN;
            EEG.RELAXProcessing.ProportionMarkedBlinks=mean(EEG.RELAX.eyeblinkmask,'omitnan');
        end

        % The following pads very brief lengths of mask periods
        % in the template (without doing this, very short periods can
        % lead to rank deficiency), and excludes extreme artifacts from the
        % cleaning template (so the MWF cleaning step just ignores extreme
        % artifacts in it's template - doesn't include them in either the
        % clean or artifact mask, but does apply cleaning to them).
        [EEG] = RELAX_pad_brief_mask_periods (EEG, RELAX_cfg, 'notblinks'); % If period has been marked as shorter than RELAX_cfg.MinimumArtifactDuration, then pad it out.
        
        EEG.RELAX.NoiseMaskFullLengthR1=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
        EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal=mean(EEG.RELAXProcessing.Details.NoiseMaskFullLength,'omitnan');
        EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR1=EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal; 
  
        %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
        RELAX_cfg.MWFDelayPeriod=RELAX_cfg.MWFDelayPeriod_for_muscle_artifacts; 
        RELAX_cfg.MWF_delay_spacing=RELAX_cfg.MWF_delay_spacing_for_muscle_artifacts;
        [EEG] = RELAX_perform_MWF_cleaning (EEG, RELAX_cfg);          

        EEG.RELAXProcessingRoundOne=EEG.RELAXProcessing; % Record MWF cleaning details from round 1 in EEG file          
        RELAXProcessingRoundOne=EEG.RELAXProcessingRoundOne; % Record MWF cleaning details from round 1 into file for all participants
        
        if isfield(RELAXProcessingRoundOne,'Details')
            RELAXProcessingRoundOne=rmfield(RELAXProcessingRoundOne,'Details');
        end
        if RELAX_cfg.KeepAllInfo==0
            if isfield(EEG.RELAXProcessingRoundOne,'Details')
                EEG.RELAXProcessingRoundOne=rmfield(EEG.RELAXProcessingRoundOne,'Details');
            end
        end
        
        % Record processing statistics for all participants in single table:
        RELAXProcessingMWFStepOneAllParticipants(FileNumber,:) = struct2table(RELAXProcessingRoundOne,'AsArray',true);
        EEG = rmfield(EEG,'RELAXProcessing');
        % Save round 1 MWF pre-processing:
        if RELAX_cfg.saveround1==1
            if ~exist([RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep '1xMWF'], 'dir')
                mkdir([RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep '1xMWF'])
            end
            SaveSetMWF1 =[RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep '1xMWF', filesep FileName '_MWF1.set'];    
            EEG = pop_saveset( EEG, SaveSetMWF1 ); 
        end
    end

    %% PERFORM A SECOND ROUND OF MWF. THIS IS HELPFUL IF THE FIRST ROUND DOESN'T SUFFICIENTLY CLEAN ARTIFACTS. 

    % This has been suggested to be useful by Somers et al (2018)
    % (particularly when used in a cascading fashion). 

    % However, I can see risks. If artifact masks fall on task relevant
    % activity in both rounds of the MWF, it may be that the task relevant data
    % is just cleaned right out of the signal.
    
    if RELAX_cfg.Do_MWF_Twice==1

        % v2.0.1 NWB added the following so that if muscle artifacts aren't
        % cleaned by MWF, the second step doesn't go back to the raw data:
        if RELAX_cfg.Do_MWF_Once==0
            EEG=continuousEEG;
        end

        EEG.RELAXProcessing.aFileName=cellstr(FileName);
        EEG.RELAXProcessing.ProportionMarkedBlinks=0;
        
        % If blinks weren't initially detected because they were 
        % disguised by the the muscle artifact, detect them here
        % (this happens in <1/200 cases, but is a good back up).
        if RELAX_cfg.ProbabilityDataHasNoBlinks==0
            if EEG.RELAX.IQRmethodDetectedBlinks(1,1)==0
                continuousEEG=EEG;
                [continuousEEG, epochedEEG] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, RELAX_cfg);
                EEG=continuousEEG;
            end
        end
        
        % If including eye blink cleaning in second round MWF, then insert
        % eye blink mask into noise mask:
        if isfield(EEG.RELAX, 'eyeblinkmask')
            if RELAX_cfg.MWFRoundToCleanBlinks==2
                EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
                EEG.RELAX.eyeblinkmask(isnan(EEG.RELAX.NaNsForExtremeOutlierPeriods))=NaN;
                EEG.RELAXProcessing.ProportionMarkedBlinks=mean(EEG.RELAX.eyeblinkmask,'omitnan');
            end
        end
  
        % The following pads very brief lengths of mask periods
        % in the template (without doing this, very short periods can
        % lead to rank deficiency), and excludes extreme artifacts from the
        % cleaning template (so the MWF cleaning step just ignores extreme
        % artifacts in it's template - doesn't include them in either the
        % clean or artifact mask, but does apply cleaning to them).
        [EEG] = RELAX_pad_brief_mask_periods (EEG, RELAX_cfg, 'blinks');
        
        EEG.RELAX.NoiseMaskFullLengthR2=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
        EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal=mean(EEG.RELAXProcessing.Details.NoiseMaskFullLength,'omitnan');
        EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR2=EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal; 

        %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
        RELAX_cfg.MWFDelayPeriod=RELAX_cfg.MWFDelayPeriod_for_eye_movements; 
        RELAX_cfg.MWF_delay_spacing=RELAX_cfg.MWF_delay_spacing_for_eye_movements; % set how sparsely the delay stacking is spread
        [EEG] = RELAX_perform_MWF_cleaning (EEG, RELAX_cfg);  
 
        EEG.RELAXProcessingRoundTwo=EEG.RELAXProcessing; % Record MWF cleaning details from round 2 in EEG file
        RELAXProcessingRoundTwo=EEG.RELAXProcessingRoundTwo; % Record MWF cleaning details from round 2 into file for all participants
        if isfield(RELAXProcessingRoundTwo,'Details')
            RELAXProcessingRoundTwo=rmfield(RELAXProcessingRoundTwo,'Details');
        end
        if RELAX_cfg.KeepAllInfo==0
            if isfield(EEG.RELAXProcessingRoundTwo,'Details')
                EEG.RELAXProcessingRoundTwo=rmfield(EEG.RELAXProcessingRoundTwo,'Details');
            end
        end
        % Record processing statistics for all participants in single table:
        RELAXProcessingMWFStepTwoAllParticipants(FileNumber,:) = struct2table(RELAXProcessingRoundTwo,'AsArray',true);
        EEG = rmfield(EEG,'RELAXProcessing');
        % Save round 2 MWF pre-processing:
        if RELAX_cfg.saveround2==1
            if ~exist([RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep '2xMWF'], 'dir')
                mkdir([RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep '2xMWF'])
            end
            SaveSetMWF2 =[RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep '2xMWF', filesep FileName '_MWF2.set'];    
            EEG = pop_saveset( EEG, SaveSetMWF2 ); 
        end     
    end
    
    %% PERFORM A THIRD ROUND OF MWF.    
    if RELAX_cfg.Do_MWF_Thrice==1

        % v2.0.1 NWB added the following so that if muscle artifacts & blinks 
        % aren't cleaned by MWF, the 3rd step doesn't go back to the raw data:
        if RELAX_cfg.Do_MWF_Once==0 && RELAX_cfg.Do_MWF_Twice==0
            EEG=continuousEEG;
        end

        EEG.RELAXProcessing.aFileName=cellstr(FileName);
        EEG.RELAXProcessing.ProportionMarkedBlinks=0;
        % If less than 5% of data was masked as eye blink cleaning in second round MWF, then insert
        % eye blink mask into noise mask in round 3:
        if isfield(EEG.RELAX,'ProportionMarkedInMWFArtifactMaskTotalR2') % NWB added to make sure function doesn't bug when trying to check this variable if it doesn't exist
            if EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR2<0.05
                if isfield(EEG.RELAX, 'eyeblinkmask')
                    EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
                    EEG.RELAX.eyeblinkmask(isnan(EEG.RELAX.NaNsForExtremeOutlierPeriods))=NaN;
                    EEG.RELAXProcessing.ProportionMarkedBlinks=mean(EEG.RELAX.eyeblinkmask,'omitnan');
                end
            end
        end

        % Epoch the data into 1 second epochs with a 500ms overlap. Outputs
        % both the ContinuousEEG (which has been filtered above by this
        % point) and the epoched data as EEG.
        [continuousEEG, epochedEEG] = RELAX_epoching(EEG, RELAX_cfg);
        
        %% THIS SECTION CONTAINS FUNCTIONS WHICH MARK ARTIFACTS

        [continuousEEG, epochedEEG] = RELAX_drift(continuousEEG, epochedEEG, RELAX_cfg); % Use epoched data to add periods showing excessive drift to the mask
        
        % Use the filtered continuous data to detect horizontal eye
        % movements and mark these in the EEG.event as well as in the mask.
        % You may want to simply reject horizontal eye movements at a later
        % stage if your task requires participants to look straight ahead
        % for the entire task. Alternatively, if your task requires
        % participants to complete horizontal eye movements time locked to
        % a stimuli, this section will mark every event with these
        % horizontal eye movements as an artifact, and should not be
        % implemented.
        
        % The output is continuous data:
        [continuousEEG] = RELAX_horizontaleye(continuousEEG, RELAX_cfg);
    
        %% Return to the "EEG" variable for MWF processing:
        EEG=continuousEEG;
        
        % If including eye blink cleaning in third round MWF, then insert
        % eye blink mask into noise mask:
        if RELAX_cfg.MWFRoundToCleanBlinks==3
            EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
            EEG.RELAX.eyeblinkmask(isnan(EEG.RELAXProcessing.Details.NaNsForNonEvents))=NaN;
            EEG.RELAXProcessing.ProportionMarkedBlinks=mean(EEG.RELAX.eyeblinkmask,'omitnan');
        end
 
        % The following pads very brief lengths of mask periods
        % in the template (without doing this, very short periods can
        % lead to rank deficiency), and excludes extreme artifacts from the
        % cleaning template (so the MWF cleaning step just ignores extreme
        % artifacts in it's template - doesn't include them in either the
        % clean or artifact mask, but does apply cleaning to them).
        [EEG] = RELAX_pad_brief_mask_periods (EEG, RELAX_cfg, 'notblinks');
        
        EEG.RELAX.NoiseMaskFullLengthR3=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
        EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal=mean(EEG.RELAXProcessing.Details.NoiseMaskFullLength,'omitnan');
        EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR3=EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal; 

        %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
        RELAX_cfg.MWFDelayPeriod=RELAX_cfg.MWFDelayPeriod_for_eye_movements; 
        RELAX_cfg.MWF_delay_spacing=RELAX_cfg.MWF_delay_spacing_for_eye_movements; % set how sparsely the delay stacking is spread
        [EEG] = RELAX_perform_MWF_cleaning (EEG, RELAX_cfg);               
        
        if isfield(EEG.RELAX, 'eyeblinkmask') % if eyeblinkmask has been created, do the following (thanks to Jane Tan for the suggested bug fix when eyeblinkmask is not created)
            EEG.RELAX=rmfield(EEG.RELAX,'eyeblinkmask'); % remove variables that are no longer necessary
        end
        
        EEG.RELAXProcessingRoundThree=EEG.RELAXProcessing; % Record MWF cleaning details from round 3 in EEG file
        RELAXProcessingRoundThree=EEG.RELAXProcessing; % Record MWF cleaning details from round 3 into file for all participants
        
        if isfield(RELAXProcessingRoundThree,'Details')
            RELAXProcessingRoundThree=rmfield(RELAXProcessingRoundThree,'Details');
        end
        if RELAX_cfg.KeepAllInfo==0
            if isfield(EEG.RELAXProcessingRoundThree,'Details')
                EEG.RELAXProcessingRoundThree=rmfield(EEG.RELAXProcessingRoundThree,'Details');
            end
        end
        % Record processing statistics for all participants in single table:
        RELAXProcessingMWFStepThreeAllParticipants(FileNumber,:) = struct2table(RELAXProcessingRoundThree,'AsArray',true);
        EEG = rmfield(EEG,'RELAXProcessing');

        if RELAX_cfg.saveround3==1
            if ~exist([RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep '3xMWF'], 'dir')
                mkdir([RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep '3xMWF'])
            end
            SaveSetMWF3 =[RELAX_cfg.foldername,filesep 'RELAXProcessed' filesep '3xMWF', filesep FileName '_MWF3.set'];    
            EEG = pop_saveset( EEG, SaveSetMWF3 ); 
        end         
    end
    
    %% Perform robust average re-referencing of the data, reject periods marked as extreme outliers    
    if RELAX_cfg.Do_MWF_Once==0 && RELAX_cfg.Do_MWF_Twice==0 && RELAX_cfg.Do_MWF_Thrice==0 % v2.0.1 NWB adjusted to only return to continuous EEG if no MWF was applied (rather than just that the first MWF wasn't applied)
        EEG=continuousEEG;
    end

    % v2.0.1 NWB 8/8/2025 - moved the low pass filtering after the MWF
    % steps to before the extreme bad period rejection steps to ensure
    % filtering is not adversely affected by boundary elements.
    if strcmp(RELAX_cfg.LowPassFilterBeforeMWF,'no') % if low pass filtering wasn't applied before MWF cleaning (recommended) apply it here
        if strcmp(RELAX_cfg.FilterType,'Butterworth')
            EEG = RELAX_filtbutter( EEG, [], RELAX_cfg.LowPassFilter, 4, 'lowpass', RELAX_cfg.causal_or_acausal_filter, RELAX_cfg.RespectBoundaries);
        end
        if strcmp(RELAX_cfg.FilterType,'pop_eegfiltnew')
            EEG = pop_eegfiltnew(EEG,[],RELAX_cfg.LowPassFilter);
        end
    end
          
    % Reject periods that were marked as NaNs in the MWF masks because they 
    % showed extreme shift within the epoch or extremely improbable data:
    EEG = eeg_eegrej( EEG, EEG.RELAX.ExtremelyBadPeriodsForDeletion);

    if ~isempty(RELAX_cfg.electrodes_2_keep_but_not_clean{1})
        % reject the same periods from the auxilary electrodes that weren't included in cleaning:
        Non_cleaned_electrodes = eeg_eegrej( Non_cleaned_electrodes, EEG.RELAX.ExtremelyBadPeriodsForDeletion);
        % store the auxilary electrodes in the EEG file:
        EEG.Non_cleaned_electrodes=Non_cleaned_electrodes;
    end
 
    [EEG] = RELAX_average_rereference(EEG);
    EEG = eeg_checkset( EEG );  

    %% Perform targeted wICA on ICLabel identified artifacts that remain:
    if RELAX_cfg.Perform_targeted_wICA==1
        % The following cleans eye movements and muscle artifacts in the
        % independent component space by a combination of wavelet enhanced
        % ICA cleaning and targeting to restrict the cleaning to only 
        % artifact periods for eye movement components, and high pass 
        % filtering muscle components at 15Hz:
        EEG.RELAXProcessing_wICA.aFileName=cellstr(FileName);
        [EEG] = RELAX_targeted_wICA(EEG,RELAX_cfg);
        % setting 'RELAX_cfg.Report_all_wICA_info' to 1 will report proportion of ICs categorized as each category, and variance explained by ICs from each category (function is ~20s slower if this is implemented)
        EEG = eeg_checkset( EEG );
        RELAXProcessing_wICA=EEG.RELAXProcessing_wICA;
        % Record processing statistics for all participants in single table:
        RELAXProcessing_wICA_AllParticipants(FileNumber,:) = struct2table(RELAXProcessing_wICA,'AsArray',true);
    end
    
    
    %% Perform wICA on ICLabel identified artifacts that remain:
    if RELAX_cfg.Perform_wICA_on_ICLabel==1
        % The following performs wICA, implemented on only the components
        % marked as artifact by ICLabel.
        EEG.RELAXProcessing_wICA.aFileName=cellstr(FileName);
        [EEG,~, ~, ~, ~] = RELAX_wICA_on_ICLabel_artifacts(EEG,RELAX_cfg.ICA_method, 1, 0, EEG.srate, 5,'coif5',RELAX_cfg.Report_all_ICA_info,RELAX_cfg.ICLabel_thresholds,RELAX_cfg.Clean_other_comps); 
        % setting 'RELAX_cfg.Report_all_wICA_info' to 1 will report proportion of ICs categorized as each category, and variance explained by ICs from each category (function is ~20s slower if this is implemented)
        EEG = eeg_checkset( EEG );
        RELAXProcessing_wICA=EEG.RELAXProcessing_wICA;
        % Record processing statistics for all participants in single table:
        RELAXProcessing_wICA_AllParticipants(FileNumber,:) = struct2table(RELAXProcessing_wICA,'AsArray',true);
    end
    
    %% Perform ICA subtract on ICLabel identified artifacts that remain:
    if RELAX_cfg.Perform_ICA_subtract==1
        % The following performs ICA sutraction, implemented on only the components
        % marked as artifact by ICLabel.
        EEG.RELAXProcessing_ICA.aFileName=cellstr(FileName);
        EEG = RELAX_ICA_subtract(EEG,RELAX_cfg);
        EEG = eeg_checkset( EEG );
        RELAXProcessing_ICA=EEG.RELAXProcessing_ICA;
        % Record processing statistics for all participants in single table:
        RELAXProcessing_ICA_AllParticipants(FileNumber,:) = struct2table(RELAXProcessing_ICA,'AsArray',true);
    end
    
    EEG.RELAX.Data_has_been_cleaned=1;
   
    %% Save ICA topoplots
    if exist('save_ICA_topo.m', 'file') && ~isempty(EEG.icaweights)
        save_ICA_topo(EEG, RELAX_cfg);
    end
    
    %% Run GEDAI pipeline - IMK added
    % Run the Generalized Eigenvalue Deartifacting Instrument (GEDAI)
    % pipeline for additional denoising
    if RELAX_cfg.Run_GEDAI==1
        % Interpolate rejected electrodes prior to running GEDAI
        % This ensures dimension matching with full-rank participant refCOV
        EEG = pop_interp(EEG, EEG.allchan, 'spherical');
        
        % Load custom refCOV for paticipant file
        %refCOV_ppt = loadGEDAI_ppt_refCOV(RELAX_cfg);

        % Run GEDAI with the custom refCOV and otherwise default params
        % Retain the cleaned continous EEG and artifact outputs for later use
        [EEG, EEGartifacts, sensai_fig] = GEDAI(EEG, 'auto', 12, 0.5, 'custom', true, false, 0.80, 'eeg');
        
        % Save figure
        sensai_plots_dir = fullfile(RELAX_cfg.myPath, 'RELAXProcessed', 'SENSAI_plots');
        if ~exist(sensai_plots_dir, 'dir')
            mkdir(sensai_plots_dir);
        end

        if ~isempty(sensai_fig)
            saveas(sensai_fig, fullfile(sensai_plots_dir, sprintf('SENSAI_vis_%s.png', FileName)));
            fprintf('Saving SENSAI plot fig to: %s\n', sensai_plots_dir);
            close(sensai_fig)
        else
            warning('No SENSAI figure generated for %s\n', FileName);
        end
   end
      
    %% COMPUTE CLEANED METRICS:
    if RELAX_cfg.computecleanedmetrics==1    
        [continuousEEG, epochedEEG] = RELAX_epoching(EEG, RELAX_cfg);
        [continuousEEG, ~] = RELAX_metrics_blinks(continuousEEG, epochedEEG);
        [continuousEEG, ~] = RELAX_metrics_muscle(continuousEEG, epochedEEG, RELAX_cfg);
        
        % Metrics are computed only on "good" chans. Interpolated chans are discarded
        % in this step and must be reinterpolated later. IMK modified function.
        %[continuousEEG] = RELAX_metrics_final_SER_and_ARR(rawEEG, continuousEEG); % this is only a good metric for testing only the cleaning of artifacts marked for cleaning by MWF, see notes in function.

        EEG=continuousEEG;
        EEG = rmfield(EEG,'RELAXProcessing');

        if isfield(EEG,'RELAX_Metrics')
            if isfield(EEG.RELAX_Metrics, 'Cleaned')
                if isfield(EEG.RELAX_Metrics.Cleaned,'BlinkAmplitudeRatio')
                    CleanedMetrics.BlinkAmplitudeRatio(1:size(EEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio,1),FileNumber)=EEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio;
                    CleanedMetrics.BlinkAmplitudeRatio(CleanedMetrics.BlinkAmplitudeRatio==0)=NaN;
                end
                if isfield(EEG.RELAX_Metrics.Cleaned,'MeanMuscleStrengthFromOnlySuperThresholdValues')
                    CleanedMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues(FileNumber)=EEG.RELAX_Metrics.Cleaned.MeanMuscleStrengthFromOnlySuperThresholdValues; 
                    CleanedMetrics.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel(FileNumber)=EEG.RELAX_Metrics.Cleaned.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel;
                end
                % if isfield(EEG.RELAX_Metrics.Cleaned,'All_SER')
                %     CleanedMetrics.All_SER(FileNumber)=EEG.RELAX_Metrics.Cleaned.All_SER;
                %     CleanedMetrics.All_ARR(FileNumber)=EEG.RELAX_Metrics.Cleaned.All_ARR;
                %end
                if isfield(EEG.etc,'GEDAI')  % IMK added - store GEDAI metrics in RELAX.cfg for export

                    g = EEG.etc.GEDAI;

                    % Core metrics
                    CleanedMetrics.GEDAI_SENSAI_score(FileNumber)           = g.SENSAI_score;
                    CleanedMetrics.GEDAI_mean_ENOVA(FileNumber)             = g.mean_ENOVA;
                    CleanedMetrics.GEDAI_epochs_rejected(FileNumber)        = g.epochs_rejected;
                    CleanedMetrics.GEDAI_total_epochs(FileNumber)           = g.total_epochs;
                    CleanedMetrics.GEDAI_percent_rejected(FileNumber)       = g.percentage_rejected;

                    if isfield(g, 'visualization_metrics')
                        CleanedMetrics.SENSAI_ssi_before_mean(FileNumber)       = g.visualization_metrics.ssi_before_mean;
                        CleanedMetrics.SENSAI_SSSI(FileNumber)                  = g.visualization_metrics.ssi_after_mean;
                        CleanedMetrics.SENSAI_NSSI(FileNumber)                  = g.visualization_metrics.ssi_artifacts_mean;
                        CleanedMetrics.SENSAI_lda_accuracy(FileNumber)          = g.visualization_metrics.lda_accuracy;
                        CleanedMetrics.SENSAI_signal_silhouette(FileNumber)     = g.visualization_metrics.signal_silhouette;
                        CleanedMetrics.SENSAI_ideal_power_target_db(FileNumber) = g.visualization_metrics.ideal_power_target_db;
                     else
                        warning('GEDAI ran but visualization_metrics not found for %s\n', FileName);
                    end
                end
            end
            if isfield(EEG.RELAX_Metrics, 'Raw')
                if isfield(EEG.RELAX_Metrics.Raw,'BlinkAmplitudeRatio')
                    RawMetrics.BlinkAmplitudeRatio(1:size(EEG.RELAX_Metrics.Raw.BlinkAmplitudeRatio,1),FileNumber)=EEG.RELAX_Metrics.Raw.BlinkAmplitudeRatio;
                    RawMetrics.BlinkAmplitudeRatio(RawMetrics.BlinkAmplitudeRatio==0)=NaN;
                end
                if isfield(EEG.RELAX_Metrics.Raw,'MeanMuscleStrengthFromOnlySuperThresholdValues')
                    RawMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues(FileNumber)=EEG.RELAX_Metrics.Raw.MeanMuscleStrengthFromOnlySuperThresholdValues; 
                    RawMetrics.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel(FileNumber)=EEG.RELAX_Metrics.Raw.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel;
                end
            end   
        end
    end

    %% Record warnings about potential issues:
    EEG.RELAX_issues_to_check.aFileName=cellstr(FileName);
    if size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1)>RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted*size(EEG.allchan,2)
        EEG.RELAX_issues_to_check.PREP_rejected_too_many_electrodes=size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1); % 1.1.4: fix dimension specification error
    else
        EEG.RELAX_issues_to_check.PREP_rejected_too_many_electrodes=0;
    end
    % The impedance screen enforces the cap itself, so its rejection count can never
    % exceed the budget. What is worth flagging is the cap having FIRED - that means more
    % channels were over the limit than were rejected, and the file needs an eye:
    EEG.RELAX_issues_to_check.Impedance_screen_unavailable=0;
    EEG.RELAX_issues_to_check.Impedance_rejected_too_many_electrodes=0;
    if isfield(EEG,'RELAX') && isfield(EEG.RELAX,'ImpedanceScreen')
        if ~EEG.RELAX.ImpedanceScreen.available
            EEG.RELAX_issues_to_check.Impedance_screen_unavailable=1;
        elseif EEG.RELAX.ImpedanceScreen.capApplied
            EEG.RELAX_issues_to_check.Impedance_rejected_too_many_electrodes=EEG.RELAX.ImpedanceScreen.nOverThreshold;
        end
    end
    % The bridge screen enforces no cap on itself (see RELAX_bridged_channels), so
    % there is nothing analogous to "rejected too many" here - just the two QC flags it
    % raises about the recording's own reconstructability:
    EEG.RELAX_issues_to_check.Bridge_screen_unavailable=0;
    EEG.RELAX_issues_to_check.Bridge_oversized_cluster=0;
    EEG.RELAX_issues_to_check.Bridge_excessive_bridging=0;
    if isfield(EEG,'RELAX') && isfield(EEG.RELAX,'BridgeScreen')
        if ~EEG.RELAX.BridgeScreen.available
            EEG.RELAX_issues_to_check.Bridge_screen_unavailable=1;
        else
            if EEG.RELAX.BridgeScreen.oversizedCluster
                EEG.RELAX_issues_to_check.Bridge_oversized_cluster=EEG.RELAX.BridgeScreen.maxClusterSize;
            end
            if EEG.RELAX.BridgeScreen.excessiveBridging
                EEG.RELAX_issues_to_check.Bridge_excessive_bridging=EEG.RELAX.BridgeScreen.nBridgedChannels;
            end
        end
    end
    % Deliberately omits ImpedanceBasedChannelToReject: it already enforces
    % MaxProportionOfElectrodesThatCanBeDeleted on itself (see the impedance block
    % above) and reports separately via Impedance_rejected_too_many_electrodes when its
    % own cap fires. This check instead totals the screens with no cap of their own
    % (muscle/noise recommendations, PREP, bridge) against the same budget - it is not a
    % true combined-deletion total across every screen, so a file where several capped
    % and uncapped screens each reject independently can still lose more of its montage
    % than this check alone would suggest.
    if (EEG.RELAXProcessingExtremeRejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete...
            +EEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete...
            +size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1)...
            +size(EEG.RELAXProcessingExtremeRejections.BridgedChannelToReject,1))...
            >=RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted*size(EEG.allchan,2)
        EEG.RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold=...
            (EEG.RELAXProcessingExtremeRejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete...
            +EEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete...
            +size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1)...
            +size(EEG.RELAXProcessingExtremeRejections.BridgedChannelToReject,1));
    else
        EEG.RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold=0;
    end
    if EEG.RELAXProcessingExtremeRejections.ProportionExcludedForExtremeOutlier>0.20
        EEG.RELAX_issues_to_check.HighProportionExcludedAsExtremeOutlier=EEG.RELAXProcessingExtremeRejections.ProportionExcludedForExtremeOutlier;
    else 
        EEG.RELAX_issues_to_check.HighProportionExcludedAsExtremeOutlier=0;
    end
    if isfield(EEG.RELAX, 'IQRmethodDetectedBlinks') % if IQRmethodDetectedBlinks has been created, do the following (thanks to Jane Tan for the suggested bug fix when IQRmethodDetectedBlinks is not created)
        EEG.RELAX_issues_to_check.NoBlinksDetected=(EEG.RELAX.IQRmethodDetectedBlinks==0);
    end
    if RELAX_cfg.Do_MWF_Once==1
        EEG.RELAX_issues_to_check.MWF_eigenvector_deficiency_R1=isa(EEG.RELAXProcessingRoundOne.RankDeficiency,'char');
    end
    if RELAX_cfg.Do_MWF_Twice==1
        EEG.RELAX_issues_to_check.MWF_eigenvector_deficiency_R2=isa(EEG.RELAXProcessingRoundTwo.RankDeficiency,'char');
    end
    if RELAX_cfg.Do_MWF_Thrice==1
        EEG.RELAX_issues_to_check.MWF_eigenvector_deficiency_R3=isa(EEG.RELAXProcessingRoundThree.RankDeficiency,'char');
    end
    if RELAX_cfg.Perform_wICA_on_ICLabel==1
        if EEG.RELAXProcessing_wICA.Proportion_artifactICs_reduced_by_wICA>0.80
            EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=EEG.RELAXProcessing_wICA.Proportion_artifactICs_reduced_by_wICA;
        else
            EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=0;
        end
        EEG.RELAX_issues_to_check.DataMaybeTooShortForValidICA = EEG.RELAXProcessing_wICA.DataMaybeTooShortForValidICA;
        EEG.RELAX_issues_to_check.fastica_symm_Didnt_Converge=EEG.RELAXProcessing_wICA.fastica_symm_Didnt_Converge(1,3);
    end
    if RELAX_cfg.Perform_ICA_subtract==1
        if EEG.RELAXProcessing_ICA.Proportion_artifactICs_reduced_by_ICA>0.80
            EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=EEG.RELAXProcessing_ICA.Proportion_artifactICs_reduced_by_ICA;
        else
            EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=0;
        end
        EEG.RELAX_issues_to_check.DataMaybeTooShortForValidICA = EEG.RELAXProcessing_ICA.DataMaybeTooShortForValidICA;
        EEG.RELAX_issues_to_check.fastica_symm_Didnt_Converge=EEG.RELAXProcessing_ICA.fastica_symm_Didnt_Converge(1,3);
    end

    %% Interpolate channels after cleaning
    if strcmp(RELAX_cfg.InterpolateRejectedElectrodesAfterCleaning,'yes')
        EEG = pop_interp(EEG, EEG.allchan, 'spherical');
    end

    %% Add MNE-readable TMS markers at trainIdx boundaries
    % MNE's Annotations model has no room for a custom 'trainIdx' struct field on a
    % boundary event, so it's silently dropped on import. Insert a companion 'TMS_train_N'
    % event at the same latency so the train identity survives into MNE as event.description.
    isBoundaryEvent = arrayfun(@(x) (ischar(x.type) || isstring(x.type)) && strcmp(x.type, 'boundary'), EEG.event);
    trainIdxBoundaries = find(isBoundaryEvent & arrayfun(@(x) isfield(x, 'trainIdx') && ~isempty(x.trainIdx), EEG.event));
    for k = 1:numel(trainIdxBoundaries)
        evIdx = trainIdxBoundaries(k);
        newEvent = EEG.event(evIdx);
        newEvent.type = sprintf('TMS_train_%d', EEG.event(evIdx).trainIdx);
        newEvent.duration = 0;
        EEG.event(end+1) = newEvent;
    end
    if ~isempty(trainIdxBoundaries)
        EEG = eeg_checkset(EEG, 'eventconsistency');
    end

    %% SAVE FILE:
    if ~exist([RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep 'Cleaned_Data'], 'dir')
        mkdir([RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep 'Cleaned_Data'])
    end
    if RELAX_cfg.Run_GEDAI==1
        % Include GEDAI in savefile name and save artifacts file. IMK added
        savingfilename = '_RELAX_GEDAI.set';
        
         if ~exist([RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep 'GEDAI_artifacts'], 'dir')
            mkdir([RELAX_cfg.foldername, filesep 'RELAXProcessed' filesep 'GEDAI_artifacts'])
         end
         
        SaveSet_ArtifactsFile =[RELAX_cfg.foldername,filesep 'RELAXProcessed' filesep 'GEDAI_artifacts', filesep FileName '_RELAX_GEDAI_artifacts.set'];
        EEGartifacts = pop_saveset( EEGartifacts, SaveSet_ArtifactsFile ); 
    else
        savingfilename = '_RELAX.set';
    end
    % Save clean EEG
    SaveSet_CleanedFile =[RELAX_cfg.foldername,filesep 'RELAXProcessed' filesep 'Cleaned_Data', filesep FileName savingfilename];
    EEG.RELAX_settings_used_to_clean_this_file=RELAX_cfg;
    EEG = pop_saveset( EEG, SaveSet_CleanedFile );
    
    % Record warnings for all participants in single table:
    % Via RELAX_append_table_row for the same reason as the extreme-rejections table
    % above: a checkpoint written before a column was added has a different variable set.
    % The try/catch below is the original fallback, which diverted the row into
    % RELAX_issues_to_check_2nd_run on any assignment failure. That path is worth NOT
    % reaching on a mere schema change, because RELAX_issues_to_check_2nd_run is saved to
    % the SAME file as RELAX_issues_to_check (see below) with a plain save() - so taking
    % it overwrites the checkpoint and loses every row the table had already accumulated.
    % Reconciling the schema keeps those rows; the catch remains for anything else.
    try
        if ~exist('RELAX_issues_to_check','var')
            RELAX_issues_to_check = [];
        end
        RELAX_issues_to_check = RELAX_append_table_row(RELAX_issues_to_check, FileNumber, ...
            struct2table(EEG.RELAX_issues_to_check,'AsArray',true));
    catch
        RELAX_issues_to_check_2nd_run(FileNumber,:) = struct2table(EEG.RELAX_issues_to_check,'AsArray',true);
        warning('The variable: "RELAX_issues_to_check" already exists and includes different settings from your current settings')
        warning('This is likely from a previous run of RELAX. Saving variable as "RELAX_issues_to_check_2nd_run" instead');
    end
    
    %% Save statistics for each participant and across participants, graph cleaning metrics:
    
    % Also set empty output variables in case these are not produced because certain
    % parameters have been switched off:

    %savefileone=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'RELAXProcessingExtremeRejectionsAllParticipants'];
    save(RELAXProcessingExtremeRejectionsAllParticipantsFile, 'RELAXProcessingExtremeRejectionsAllParticipants')
    if RELAX_cfg.Do_MWF_Once==1
        %savefileone=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatisticsRoundOne'];
        save(ProcessingStatisticsFile,'RELAXProcessingMWFStepOneAllParticipants')
    else
        RELAXProcessingMWFStepOneAllParticipants={};
    end
    if RELAX_cfg.Do_MWF_Twice==1
        %savefiletwo=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatisticsRoundTwo'];
        save(ProcessingStatisticsFile,'RELAXProcessingMWFStepTwoAllParticipants')
    else
        RELAXProcessingMWFStepTwoAllParticipants={};
    end
    if RELAX_cfg.Do_MWF_Thrice==1
        %savefilethree=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatisticsRoundThree'];
        save(ProcessingStatisticsFile,'RELAXProcessingMWFStepThreeAllParticipants')
    else
        RELAXProcessingMWFStepThreeAllParticipants={};
    end
    if RELAX_cfg.Perform_wICA_on_ICLabel==1 || RELAX_cfg.Perform_targeted_wICA==1
        %savefilefour=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatistics_wICA'];
        save(ProcessingStatistics_wICAFile,'RELAXProcessing_wICA_AllParticipants')
    else
        RELAXProcessing_wICA_AllParticipants={}; 
    end
    if RELAX_cfg.Perform_ICA_subtract==1
        %savefilefour=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatistics_ICA'];
        save(ProcessingStatistics_ICAFile,'RELAXProcessing_ICA_AllParticipants')
    else
        RELAXProcessing_ICA_AllParticipants={}; 
    end
    if exist('CleanedMetrics','var')
        %savemetrics=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'CleanedMetrics'];
        save(CleanedMetricsFile,'CleanedMetrics')
    else
        CleanedMetrics={};
    end
    if exist('RawMetrics','var')
        %savemetrics=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'RawMetrics'];
        save(RawMetricsFile,'RawMetrics')
    else
        RawMetrics={};
    end
    if exist('RELAX_issues_to_check','var')
        %savemetrics=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'RELAX_issues_to_check'];
        save(RELAX_issues_to_checkFile,'RELAX_issues_to_check')
    end
    if exist('RELAX_issues_to_check_2nd_run','var')
        %savemetrics=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'RELAX_issues_to_check_2nd_run'];
        save(RELAX_issues_to_checkFile,'RELAX_issues_to_check_2nd_run')
    end
    % This file is cleaned: if an earlier run recorded it as unprocessed, retract that.
    if exist('RELAX_unprocessed_files','var') && istable(RELAX_unprocessed_files) && ...
            any(strcmp(RELAX_unprocessed_files.File, RELAX_cfg.filename))
        RELAX_unprocessed_files(strcmp(RELAX_unprocessed_files.File, RELAX_cfg.filename), :) = [];
        save(RELAX_unprocessed_filesFile, 'RELAX_unprocessed_files');
        writetable(RELAX_unprocessed_files, strrep(RELAX_unprocessed_filesFile, '.mat', '.csv'));
    end

    RELAX_cfg.filename=[];
    RELAX_cfg_file = fullfile(statsDir, sprintf('RELAX_cfg_%s.mat', tag));  % <- Add this to use tagged filename
    save(RELAX_cfg_file,'RELAX_cfg')      

    catch fileErr
        % A file this run cannot clean. Two kinds arrive here: a deliberate exclusion
        % raised with identifier RELAX:excludeRecording (a recording judged unusable,
        % e.g. mostly-bridged - see the bridge screen above), and any other MATLAB
        % error. Both are logged in full, recorded in the per-participant
        % RELAX_unprocessed_files table (.mat and .csv next to the other
        % PerParticipantStats), and the loop moves on to the participant's next file.
        % No cleaned .set and no aggregate rows are written for this file, so it is
        % absent from the cleaned data rather than present and wrong.
        if strcmp(fileErr.identifier, 'RELAX:excludeRecording')
            outcome = 'excluded';
        else
            outcome = 'error';
        end
        fprintf(2, '\n==== %s: %s (%s) ====\n%s\n', RELAX_cfg.filename, outcome, ...
            fileErr.identifier, getReport(fileErr, 'extended', 'hyperlinks', 'off'));
        warning('RELAX:fileNotProcessed', '%s not processed (%s): %s', ...
            RELAX_cfg.filename, outcome, fileErr.message);
        failedRow = table({RELAX_cfg.filename}, {outcome}, {fileErr.identifier}, ...
            {fileErr.message}, datetime('now'), ...
            'VariableNames', {'File', 'Outcome', 'Identifier', 'Message', 'When'});
        if ~exist('RELAX_unprocessed_files','var') || ~istable(RELAX_unprocessed_files)
            RELAX_unprocessed_files = failedRow;
        else
            RELAX_unprocessed_files(strcmp(RELAX_unprocessed_files.File, RELAX_cfg.filename), :) = [];
            RELAX_unprocessed_files = [RELAX_unprocessed_files; failedRow];
        end
        save(RELAX_unprocessed_filesFile, 'RELAX_unprocessed_files');
        writetable(RELAX_unprocessed_files, strrep(RELAX_unprocessed_filesFile, '.mat', '.csv'));
        continue
    end
end

% If every file of this participant ended in the catch above, none of the per-file
% tables exists; return them empty rather than failing on an unassigned output after
% the failures have already been recorded.
for emptyOut = {'CleanedMetrics', 'RawMetrics', 'RELAXProcessingMWFStepOneAllParticipants', ...
        'RELAXProcessingMWFStepTwoAllParticipants', 'RELAXProcessing_wICA_AllParticipants', ...
        'RELAXProcessing_ICA_AllParticipants', 'RELAXProcessingMWFStepThreeAllParticipants', ...
        'RELAXProcessingExtremeRejectionsAllParticipants'}
    if ~exist(emptyOut{1}, 'var'); eval([emptyOut{1} ' = {};']); end
end
if ~exist('RELAX_issues_to_check', 'var'); RELAX_issues_to_check = table(); end

set(groot, 'defaultAxesTickLabelInterpreter','none');
if RELAX_cfg.computecleanedmetrics==1
    % Create output directory for figures
    fig_dir = fullfile(RELAX_cfg.myPath, 'RELAXProcessed', 'QualityMetrics_Figures');
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
    % Plot QC Metrics
    try
        fig1 = figure('Name','BlinkAmplitudeRatio','units','normalized','outerposition',[0.05 0.05 0.95 0.95]);
        boxplot(CleanedMetrics.BlinkAmplitudeRatio);
        xticklabels(RELAX_cfg.files); xtickangle(90);
        set(gca,'FontSize',16, 'FontWeight', 'bold') % Creates an axes and sets its FontSize to 21
        % Save figure
        saveas(fig1, fullfile(fig_dir, 'BlinkAmplitudeRatio.png'));
        %saveas(fig1, fullfile(fig_dir, 'BlinkAmplitudeRatio.fig'));
    catch
    end
    try
        fig2 = figure('Name','MeanMuscleStrengthFromOnlySuperThresholdValues','units','normalized','outerposition',[0.05 0.05 0.95 0.95]);
        b=bar(CleanedMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues); 
        xtickangle(90); xticks([1:1:size(RELAX_cfg.files,2)]); b(1).BaseValue = RELAX_cfg.MuscleSlopeThreshold;
        xticklabels(RELAX_cfg.files); ylim([RELAX_cfg.MuscleSlopeThreshold max(CleanedMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues)+1]);b.ShowBaseLine='off';
        set(gca,'FontSize',16, 'FontWeight', 'bold') % Creates an axes and sets its FontSize to 21
        % Save figure
        saveas(fig2, fullfile(fig_dir, 'MeanMuscleStrengthFromOnlySuperThresholdValues.png'));
        %saveas(fig2, fullfile(fig_dir, 'MeanMuscleStrengthFromOnlySuperThresholdValues.fig'));
    catch
    end
    try
        fig3 = figure('Name','ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel','units','normalized','outerposition',[0.05 0.05 0.95 0.95]);
        bar(CleanedMetrics.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel);
        set(gca,'FontSize',16, 'FontWeight', 'bold') % Creates an axes and sets its FontSize to 21
        xtickangle(90); xticks([1:1:size(RELAX_cfg.files,2)]);
        xticklabels(RELAX_cfg.files);
        % Save figure
        saveas(fig3, fullfile(fig_dir, 'ProportionEpochsWithMuscleAboveThresholdAnyChannel.png'));
        %saveas(fig3, fullfile(fig_dir, 'ProportionEpochsWithMuscleAboveThresholdAnyChannel.fig'));
    catch
    end
end

clearvars -except 'RELAX_cfg' 'FileNumber' 'CleanedMetrics' 'RawMetrics' 'RELAXProcessingMWFStepOneAllParticipants' 'RELAXProcessingMWFStepTwoAllParticipants' 'RELAXProcessing_wICA_AllParticipants'...
        'RELAXProcessing_ICA_AllParticipants' 'RELAXProcessingMWFStepThreeAllParticipants' 'Warning' 'RELAX_issues_to_check' 'RELAX_issues_to_check_2nd_run'...
        'RELAXProcessingExtremeRejectionsAllParticipants' 'WarningAboutFileNumber' 'RELAX_unprocessed_files';

if ~exist('RELAX_issues_to_check_2nd_run','var')    
    warning('Check "RELAX_issues_to_check" to see if any issues were noted for specific files');
    RELAX_issues_to_check_2nd_run=['This variable is only here as a placeholder in case you have already run RELAX once,...' ...
        ' and had previously left the output variables in the same folder as the folder where you have saved the currently cleaned data.' ...
        'This was not an issue with the current run, so this placeholder variable was not filled'];
elseif exist('RELAX_issues_to_check_2nd_run','var')
    warning('Check "RELAX_issues_to_check_2nd_run" to see if any issues were noted for specific files');
end

if WarningAboutFileNumber==1
    warning('You instructed RELAX to clean more files than were in your data folder. Check all your expected files were there?');
end

if RELAX_cfg.ProbabilityDataHasNoBlinks<2 && ismember('NoBlinksDetected', RELAX_issues_to_check.Properties.VariableNames) && sum(RELAX_issues_to_check.NoBlinksDetected)>1
    f = msgbox('RELAX did not detect any blinks for some files. Open the "RELAX_issues_to_check" struct in the workspace to check which files. We recommend visually inspecting these files to ensure there has not been an error.'...
    ,'No blinks detected for some files');    
    set(f,'Position',[500,500,450,100]);
    ah = get( f, 'CurrentAxes' );
    ch = get( ah, 'Children' );
    set( ch, 'FontSize', 12 ); %makes text bigger
end

if ismember('ElectrodeRejectionRecommendationsMetOrExceededThreshold', RELAX_issues_to_check.Properties.VariableNames) && any(RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold > 0)
    f = msgbox('Some files met or exceeded the electrode rejection threshold. We recommend visually inspecting the raw and cleaned files where this is the case. Open the "RELAX_issues_to_check" struct in the workspace, and check the third column. Files that exceeded the threshold will show a value above 0. Exclude files where raw data seems irretrievably noisy, or cleaned data still contains excessive noise.'...
    ,'Some files met or exceeded the electrode rejection threshold');    
    set(f,'Position',[300,300,450,150]);
    ah = get( f, 'CurrentAxes' );
    ch = get( ah, 'Children' );
    set( ch, 'FontSize', 12 ); %makes text bigger
end

toc

% Turn auto fig display back on
set(0, 'DefaultFigureVisible', 'on')

%% POTENTIAL IMPROVEMENTS THAT COULD BE MADE:
% 1) work out a way to threshold horizontal eye movements so the script
% catches the onset and offset, rather than just the absolute +/-2MAD on
% opposite sides of the head?
% 2) adding a requirement that the IQR blink detection method detects that
% positive amplitude shifts are biased towards frontal electrodes?
