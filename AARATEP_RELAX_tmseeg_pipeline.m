function EEG = AARATEP_RELAX_tmseeg_pipeline(EEG, config)
% AARATEP_RELAX_tmseeg_pipeline - Combined AARATEP and RELAX pipeline for TMS-EEG preprocessing
%
% Inputs:
%   EEG: EEGLAB format EEG structure
%   config: Configuration structure with pipeline settings
%
% Output:
%   EEG: Preprocessed EEG structure
%
% This function combines steps from the AARATEP pipeline (c_TMSEEG_Preprocess_AARATEPPipeline)
% and the RELAX pipeline (RELAX_Wrapper) for preprocessing continuous EEG data recorded
% during TMS treatments.

%% Add dependencies to path
persistent pathModified
if isempty(pathModified)
    mfilepath = fileparts(which(mfilename));
    addpath(fullfile(mfilepath, './AARATEPPipeline/Common'));
    addpath(fullfile(mfilepath, './AARATEPPipeline/Common/EEGAnalysisCode'));
    pathModified = true;
end
PrepFolderLocation = "/home/nga4004/gitRepos/HAPPE/Packages/eeglab2022.0/plugins/PrepPipeline0.57.0";

%% Set default AARATEP parameters if not provided
if ~isfield(config, 'artifactTimespan')
    config.artifactTimespan = [-0.004, 0.012]; % in seconds
end
if ~isfield(config, 'filterPrePostExtrapolationDurations')
    config.filterPrePostExtrapolationDurations = [0 0];
end
if ~isfield(config, 'doDecayRemovalPerTrial')
    config.doDecayRemovalPerTrial = false; % For continuous data
end
if ~isfield(config, 'burstMaxIPI')
    config.burstMaxIPI = 0.022; % 22 ms
end
if ~isfield(config, 'prePostFitDurations')
    config.prePostFitDurations = [30 30]*1e-3; % 30 ms
end
if ~isfield(config, 'maxTau')
    config.maxTau = 20e-3; % 20 ms
end
if ~isfield(config, 'epochWindow')
    config.epochWindow = [-50 50]*1e-3;
end
if ~isfield(config, 'lastPulse_epochWindow')
    config.lastPulse_epochWindow = [-50 150]*1e-3;
end

if ~isfield(config, 'doDebug')
    config.doDebug = false; % Default to no debug output
end

if ~isfield(config, 'eventsToRemove')
    config.eventsToRemove = []; % Default: no events to remove. 
    % Format: struct array with fields 'train' and 'triplet', e.g.:
    %   config.eventsToRemove(1).train = 60; config.eventsToRemove(1).triplet = 11;
end

%% Set default RELAX parameters if not provided
if isfield(EEG, 'filename') && ~isempty(EEG.filename)
    EEG.RELAX.ProcessedFile = EEG.filename;
end

%% Check for TMS events and detect if missing
c_say('Checking for TMS events');
hasTMSEvents = false;
if ~isempty(EEG.event)
    eventTypes = {EEG.event.type};
    hasTMSEvents = any(strcmp(eventTypes, 'TMS'));
end

if ~hasTMSEvents
    c_saySingle('No TMS events found. Attempting automatic detection...');
    
    % Set default parameters for TMS pulse detection if not provided
    if ~isfield(config, 'maxPulseRate')
        config.maxPulseRate = 53; % Hz
    end
    if ~isfield(config, 'minNumPulses')
        config.minNumPulses = 1700; % Set to expected number if known
    end
    if ~isfield(config, 'filterMethod')
        config.filterMethod = 'median';
    end
    if ~isfield(config, 'minPulseThreshold')
        config.minPulseThreshold = 60;
    end
    if ~isfield(config, 'maxPulseThreshold')
        config.maxPulseThreshold = 300;
    end
    
    % Attempt to find TMS pulses
    EEG = c_TMSEEG_findTMSPulses(EEG,...
        'filterMethod', config.filterMethod,...
        'maxPulseRate', config.maxPulseRate,...
        'addEventsOfType', 'TMS',...
        'doRefineOnsetTimes', false,...
        'minNumPulses', config.minNumPulses,...
        'minPulseThreshold', config.minPulseThreshold,...
        'maxPulseThreshold', config.maxPulseThreshold);
    
    % Verify that pulses were found
    if ~isempty(EEG.event)
        eventTypes = {EEG.event.type};
        numTMSEvents = sum(strcmp(eventTypes, 'TMS'));
        if numTMSEvents > 0
            c_saySingle('Successfully detected %d TMS pulses', numTMSEvents);
        else
            error('Failed to detect TMS pulses. Please check your data or add TMS events manually.');
        end
    else
        error('Failed to detect TMS pulses. Please check your data or add TMS events manually.');
    end
else
    eventTypes = {EEG.event.type};
    numTMSEvents = sum(strcmp(eventTypes, 'TMS'));
    c_saySingle('Found %d existing TMS events', numTMSEvents);
end
c_sayDone();

% Pre-step: Label TMS events with train and triplet indices
c_say('Labelling TMS events with train and triplet indices');

tmsEventIndices = find(strcmp({EEG.event.type}, 'TMS'));
numTMSEvents = length(tmsEventIndices);
tmsLatencies = [EEG.event(tmsEventIndices).latency] / EEG.srate; % in seconds

% Define IPI thresholds (in seconds)
withinTripletMaxIPI = 0.025;  % < 25 ms = within triplet (actual is ~20 ms)
withinTrainMaxIPI  = 0.250;   % < 250 ms = within train (actual is ~200 ms)
% Between trains is ~10 s

iPIs = [Inf, diff(tmsLatencies)]; % Inf for first event so it starts a new train

trainIdx   = 0;
tripletIdx = 0;

for ii = 1:numTMSEvents
    ipi = iPIs(ii);
    
    if (ipi > withinTrainMaxIPI)
        % New train
        trainIdx   = trainIdx + 1;
        tripletIdx = 1;
    elseif (ipi > withinTripletMaxIPI)
        % New triplet within same train
        tripletIdx = tripletIdx + 1;
    end
    % else: within-triplet pulse, same train and triplet indices
    
    EEG.event(tmsEventIndices(ii)).train   = trainIdx;
    EEG.event(tmsEventIndices(ii)).triplet = tripletIdx;
end

c_saySingle('Labelled %d TMS events across %d trains', numTMSEvents, trainIdx);
c_sayDone();

%% Remove manually specified TMS events by train/triplet index
if ~isempty(config.eventsToRemove)
    c_say('Removing manually specified TMS events');
    indicesToRemove = false(1, length(EEG.event));
    for iRemove = 1:length(config.eventsToRemove)
        targetTrain   = config.eventsToRemove(iRemove).train;
        targetTriplet = config.eventsToRemove(iRemove).triplet;
        for iEv = tmsEventIndices
            if isfield(EEG.event(iEv), 'train') && isfield(EEG.event(iEv), 'triplet')
                if (EEG.event(iEv).train == targetTrain) && (EEG.event(iEv).triplet == targetTriplet)
                    indicesToRemove(iEv) = true;
                end
            end
        end
    end
    numRemoved = sum(indicesToRemove);
    
    % Before removing events, interpolate the artifact timespan around 
    % the to-be-removed events so TMS artifacts don't remain in the data.
    % Temporarily mark these events with a unique type for targeted interpolation.
    eventsToRemoveIndices = find(indicesToRemove);
    for iEv = eventsToRemoveIndices
        EEG.event(iEv).type = 'TMS_toRemove';
    end
    EEG = eeg_checkset(EEG);
    
    % Interpolate artifact timespan around the marked events
    EEG = c_EEG_ReplaceEpochTimeSegment(EEG,...
        'timespanToReplace', [config.artifactTimespan(1), 3*config.artifactTimespan(2)],...
        'method', 'ARExtrapolation',...
        'prePostFitDurations', config.prePostFitDurations,...
        'eventType', 'TMS_toRemove');
    
    % Now remove the marked events
    indicesToRemove = strcmp({EEG.event.type}, 'TMS_toRemove');
    EEG.event(indicesToRemove) = [];
    EEG = eeg_checkset(EEG);
    c_saySingle('Removed %d TMS events matching specified train/triplet criteria (artifacts interpolated)', numRemoved);
    
    % Renumber trains and triplets to remove gaps
    tmsEventIndices = find(strcmp({EEG.event.type}, 'TMS'));
    if ~isempty(tmsEventIndices)
        oldTrains = [EEG.event(tmsEventIndices).train];
        uniqueTrains = unique(oldTrains, 'stable');
        for iNewTrain = 1:length(uniqueTrains)
            oldTrain = uniqueTrains(iNewTrain);
            inThisTrain = tmsEventIndices(oldTrains == oldTrain);
            oldTriplets = [EEG.event(inThisTrain).triplet];
            uniqueTriplets = unique(oldTriplets, 'stable');
            for iNewTriplet = 1:length(uniqueTriplets)
                oldTriplet = uniqueTriplets(iNewTriplet);
                inThisTriplet = inThisTrain(oldTriplets == oldTriplet);
                for iEv = inThisTriplet
                    EEG.event(iEv).train   = iNewTrain;
                    EEG.event(iEv).triplet = iNewTriplet;
                end
            end
        end
        c_saySingle('Renumbered TMS events: %d trains remaining', length(uniqueTrains));
    end
    
    c_sayDone();
end

%% Initialize debug tracking
if config.doDebug
    intermediateEEGs = {};
    intermediateEEGs{1} = EEG; % Store original EEG as first step
    intermediateLabels = {};
    intermediateLabels{1} = 'Original EEG';
end

%% AARATEP Pipeline Steps

% Step 1: Handle burst events - remove pulses that are too close together (~20ms)
c_say('Handling burst events');
EEG = c_TMSEEG_handleBurstEvents(EEG,...
    'pulseEvents', {'TMS'},...
    'method', 'cutIPIKeepFirst',...
    'burstMaxIPI', config.burstMaxIPI);
c_sayDone();

if config.doDebug
    intermediateEEGs{end+1} = EEG;
    intermediateLabels{end+1} = 'After burst handling';
end

% Step 2: Interpolate artifact timespan (first time)
c_say('Interpolating artifact timespan (pre-filter)');
EEG = c_EEG_ReplaceEpochTimeSegment(EEG,...
    'timespanToReplace', config.artifactTimespan,...
    'method', 'ARExtrapolation',...
    'prePostFitDurations', config.prePostFitDurations,...
    'eventType', 'TMS');
c_sayDone();

if config.doDebug
    intermediateEEGs{end+1} = EEG;
    intermediateLabels{end+1} = 'Artifact interpolated (pre-filter)';
end

% Step 3: Epoch data, simultaneously fit and remove exponential decay + DC offset per channel/trial
c_say('Epoching and fitting/removing exponential decay + offset (a*exp(-b*x) + c)');

% Store original continuous EEG for reconstruction
originalEEG = EEG;

% Epoch the data around TMS events
EEG = pop_epoch(EEG, {'TMS'}, config.epochWindow, 'epochinfo', 'yes');

% Create a separate epoched struct for the 10th triplet in each train
tmsEventsAll = find(strcmp({originalEEG.event.type}, 'TMS'));
lastPulseLatencies = [];
lastPulseTrains = [];
for iEv = tmsEventsAll
    if isfield(originalEEG.event(iEv), 'triplet') && isfield(originalEEG.event(iEv), 'train')
        if originalEEG.event(iEv).triplet == 10
            lastPulseLatencies(end+1) = originalEEG.event(iEv).latency;
            lastPulseTrains(end+1) = originalEEG.event(iEv).train;
        end
    end
end

if ~isempty(lastPulseLatencies)
    % Temporarily add a unique marker for the 10th triplet events
    originalEEG_tmp = originalEEG;
    for iLat = 1:length(lastPulseLatencies)
        newIdx = length(originalEEG_tmp.event) + 1;
        originalEEG_tmp.event(newIdx).type    = 'TMS_triplet10';
        originalEEG_tmp.event(newIdx).latency = lastPulseLatencies(iLat);
        originalEEG_tmp.event(newIdx).duration = 0;
        originalEEG_tmp.event(newIdx).train   = lastPulseTrains(iLat);
    end
    originalEEG_tmp = eeg_checkset(originalEEG_tmp);
    EEG_lastpulse = pop_epoch(originalEEG_tmp, {'TMS_triplet10'}, config.lastPulse_epochWindow, 'epochinfo', 'yes');
    c_saySingle('Created EEG_lastpulse with %d epochs (10th triplet per train)', EEG_lastpulse.trials);
else
    EEG_lastpulse = [];
    c_saySingle('No 10th-triplet events found; EEG_lastpulse is empty');
end

% %% For Testing: Take only the last 100 epochs to speed up processing during development
% if EEG.trials > 100
%     c_saySingle('Keeping only the last 100 epochs for testing');
%     EEG = pop_select(EEG, 'trial', (EEG.trials-99):EEG.trials);
% end
% if ~isempty(EEG_lastpulse) && EEG_lastpulse.trials > 10
%     c_saySingle('Keeping only the last 10 epochs of EEG_lastpulse for testing');
%     EEG_lastpulse = pop_select(EEG_lastpulse, 'trial', (EEG_lastpulse.trials-9):EEG_lastpulse.trials);
    
%     % Subsample parallel arrays that track original last pulse metadata
%     lastPulseLatencies = lastPulseLatencies(end-9:end);
%     lastPulseTrains = lastPulseTrains(end-9:end);
% end

if config.doDebug
    intermediateEEGs{end+1} = {EEG, EEG_lastpulse};
    intermediateLabels{end+1} = 'Epoched EEG';
end

% Fit and remove exponential decay + offset simultaneously
[EEG_lastpulse, decayMisc] = c_TMSEEG_fitAndRemoveDecayAndOffset(EEG_lastpulse, ...
    'artifactTimespan',   config.artifactTimespan, ...
    'fitTimespan',        config.artifactTimespan(2) + [0, 2.0], ...
    'doPlot',             false);

fitParamsToPass = nan(EEG.nbchan, EEG.trials, 7);
if ~isempty(EEG_lastpulse)
    for iTrial = 1:EEG.trials
        tr = [];
        latencies = EEG.epoch(iTrial).eventlatency;
        if iscell(latencies)
            zeroIdx = find(cellfun(@(x) isequal(x, 0), latencies), 1);
        else
            zeroIdx = find(latencies == 0, 1);
        end
        if ~isempty(zeroIdx) && isfield(EEG.epoch(iTrial), 'eventtrain')
            if iscell(EEG.epoch(iTrial).eventtrain)
                tr = EEG.epoch(iTrial).eventtrain{zeroIdx};
            else
                tr = EEG.epoch(iTrial).eventtrain(zeroIdx);
            end
        end
        
        if ~isempty(tr)
            matchIdx = find(lastPulseTrains == tr, 1);
            if ~isempty(matchIdx) && matchIdx <= size(decayMisc.fitParams, 2)
                fitParamsToPass(:, iTrial, :) = decayMisc.fitParams(:, matchIdx, :);
            end
        end
    end
end

[EEG, ~] = c_TMSEEG_fitAndRemoveDecayAndOffset(EEG, ...
    'artifactTimespan',   config.artifactTimespan, ...
    'fitParams',          fitParamsToPass, ...
    'steadyStateLength',  30, ...
    'expWin',            0, ...
    'doPlot',             false);

c_sayDone();

if config.doDebug
    intermediateEEGs{end+1} = {EEG, EEG_lastpulse};
    intermediateLabels{end+1} = 'After decay+offset removal';
end

% Step 4: Apply modified highpass filter to remove any residual low-frequency artifact (using pre/post extrapolation to prevent edge artifacts)
EEG_lastpulse = c_TMSEEG_applyModifiedBandpassFilter(EEG_lastpulse,...
	'lowCutoff', 1.0,...
    'prePostFitDurations', [50, 50]*1e-3,...
    'piecewiseTimeToExtend', 2.0,...
    'prePostExtrapolationDurations',  [2.0, 2.0],...
	'artifactTimespan', config.artifactTimespan*2);
c_sayDone();
if config.doDebug
	intermediateEEGs{end+1} = EEG_lastpulse;
	intermediateLabels{end+1} = 'High-pass filtered';
end

EEG = c_EEG_epochs2continuous(EEG, EEG_lastpulse, originalEEG, config.artifactTimespan(1));
c_sayDone();

if config.doDebug
	intermediateEEGs{end+1} = EEG;
	intermediateLabels{end+1} = 'reconstructed continuous';
end

c_say('Interpolating artifact timespan (post-decay)');
EEG = c_EEG_ReplaceEpochTimeSegment(EEG,...
    'timespanToReplace', [config.artifactTimespan(1),...
                          config.artifactTimespan(2) + 1e-3],...
    'method', 'ARExtrapolation',...
    'prePostFitDurations', config.prePostFitDurations,...
    'eventType', 'TMS');

c_sayDone();

if config.doDebug
    intermediateEEGs{end+1} = EEG;
    intermediateLabels{end+1} = 'Artifact interpolated (post-decay)';
end


%% RELAX Pipeline Steps (starting from notch filter)

config.ms_per_sample = (1000/EEG.srate);

% Initialize minimal RELAX fields needed for processing
EEG.RELAX.Data_has_been_averagerereferenced = 0;
EEG.RELAX.Data_has_been_cleaned = 0;

% Store all channels before any rejections
EEG.allchan = EEG.chanlocs;

% Notch filter data (lines 157-175 from RELAX_Wrapper)
if strcmp(config.NotchFilterType, 'Butterworth')
    if length(config.LineNoiseFrequency) > 1
        for freq = config.LineNoiseFrequency
            EEG = RELAX_filtbutter(EEG, freq-5, freq+5, 4, 'bandstop', 'acausal');
        end
    else
        EEG = RELAX_filtbutter(EEG, config.LineNoiseFrequency-5, config.LineNoiseFrequency+5, 4, 'bandstop', 'acausal');
    end
end

if strcmp(config.NotchFilterType, 'PMnotch')
    EEG = pop_basicfilter(EEG, 1:EEG.nbchan, 'Boundary', 'boundary', 'Cutoff', config.LineNoiseFrequency, 'Design', 'notch', 'Filter', 'PMnotch', 'Order', 180, 'RemoveDC', 'on');
end

if config.doDebug && (strcmp(config.NotchFilterType, 'Butterworth') || strcmp(config.NotchFilterType, 'PMnotch'))
    intermediateEEGs{end+1} = EEG;
    intermediateLabels{end+1} = 'After notch filter';
end

% Downsample if requested
if strcmp(config.DownSample, 'yes')
    EEG = pop_resample(EEG, config.DownSample_to_X_Hz);
    config.ms_per_sample = (1000/EEG.srate);
    
    if config.doDebug
        intermediateEEGs{end+1} = EEG;
        intermediateLabels{end+1} = 'Downsampled';
    end
end

if strcmp(config.NotchFilterType, 'ZaplinePlus')
    % Put the vendored, boundary-fixed zapline-plus ahead of the EEGLAB plugin copy.
    % Must come after eeglab(), which prepends its own plugin directories.
    RELAX_add_vendored_zapline();

    [EEG] = clean_data_with_zapline_plus_eeglab_wrapper(EEG, struct('plotResults', 0));
    
    if config.doDebug
        intermediateEEGs{end+1} = EEG;
        intermediateLabels{end+1} = 'After ZaplinePlus';
    end
end

if strcmp(config.LowPassFilterBeforeMWF,'yes') % original implementation, not recommended before MWF unless downsampling, as increases chances of rank deficiencies
    if strcmp(config.FilterType,'Butterworth')
        EEG = RELAX_filtbutter( EEG, config.HighPassFilter, config.LowPassFilter, 4, 'bandpass', config.causal_or_acausal_filter);
        if ~isempty(config.electrodes_2_keep_but_not_clean{1})
            if (config.HighPassFilter_aux_elecs~=config.HighPassFilter) || (config.LowPassFilter_aux_elecs~=config.LowPassFilter)
                disp('Applying filter settings to auxilary electrodes:')
            end
            Non_cleaned_electrodes = RELAX_filtbutter( Non_cleaned_electrodes, config.HighPassFilter_aux_elecs, config.LowPassFilter_aux_elecs, 4, 'bandpass', config.causal_or_acausal_filter);
        end
    end
    if strcmp(config.FilterType,'pop_eegfiltnew')
        EEG = pop_eegfiltnew(EEG,config.HighPassFilter,config.LowPassFilter_aux_elecs);
        if ~isempty(config.electrodes_2_keep_but_not_clean{1})
            if (config.HighPassFilter_aux_elecs~=config.HighPassFilter) || (config.LowPassFilter_aux_elecs~=config.LowPassFilter)
                disp('Applying filter settings to auxilary electrodes:')
            end
            Non_cleaned_electrodes = pop_eegfiltnew(Non_cleaned_electrodes,config.HighPassFilter_aux_elecs,config.LowPassFilter_aux_elecs);
        end
    end
end

if strcmp(config.DownSample,'yes')
    EEG = pop_resample(EEG,config.DownSample_to_X_Hz); % downsample data (if applied, should always be applied after low pass filtering)
    config.ms_per_sample=(1000/EEG.srate);
    if ~isempty(config.electrodes_2_keep_but_not_clean{1})
        Non_cleaned_electrodes = pop_resample(Non_cleaned_electrodes,config.DownSample_to_X_Hz);
    end
end

if config.ms_per_sample<0.7
    warning('The sampling rate for this file is quite high. Depending on your processing power, RELAX may run slowly or even stall, especially if applying MWF cleaning');
    warning('RELAX was validated using 1000Hz sampling rates, which is still a high sample rate for most analyses. You could downsample your data by setting the relevant options in RELAX');
end

if strcmp(config.NotchFilterType,'ZaplinePlus')
    [EEG ] = clean_data_with_zapline_plus_eeglab_wrapper(EEG,struct('plotResults',0)); % requires the zapline plus plugin. Best applied after downsampling to 250Hz or 500Hz.
    if ~isempty(config.electrodes_2_keep_but_not_clean{1})
        % zapline plus won't fix the auxilary electrodes, so Butterworth filtering them instead:
        disp('Applying notch filter to auxilary electrodes as zapline will not work on them:')
        Non_cleaned_electrodes = RELAX_filtbutter( Non_cleaned_electrodes, config.LineNoiseFrequency-3, config.LineNoiseFrequency+3, 4, 'bandstop','acausal');
    end
end

%% Clean flat channels and bad channels showing improbable data:
% PREP pipeline: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4471356/
addpath(genpath(PrepFolderLocation{1,1})); %  1.1.4: fix error where PREP seems to be removed from the path after an EEGLAB update
fncOpts.badTimeThreshold = 0.02;
noisyOut = findNoisyChannels(EEG, fncOpts);  
EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject={};
for x=1:size(noisyOut.noisyChannels.all,2) % loop through output of PREP's findNoisyChannels and take a record of noisy electrodes for deletion:
    PREPBasedChannelToReject{x}=EEG.chanlocs(noisyOut.noisyChannels.all(x)).labels;
    EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject = PREPBasedChannelToReject';
end

% Print channel names to be removed by PREP Noisy Channels
if ~isempty(noisyOut.noisyChannels.all)
    noisy_labels = {EEG.chanlocs(noisyOut.noisyChannels.all).labels};
    fprintf('Removing noisy channels: %s\n', strjoin(noisy_labels, ', '));
else
    fprintf('No noisy channels detected by PREP\n');
end

EEG=pop_select(EEG,'nochannel',noisyOut.noisyChannels.all); % delete noisy electrodes detected by PREP

if config.doDebug
    intermediateEEGs{end+1} = EEG;
    intermediateLabels{end+1} = 'After PREP channel rejection';
end

continuousEEG=EEG;

[continuousEEG, epochedEEG] = RELAX_excluding_channels_and_epoching_iTBS(continuousEEG, config); % Epoch data, detect extremely bad data, delete channels if over the set threshold for proportion of data affected by extreme outlier for each electrode
[continuousEEG, epochedEEG] = RELAX_excluding_extreme_values(continuousEEG, epochedEEG, config); % Mark extreme periods for exclusion from MWF cleaning, and deletion before wICA cleaning

if config.doDebug
    intermediateEEGs{end+1} = continuousEEG;
    intermediateLabels{end+1} = 'After extreme value marking';
end

% Use the continuous data to detect eye blinks and mark
% these in the EEG.event as well as in the mask. The output is
% continuous data but includes all the previous extreme period 
% markings from the epoched data.
if config.ProbabilityDataHasNoBlinks<2
    [continuousEEG, epochedEEG] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, config); % use an IQR threshold method to detect and mark blinks
    if continuousEEG.RELAX.IQRmethodDetectedBlinks(1,1)==0 % If a participants doesn't show any blinks, make a note
        NoBlinksDetected{FileNumber,1}=config.filename; 
        warning('No blinks were detected - if blinks are expected then you should visually inspect the file');
    end
    if config.computerawmetrics==1
    [continuousEEG, epochedEEG] = RELAX_metrics_blinks(continuousEEG, epochedEEG); % record blink amplitude ratio from raw data for comparison.
    end
end

if config.doDebug && config.ProbabilityDataHasNoBlinks<2
    intermediateEEGs{end+1} = continuousEEG;
    intermediateLabels{end+1} = 'After blink detection';
end

% Record extreme artifact rejection details for all participants in single table:
RELAXProcessingExtremeRejectionsAllParticipants = struct2table(epochedEEG.RELAXProcessingExtremeRejections,'AsArray',true);

rawEEG=continuousEEG; % Take a copy of the not yet cleaned data for calculation of all cleaning SER and ARR at the end

%% Mark artifacts for calculating SER and ARR, regardless of whether MWF is performed (RELAX v1.1.3 update): 
if config.computecleanedmetrics==1 && (config.Do_MWF_Once==0 || config.Do_MWF_Twice==0 || config.Do_MWF_Thrice==0)
    [Marking_artifacts_for_SER_ARR, ~] = RELAX_muscle(continuousEEG, epochedEEG, config); 
    Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(Marking_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength==1)=1; 
    [Marking_artifacts_for_SER_ARR] = RELAX_horizontaleye(continuousEEG, config); 
    Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(Marking_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength==1)=1; 
    [Marking_artifacts_for_SER_ARR, ~] = RELAX_drift(continuousEEG, epochedEEG, config); % Use epoched data to add periods showing excessive drift to the mask 
    Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(Marking_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength==1)=1; 
    Marking_all_artifacts_for_SER_ARR.RELAX.NaNsForExtremeOutlierPeriods=continuousEEG.RELAX.NaNsForExtremeOutlierPeriods; 
    [Marking_all_artifacts_for_SER_ARR] = RELAX_pad_brief_mask_periods (Marking_all_artifacts_for_SER_ARR, config, 'notblinks'); % If period has been marked as shorter than config.MinimumArtifactDuration, then pad it out. 
    if isfield(continuousEEG.RELAX,'eyeblinkmask')
            Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength(continuousEEG.RELAX.eyeblinkmask==1)=1; 
            [Marking_all_artifacts_for_SER_ARR] = RELAX_pad_brief_mask_periods (Marking_all_artifacts_for_SER_ARR, config, 'blinks'); 
    end
    continuousEEG.RELAX.NoiseMaskFullLengthR1=Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength; 
    rawEEG.RELAX.NoiseMaskFullLengthR1=Marking_all_artifacts_for_SER_ARR.RELAXProcessing.Details.NoiseMaskFullLength; 
end

if config.saveextremesrejected==1
    if ~exist([config.foldername, filesep 'RELAXProcessed' filesep 'Extremes_Rejected'], 'dir')
        mkdir([config.foldername, filesep 'RELAXProcessed' filesep 'Extremes_Rejected'])
    end
    SaveSetExtremes_Rejected =[config.foldername, filesep 'RELAXProcessed' filesep 'Extremes_Rejected', filesep config.filename '_Extremes_Rejected.set'];    
    EEG = pop_saveset( rawEEG, SaveSetExtremes_Rejected ); % If desired, save data here with bad channels deleted, filtering applied, extreme outlying data periods marked
end

%% THIS SECTION CONTAINS FUNCTIONS WHICH MARK AND CLEAN MUSCLE ARTIFACTS
% Any one of these functions can be commented out to ignore those artifacts
% when creating the mask    
if config.Do_MWF_Once==1

    % Use epoched data and FFT to detect slope of log frequency log
    % power, add periods exceeding muscle threshold to mask:
    [continuousEEG, epochedEEG] = RELAX_muscle(continuousEEG, epochedEEG, config);  
    if config.computerawmetrics==1
        [continuousEEG, epochedEEG] = RELAX_metrics_muscle(continuousEEG, epochedEEG, config); % record muscle contamination metrics from raw data for comparison.
    end

    EEG=continuousEEG; % Return continuousEEG to the "EEG" variable for MWF processing

    % If including eye blink cleaning in first round MWF, then insert
    % eye blink mask into noise mask:
    if config.MWFRoundToCleanBlinks==1
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
    [EEG] = RELAX_pad_brief_mask_periods (EEG, config, 'notblinks'); % If period has been marked as shorter than config.MinimumArtifactDuration, then pad it out.
    
    EEG.RELAX.NoiseMaskFullLengthR1=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
    EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal=mean(EEG.RELAXProcessing.Details.NoiseMaskFullLength,'omitnan');
    EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR1=EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal; 

    %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
    config.MWFDelayPeriod=config.MWFDelayPeriod_for_muscle_artifacts; 
    config.MWF_delay_spacing=config.MWF_delay_spacing_for_muscle_artifacts;
    [EEG] = RELAX_perform_MWF_cleaning (EEG, config);          

    if config.doDebug
        intermediateEEGs{end+1} = EEG;
        intermediateLabels{end+1} = 'After MWF Round 1';
    end

    EEG.RELAXProcessingRoundOne=EEG.RELAXProcessing; % Record MWF cleaning details from round 1 in EEG file          
    RELAXProcessingRoundOne=EEG.RELAXProcessingRoundOne; % Record MWF cleaning details from round 1 into file for all participants
    
    if isfield(RELAXProcessingRoundOne,'Details')
        RELAXProcessingRoundOne=rmfield(RELAXProcessingRoundOne,'Details');
    end
    if config.KeepAllInfo==0
        if isfield(EEG.RELAXProcessingRoundOne,'Details')
            EEG.RELAXProcessingRoundOne=rmfield(EEG.RELAXProcessingRoundOne,'Details');
        end
    end
    
    % Record processing statistics for all participants in single table:
    RELAXProcessingMWFStepOneAllParticipants(FileNumber,:) = struct2table(RELAXProcessingRoundOne,'AsArray',true);
    EEG = rmfield(EEG,'RELAXProcessing');
    % Save round 1 MWF pre-processing:
    if config.saveround1==1
        if ~exist([config.foldername, filesep 'RELAXProcessed' filesep '1xMWF'], 'dir')
            mkdir([config.foldername, filesep 'RELAXProcessed' filesep '1xMWF'])
        end
        SaveSetMWF1 =[config.foldername, filesep 'RELAXProcessed' filesep '1xMWF', filesep config.filename '_MWF1.set'];    
        EEG = pop_saveset( EEG, SaveSetMWF1 ); 
    end
end

%% PERFORM A SECOND ROUND OF MWF. THIS IS HELPFUL IF THE FIRST ROUND DOESN'T SUFFICIENTLY CLEAN ARTIFACTS. 

% This has been suggested to be useful by Somers et al (2018)
% (particularly when used in a cascading fashion). 

% However, I can see risks. If artifact masks fall on task relevant
% activity in both rounds of the MWF, it may be that the task relevant data
% is just cleaned right out of the signal.

if config.Do_MWF_Twice==1

    % v2.0.1 NWB added the following so that if muscle artifacts aren't
    % cleaned by MWF, the second step doesn't go back to the raw data:
    if config.Do_MWF_Once==0
        EEG=continuousEEG;
    end

    EEG.RELAXProcessing.aFileName=cellstr(config.filename);
    EEG.RELAXProcessing.ProportionMarkedBlinks=0;
    
    % If blinks weren't initially detected because they were 
    % disguised by the the muscle artifact, detect them here
    % (this happens in <1/200 cases, but is a good back up).
    if config.ProbabilityDataHasNoBlinks==0
        if EEG.RELAX.IQRmethodDetectedBlinks(1,1)==0
            continuousEEG=EEG;
            [continuousEEG, epochedEEG] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, config);
            EEG=continuousEEG;
        end
    end
    
    % If including eye blink cleaning in second round MWF, then insert
    % eye blink mask into noise mask:
    if isfield(EEG.RELAX, 'eyeblinkmask')
        if config.MWFRoundToCleanBlinks==2
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
    [EEG] = RELAX_pad_brief_mask_periods (EEG, config, 'blinks');
    
    EEG.RELAX.NoiseMaskFullLengthR2=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
    EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal=mean(EEG.RELAXProcessing.Details.NoiseMaskFullLength,'omitnan');
    EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR2=EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal; 

    %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
    config.MWFDelayPeriod=config.MWFDelayPeriod_for_eye_movements; 
    config.MWF_delay_spacing=config.MWF_delay_spacing_for_eye_movements; % set how sparsely the delay stacking is spread
    [EEG] = RELAX_perform_MWF_cleaning (EEG, config);  

    if config.doDebug
        intermediateEEGs{end+1} = EEG;
        intermediateLabels{end+1} = 'After MWF Round 2';
    end

    EEG.RELAXProcessingRoundTwo=EEG.RELAXProcessing; % Record MWF cleaning details from round 2 in EEG file
    RELAXProcessingRoundTwo=EEG.RELAXProcessingRoundTwo; % Record MWF cleaning details from round 2 into file for all participants
    if isfield(RELAXProcessingRoundTwo,'Details')
        RELAXProcessingRoundTwo=rmfield(RELAXProcessingRoundTwo,'Details');
    end
    if config.KeepAllInfo==0
        if isfield(EEG.RELAXProcessingRoundTwo,'Details')
            EEG.RELAXProcessingRoundTwo=rmfield(EEG.RELAXProcessingRoundTwo,'Details');
        end
    end
    % Record processing statistics for all participants in single table:
    RELAXProcessingMWFStepTwoAllParticipants(FileNumber,:) = struct2table(RELAXProcessingRoundTwo,'AsArray',true);
    EEG = rmfield(EEG,'RELAXProcessing');
    % Save round 2 MWF pre-processing:
    if config.saveround2==1
        if ~exist([config.foldername, filesep 'RELAXProcessed' filesep '2xMWF'], 'dir')
            mkdir([config.foldername, filesep 'RELAXProcessed' filesep '2xMWF'])
        end
        SaveSetMWF2 =[config.foldername, filesep 'RELAXProcessed' filesep '2xMWF', filesep config.filename '_MWF2.set'];    
        EEG = pop_saveset( EEG, SaveSetMWF2 ); 
    end     
end

%% PERFORM A THIRD ROUND OF MWF.    
if config.Do_MWF_Thrice==1

    % v2.0.1 NWB added the following so that if muscle artifacts & blinks 
    % aren't cleaned by MWF, the 3rd step doesn't go back to the raw data:
    if config.Do_MWF_Once==0 && config.Do_MWF_Twice==0
        EEG=continuousEEG;
    end

    EEG.RELAXProcessing.aFileName=cellstr(config.filename);
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
    [continuousEEG, epochedEEG] = RELAX_epoching(EEG, config);
    
    %% THIS SECTION CONTAINS FUNCTIONS WHICH MARK ARTIFACTS

    [continuousEEG, epochedEEG] = RELAX_drift(continuousEEG, epochedEEG, config); % Use epoched data to add periods showing excessive drift to the mask
    
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
    [continuousEEG] = RELAX_horizontaleye(continuousEEG, config);

    %% Return to the "EEG" variable for MWF processing:
    EEG=continuousEEG;
    
    % If including eye blink cleaning in third round MWF, then insert
    % eye blink mask into noise mask:
    if config.MWFRoundToCleanBlinks==3
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
    [EEG] = RELAX_pad_brief_mask_periods (EEG, config, 'notblinks');
    
    EEG.RELAX.NoiseMaskFullLengthR3=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
    EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal=mean(EEG.RELAXProcessing.Details.NoiseMaskFullLength,'omitnan');
    EEG.RELAX.ProportionMarkedInMWFArtifactMaskTotalR3=EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal; 

    %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
    config.MWFDelayPeriod=config.MWFDelayPeriod_for_eye_movements; 
    config.MWF_delay_spacing=config.MWF_delay_spacing_for_eye_movements; % set how sparsely the delay stacking is spread
    [EEG] = RELAX_perform_MWF_cleaning (EEG, config);               
    
    if config.doDebug
        intermediateEEGs{end+1} = EEG;
        intermediateLabels{end+1} = 'After MWF Round 3';
    end

    if isfield(EEG.RELAX, 'eyeblinkmask') % if eyeblinkmask has been created, do the following (thanks to Jane Tan for the suggested bug fix when eyeblinkmask is not created)
        EEG.RELAX=rmfield(EEG.RELAX,'eyeblinkmask'); % remove variables that are no longer necessary
    end
    
    EEG.RELAXProcessingRoundThree=EEG.RELAXProcessing; % Record MWF cleaning details from round 3 in EEG file
    RELAXProcessingRoundThree=EEG.RELAXProcessing; % Record MWF cleaning details from round 3 into file for all participants
    
    if isfield(RELAXProcessingRoundThree,'Details')
        RELAXProcessingRoundThree=rmfield(RELAXProcessingRoundThree,'Details');
    end
    if config.KeepAllInfo==0
        if isfield(EEG.RELAXProcessingRoundThree,'Details')
            EEG.RELAXProcessingRoundThree=rmfield(EEG.RELAXProcessingRoundThree,'Details');
        end
    end
    % Record processing statistics for all participants in single table:
    RELAXProcessingMWFStepThreeAllParticipants(FileNumber,:) = struct2table(RELAXProcessingRoundThree,'AsArray',true);
    EEG = rmfield(EEG,'RELAXProcessing');

    if config.saveround3==1
        if ~exist([config.foldername, filesep 'RELAXProcessed' filesep '3xMWF'], 'dir')
            mkdir([config.foldername, filesep 'RELAXProcessed' filesep '3xMWF'])
        end
        SaveSetMWF3 =[config.foldername,filesep 'RELAXProcessed' filesep '3xMWF', filesep config.filename '_MWF3.set'];    
        EEG = pop_saveset( EEG, SaveSetMWF3 ); 
    end         
end

%% Perform robust average re-referencing of the data, reject periods marked as extreme outliers    
if config.Do_MWF_Once==0 && config.Do_MWF_Twice==0 && config.Do_MWF_Thrice==0 % v2.0.1 NWB adjusted to only return to continuous EEG if no MWF was applied (rather than just that the first MWF wasn't applied)
    EEG=continuousEEG;
end

% v2.0.1 NWB 8/8/2025 - moved the low pass filtering after the MWF
% steps to before the extreme bad period rejection steps to ensure
% filtering is not adversely affected by boundary elements.
if strcmp(config.LowPassFilterBeforeMWF,'no') % if low pass filtering wasn't applied before MWF cleaning (recommended) apply it here
    if strcmp(config.FilterType,'Butterworth')
        EEG = RELAX_filtbutter( EEG, [], config.LowPassFilter, 4, 'lowpass', config.causal_or_acausal_filter);
    end
    if strcmp(config.FilterType,'pop_eegfiltnew')
        EEG = pop_eegfiltnew(EEG,[],config.LowPassFilter);
    end
end
      
% Reject periods that were marked as NaNs in the MWF masks because they 
% showed extreme shift within the epoch or extremely improbable data:
EEG = eeg_eegrej( EEG, EEG.RELAX.ExtremelyBadPeriodsForDeletion);

if config.doDebug
    intermediateEEGs{end+1} = EEG;
    intermediateLabels{end+1} = 'After extreme period rejection';
end

if ~isempty(config.electrodes_2_keep_but_not_clean{1})
    % reject the same periods from the auxilary electrodes that weren't included in cleaning:
    Non_cleaned_electrodes = eeg_eegrej( Non_cleaned_electrodes, EEG.RELAX.ExtremelyBadPeriodsForDeletion);
    % store the auxilary electrodes in the EEG file:
    EEG.Non_cleaned_electrodes=Non_cleaned_electrodes;
end

[EEG] = RELAX_average_rereference(EEG);
EEG = eeg_checkset( EEG );  

if config.doDebug
    intermediateEEGs{end+1} = EEG;
    intermediateLabels{end+1} = 'After average rereferencing';
end

%% Perform wICA on ICLabel identified artifacts that remain:
if config.Perform_targeted_wICA==1
    % The following cleans eye movements and muscle artifacts in the
    % independent component space by a combination of wavelet enhanced
    % ICA cleaning and targeting to restrict the cleaning to only 
    % artifact periods for eye movement components, and high pass 
    % filtering muscle components at 15Hz:
    EEG.RELAXProcessing_wICA.aFileName=cellstr(config.filename);
    [EEG] = RELAX_targeted_wICA(EEG,config);
    % setting 'config.Report_all_wICA_info' to 1 will report proportion of ICs categorized as each category, and variance explained by ICs from each category (function is ~20s slower if this is implemented)
    EEG = eeg_checkset( EEG );
    RELAXProcessing_wICA=EEG.RELAXProcessing_wICA;
    % Record processing statistics for all participants in single table:
    RELAXProcessing_wICA_AllParticipants = struct2table(RELAXProcessing_wICA,'AsArray',true);
    
    if config.doDebug
        intermediateEEGs{end+1} = EEG;
        intermediateLabels{end+1} = 'After targeted wICA';
    end
end


%% Perform wICA on ICLabel identified artifacts that remain:
if config.Perform_wICA_on_ICLabel==1
    % The following performs wICA, implemented on only the components
    % marked as artifact by ICLabel.
    EEG.RELAXProcessing_wICA.aFileName=cellstr(config.filename);
    [EEG,~, ~, ~, ~] = RELAX_wICA_on_ICLabel_artifacts(EEG,config.ICA_method, 1, 0, EEG.srate, 5,'coif5',config.Report_all_ICA_info,config.ICLabel_thresholds,config.Clean_other_comps); 
    % setting 'config.Report_all_wICA_info' to 1 will report proportion of ICs categorized as each category, and variance explained by ICs from each category (function is ~20s slower if this is implemented)
    EEG = eeg_checkset( EEG );
    RELAXProcessing_wICA=EEG.RELAXProcessing_wICA;
    % Record processing statistics for all participants in single table:
    RELAXProcessing_wICA_AllParticipants(FileNumber,:) = struct2table(RELAXProcessing_wICA,'AsArray',true);
    
    if config.doDebug
        intermediateEEGs{end+1} = EEG;
        intermediateLabels{end+1} = 'After wICA on ICLabel';
    end
end

%% Perform ICA subtract on ICLabel identified artifacts that remain:
if config.Perform_ICA_subtract==1
    % The following performs ICA sutraction, implemented on only the components
    % marked as artifact by ICLabel.
    EEG.RELAXProcessing_ICA.aFileName=cellstr(config.filename);
    EEG = RELAX_ICA_subtract(EEG,config);
    EEG = eeg_checkset( EEG );
    RELAXProcessing_ICA=EEG.RELAXProcessing_ICA;
    % Record processing statistics for all participants in single table:
    RELAXProcessing_ICA_AllParticipants(FileNumber,:) = struct2table(RELAXProcessing_ICA,'AsArray',true);
    
    if config.doDebug
        intermediateEEGs{end+1} = EEG;
        intermediateLabels{end+1} = 'After ICA subtract';
    end
end

EEG.RELAX.Data_has_been_cleaned=1;

%% Save ICA topoplots
if exist('save_ICA_topo.m', 'file')
    save_ICA_topo(EEG, config);  % Pass EEG structure and config
end

%% COMPUTE CLEANED METRICS:
if config.computecleanedmetrics==1    
    [continuousEEG, epochedEEG] = RELAX_epoching(EEG, config);
    [continuousEEG, ~] = RELAX_metrics_blinks(continuousEEG, epochedEEG);
    [continuousEEG, ~] = RELAX_metrics_muscle(continuousEEG, epochedEEG, config);

    [continuousEEG] = RELAX_metrics_final_SER_and_ARR(rawEEG, continuousEEG); % this is only a good metric for testing only the cleaning of artifacts marked for cleaning by MWF, see notes in function.

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
            if isfield(EEG.RELAX_Metrics.Cleaned,'All_SER')
                CleanedMetrics.All_SER(FileNumber)=EEG.RELAX_Metrics.Cleaned.All_SER;
                CleanedMetrics.All_ARR(FileNumber)=EEG.RELAX_Metrics.Cleaned.All_ARR;
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
EEG.RELAX_issues_to_check.aFileName=cellstr(config.filename);
if size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1)>config.MaxProportionOfElectrodesThatCanBeDeleted*size(EEG.allchan,2)
    EEG.RELAX_issues_to_check.PREP_rejected_too_many_electrodes=size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1); % 1.1.4: fix dimension specification error
else
    EEG.RELAX_issues_to_check.PREP_rejected_too_many_electrodes=0;
end
if (EEG.RELAXProcessingExtremeRejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete...
        +EEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete...
        +size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1))...
        >=config.MaxProportionOfElectrodesThatCanBeDeleted*size(EEG.allchan,2)
    EEG.RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold=...
        (EEG.RELAXProcessingExtremeRejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete...
        +EEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete...
        +size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1));
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
if config.Do_MWF_Once==1
    EEG.RELAX_issues_to_check.MWF_eigenvector_deficiency_R1=isa(EEG.RELAXProcessingRoundOne.RankDeficiency,'char');
end
if config.Do_MWF_Twice==1
    EEG.RELAX_issues_to_check.MWF_eigenvector_deficiency_R2=isa(EEG.RELAXProcessingRoundTwo.RankDeficiency,'char');
end
if config.Do_MWF_Thrice==1
    EEG.RELAX_issues_to_check.MWF_eigenvector_deficiency_R3=isa(EEG.RELAXProcessingRoundThree.RankDeficiency,'char');
end
if config.Perform_wICA_on_ICLabel==1
    if EEG.RELAXProcessing_wICA.Proportion_artifactICs_reduced_by_wICA>0.80
        EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=EEG.RELAXProcessing_wICA.Proportion_artifactICs_reduced_by_wICA;
    else
        EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=0;
    end
    EEG.RELAX_issues_to_check.DataMaybeTooShortForValidICA = EEG.RELAXProcessing_wICA.DataMaybeTooShortForValidICA;
    EEG.RELAX_issues_to_check.fastica_symm_Didnt_Converge=EEG.RELAXProcessing_wICA.fastica_symm_Didnt_Converge(1,3);
end
if config.Perform_ICA_subtract==1
    if EEG.RELAXProcessing_ICA.Proportion_artifactICs_reduced_by_ICA>0.80
        EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=EEG.RELAXProcessing_ICA.Proportion_artifactICs_reduced_by_ICA;
    else
        EEG.RELAX_issues_to_check.HighProportionOfArtifact_ICs=0;
    end
    EEG.RELAX_issues_to_check.DataMaybeTooShortForValidICA = EEG.RELAXProcessing_ICA.DataMaybeTooShortForValidICA;
    EEG.RELAX_issues_to_check.fastica_symm_Didnt_Converge=EEG.RELAXProcessing_ICA.fastica_symm_Didnt_Converge(1,3);
end

if strcmp(config.InterpolateRejectedElectrodesAfterCleaning,'yes')
    EEG = pop_interp(EEG, EEG.allchan, 'spherical');
    
    if config.doDebug
        intermediateEEGs{end+1} = EEG;
        intermediateLabels{end+1} = 'After electrode interpolation';
    end
end

%% Store debug information in EEG structure
if config.doDebug
    c_saySingle('Debug mode: Stored %d intermediate processing stages', length(intermediateLabels));
end

%% SAVE FILE:
if ~exist([config.foldername, filesep 'RELAXProcessed' filesep 'Cleaned_Data'], 'dir')
    mkdir([config.foldername, filesep 'RELAXProcessed' filesep 'Cleaned_Data'])
end
SaveSet_CleanedFile =[config.foldername,filesep 'RELAXProcessed' filesep 'Cleaned_Data', filesep config.filename '_RELAX.set'];  
EEG.RELAX_settings_used_to_clean_this_file=config;
EEG = pop_saveset( EEG, SaveSet_CleanedFile ); 

% Record warnings for all participants in single table:
try
    RELAX_issues_to_check = struct2table(EEG.RELAX_issues_to_check,'AsArray',true);
catch
    RELAX_issues_to_check_2nd_run = struct2table(EEG.RELAX_issues_to_check,'AsArray',true);
    warning('The variable: "RELAX_issues_to_check" already exists and includes different settings from your current settings')
    warning('This is likely from a previous run of RELAX. Saving variable as "RELAX_issues_to_check_2nd_run" instead');
end

%% Save statistics for each participant and across participants, graph cleaning metrics:

% Also set empty output variables in case these are not produced because certain
% parameters have been switched off:

savefileone=[config.myPath filesep 'RELAXProcessed' filesep 'RELAXProcessingExtremeRejectionsAllParticipants'];
save(savefileone,'RELAXProcessingExtremeRejectionsAllParticipants')
if config.Do_MWF_Once==1
    savefileone=[config.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatisticsRoundOne'];
    save(savefileone,'RELAXProcessingMWFStepOneAllParticipants')
else
    RELAXProcessingMWFStepOneAllParticipants={};
end
if config.Do_MWF_Twice==1
    savefiletwo=[config.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatisticsRoundTwo'];
    save(savefiletwo,'RELAXProcessingMWFStepTwoAllParticipants')
else
    RELAXProcessingMWFStepTwoAllParticipants={};
end
if config.Do_MWF_Thrice==1
    savefilethree=[config.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatisticsRoundThree'];
    save(savefilethree,'RELAXProcessingMWFStepThreeAllParticipants')
else
    RELAXProcessingMWFStepThreeAllParticipants={};
end
if config.Perform_wICA_on_ICLabel==1 || config.Perform_targeted_wICA==1
    savefilefour=[config.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatistics_wICA'];
    save(savefilefour,'RELAXProcessing_wICA_AllParticipants')
else
    RELAXProcessing_wICA_AllParticipants={}; 
end
if config.Perform_ICA_subtract==1
    savefilefour=[config.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatistics_ICA'];
    save(savefilefour,'RELAXProcessing_ICA_AllParticipants')
else
    RELAXProcessing_ICA_AllParticipants={}; 
end
if exist('CleanedMetrics','var')
    savemetrics=[config.myPath filesep 'RELAXProcessed' filesep 'CleanedMetrics'];
    save(savemetrics,'CleanedMetrics')
else
    CleanedMetrics={};
end
if exist('RawMetrics','var')
    savemetrics=[config.myPath filesep 'RELAXProcessed' filesep 'RawMetrics'];
    save(savemetrics,'RawMetrics')
else
    RawMetrics={};
end
if exist('RELAX_issues_to_check','var')
    savemetrics=[config.myPath filesep 'RELAXProcessed' filesep 'RELAX_issues_to_check'];
    save(savemetrics,'RELAX_issues_to_check')
end
if exist('RELAX_issues_to_check_2nd_run','var')
    savemetrics=[config.myPath filesep 'RELAXProcessed' filesep 'RELAX_issues_to_check_2nd_run'];
    save(savemetrics,'RELAX_issues_to_check_2nd_run')
end
config.filename=[];
savefileone=[config.myPath filesep 'RELAXProcessed' filesep 'config'];
save(savefileone,'config')

set(groot, 'defaultAxesTickLabelInterpreter','none');
if config.computecleanedmetrics==1
    % Create output directory for figures
    fig_dir = fullfile(config.myPath, 'RELAXProcessed', 'QualityMetrics_Figures');
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
    % Plot QC Metrics
    try
        fig1 = figure('Name','BlinkAmplitudeRatio','units','normalized','outerposition',[0.05 0.05 0.95 0.95]);
        boxplot(CleanedMetrics.BlinkAmplitudeRatio);
        xticklabels(config.files); xtickangle(90);
        set(gca,'FontSize',16, 'FontWeight', 'bold') % Creates an axes and sets its FontSize to 21
        % Save figure
        saveas(fig1, fullfile(fig_dir, 'BlinkAmplitudeRatio.png'));
        saveas(fig1, fullfile(fig_dir, 'BlinkAmplitudeRatio.fig'));
    catch
    end
    try
        fig2 = figure('Name','MeanMuscleStrengthFromOnlySuperThresholdValues','units','normalized','outerposition',[0.05 0.05 0.95 0.95]);
        b=bar(CleanedMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues); 
        xtickangle(90); xticks([1:1:size(config.files,2)]); b(1).BaseValue = config.MuscleSlopeThreshold;
        xticklabels(config.files); ylim([config.MuscleSlopeThreshold max(CleanedMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues)+1]);b.ShowBaseLine='off';
        set(gca,'FontSize',16, 'FontWeight', 'bold') % Creates an axes and sets its FontSize to 21
        % Save figure
        saveas(fig2, fullfile(fig_dir, 'MeanMuscleStrengthFromOnlySuperThresholdValues.png'));
        saveas(fig2, fullfile(fig_dir, 'MeanMuscleStrengthFromOnlySuperThresholdValues.fig'));
    catch
    end
    try
        fig3 = figure('Name','ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel','units','normalized','outerposition',[0.05 0.05 0.95 0.95]);
        bar(CleanedMetrics.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel);
        set(gca,'FontSize',16, 'FontWeight', 'bold') % Creates an axes and sets its FontSize to 21
        xtickangle(90); xticks([1:1:size(config.files,2)]);
        xticklabels(config.files);
        % Save figure
        saveas(fig3, fullfile(fig_dir, 'ProportionEpochsWithMuscleAboveThresholdAnyChannel.png'));
        saveas(fig3, fullfile(fig_dir, 'ProportionEpochsWithMuscleAboveThresholdAnyChannel.fig'));
    catch
    end
end

clearvars -except 'config' 'FileNumber' 'CleanedMetrics' 'RawMetrics' 'RELAXProcessingMWFStepOneAllParticipants' 'RELAXProcessingMWFStepTwoAllParticipants' 'RELAXProcessing_wICA_AllParticipants'...
    'RELAXProcessing_ICA_AllParticipants' 'RELAXProcessingMWFStepThreeAllParticipants' 'Warning' 'RELAX_issues_to_check' 'RELAX_issues_to_check_2nd_run'...
    'RELAXProcessingExtremeRejectionsAllParticipants' 'WarningAboutFileNumber' 'EEG';

if ~exist('RELAX_issues_to_check_2nd_run','var')    
warning('Check "RELAX_issues_to_check" to see if any issues were noted for specific files');
RELAX_issues_to_check_2nd_run=['This variable is only here as a placeholder in case you have already run RELAX once,...' ...
    ' and had previously left the output variables in the same folder as the folder where you have saved the currently cleaned data.' ...
    'This was not an issue with the current run, so this placeholder variable was not filled'];
elseif exist('RELAX_issues_to_check_2nd_run','var')
warning('Check "RELAX_issues_to_check_2nd_run" to see if any issues were noted for specific files');
end

if config.ProbabilityDataHasNoBlinks<2 && sum(RELAX_issues_to_check.NoBlinksDetected)>1
f = msgbox('RELAX did not detect any blinks for some files. Open the "RELAX_issues_to_check" struct in the workspace to check which files. We recommend visually inspecting these files to ensure there has not been an error.'...
,'No blinks detected for some files');    
set(f,'Position',[500,500,450,100]);
ah = get( f, 'CurrentAxes' );
ch = get( ah, 'Children' );
set( ch, 'FontSize', 12 ); %makes text bigger
end

if find(RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold>0)>0
f = msgbox('Some files met or exceeded the electrode rejection threshold. We recommend visually inspecting the raw and cleaned files where this is the case. Open the "RELAX_issues_to_check" struct in the workspace, and check the third column. Files that exceeded the threshold will show a value above 0. Exclude files where raw data seems irretrievably noisy, or cleaned data still contains excessive noise.'...
,'Some files met or exceeded the electrode rejection threshold');    
set(f,'Position',[300,300,450,150]);
ah = get( f, 'CurrentAxes' );
ch = get( ah, 'Children' );
set( ch, 'FontSize', 12 ); %makes text bigger
end

c_saySingle('AARATEP-RELAX pipeline preprocessing complete');

end

% Helper function for conditional statements
function result = c_if(condition, trueVal, falseVal)
    if condition
        result = trueVal;
    else
        result = falseVal;
    end
end

% Helper function for printing messages
function c_say(varargin)
    fprintf([varargin{1} '...\n'], varargin{2:end});
end

function c_sayDone(varargin)
    if nargin > 0
        fprintf([varargin{1} '\n'], varargin{2:end});
    else
        fprintf('Done.\n');
    end
end

function c_saySingle(varargin)
    fprintf([varargin{1} '\n'], varargin{2:end});
end
