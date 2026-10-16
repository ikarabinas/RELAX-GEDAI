% RELAX_QC_load_tables() - Join the wrapper's per-participant aggregate outputs into one
%                          table keyed by file base name, for the QC figure set.
%
% Usage:
%   >> [T, cfgs] = RELAX_QC_load_tables(statsDir, cleanedDir);
%
% Inputs:
%   statsDir   - RELAXProcessed/PerParticipantStats
%   cleanedDir - RELAXProcessed/Cleaned_Data
%
% Outputs:
%   T    - one row per file named in any RELAX_cfg_<participant>.mat. RawMetrics and
%          CleanedMetrics are indexed by FileNumber, i.e. by position in that
%          participant's RELAX_cfg.files; the covariate tables are keyed by aFileName.
%          Outcome is 'cleaned' if <File>_RELAX.set exists, else the Outcome recorded in
%          RELAX_unprocessed_files_<participant>.mat, else 'notProcessed'.
%   cfgs - containers.Map, participant tag -> that participant's RELAX_cfg
%
% A missing per-participant file warns (RELAX_QC:missingStats) and leaves NaN/empty
% values; an absent covariate column gives NaN, never 0.

function [T, cfgs] = RELAX_QC_load_tables(statsDir, cleanedDir)

    cfgFiles = dir(fullfile(statsDir, 'RELAX_cfg_*.mat'));
    cfgs = containers.Map();
    perParticipant = cell(numel(cfgFiles), 1);

    for k = 1:numel(cfgFiles)
        participant = regexprep(cfgFiles(k).name, '^RELAX_cfg_(.*)\.mat$', '$1');
        cfg         = loadStatsVariable(statsDir, 'RELAX_cfg', participant, 'RELAX_cfg', false);
        cfgs(participant) = cfg;
        raw         = loadStatsVariable(statsDir, 'RawMetrics', participant, 'RawMetrics', false);
        cleaned     = loadStatsVariable(statsDir, 'CleanedMetrics', participant, 'CleanedMetrics', false);
        extreme     = loadStatsVariable(statsDir, 'RELAXProcessingExtremeRejectionsAllParticipants', participant, ...
                                        'RELAXProcessingExtremeRejectionsAllParticipants', false);
        wica        = loadStatsVariable(statsDir, 'ProcessingStatistics_wICA', participant, ...
                                        'RELAXProcessing_wICA_AllParticipants', false);
        unprocessed = loadStatsVariable(statsDir, 'RELAX_unprocessed_files', participant, ...
                                        'RELAX_unprocessed_files', true);

        files = string(cellfun(@baseName, cfg.files(:), 'UniformOutput', false));
        n = numel(files);
        R = table(files, repmat(string(participant), n, 1), (1:n)', repmat("notProcessed", n, 1), repmat("", n, 1), ...
                  'VariableNames', {'File', 'Participant', 'FileNumber', 'Outcome', 'Message'});

        for i = 1:n
            if isfile(fullfile(cleanedDir, files(i) + "_RELAX.set"))
                R.Outcome(i) = "cleaned";
            elseif istable(unprocessed)
                hit = find(string(cellfun(@baseName, unprocessed.File, 'UniformOutput', false)) == files(i), 1, 'last');
                if ~isempty(hit)
                    R.Outcome(i) = string(unprocessed.Outcome{hit});
                    R.Message(i) = string(unprocessed.Message{hit});
                end
            end
        end

        srate = cfg.DownSample_to_X_Hz;
        covariates = { ...
            'NChanImpedance',  'ImpedanceBasedChannelToReject',  @countEntries; ...
            'NChanPop',        'PopBasedChannelToReject',        @countEntries; ...
            'NChanBridge',     'BridgedChannelToReject',         @countEntries; ...
            'NChanPREP',       'PREPBasedChannelToReject',       @countEntries; ...
            'NChanMuscle',     'MuscleBasedElectrodesToReject',  @countEntries; ...
            'NChanExtreme',    'ExtremeDataBasedChannelToReject', @countEntries; ...
            'ExtremeExcisedSeconds', 'ExtremelyBadPeriodsForDeletion', @(v) periodSamples(v) / srate; ...
            'PopExcisedSeconds',              'PopExcisedSeconds',              @scalarOf; ...
            'PopFractionSamplesInterpolated', 'PopFractionSamplesInterpolated', @scalarOf; ...
            'PopNSteps',                      'PopNSteps',                      @scalarOf; ...
            'BlinkElectrodeFallback',         'BlinkElectrodeFallback',         @scalarOf};
        for c = 1:size(covariates, 1)
            R.(covariates{c, 1}) = arrayfun(@(f) covariate(keyedRow(extreme, f), covariates{c, 2}, covariates{c, 3}), files);
        end
        R.WicaArtifactICsReduced = arrayfun(@(f) covariate(keyedRow(wica, f), ...
            'Proportion_artifactICs_reduced_by_wICA', @scalarOf), files);

        R.TableRawBlinkRatio           = metricColumns(raw, 'BlinkAmplitudeRatio', n);
        R.TableCleanedBlinkRatio       = metricColumns(cleaned, 'BlinkAmplitudeRatio', n);
        R.TableCleanedMuscleProportion = metricScalars(cleaned, 'ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel', n);
        R.TableCleanedMuscleStrength   = metricScalars(cleaned, 'MeanMuscleStrengthFromOnlySuperThresholdValues', n);
        % The wrapper grows these arrays by FileNumber, so a file that never reached the
        % metrics stage reads back as 0, not as missing:
        notCleaned = R.Outcome ~= "cleaned";
        R.TableCleanedMuscleProportion(notCleaned) = NaN;
        R.TableCleanedMuscleStrength(notCleaned)   = NaN;
        R.TableCleanedBlinkRatio(notCleaned)       = {NaN};

        perParticipant{k} = R;
    end
    T = vertcat(perParticipant{:});
end

function name = baseName(f)
    [~, name] = fileparts(char(f));
end

function v = loadStatsVariable(statsDir, prefix, participant, varName, quiet)
    v = [];
    f = fullfile(statsDir, sprintf('%s_%s.mat', prefix, participant));
    if ~isfile(f)
        if ~quiet, warning('RELAX_QC:missingStats', 'Missing %s', f); end
        return
    end
    s = load(f, varName);
    if isfield(s, varName)
        v = s.(varName);
    else
        warning('RELAX_QC:missingStats', '%s has no variable %s', f, varName);
    end
end

function row = keyedRow(t, file)
    row = [];
    if ~istable(t) || ~ismember('aFileName', t.Properties.VariableNames), return, end
    hit = find(string(t.aFileName) == file, 1, 'last');
    if ~isempty(hit), row = t(hit, :); end
end

function x = covariate(row, name, reduce)
    x = NaN;
    if isempty(row) || ~ismember(name, row.Properties.VariableNames), return, end
    v = row.(name);
    while iscell(v) && isscalar(v), v = v{1}; end
    x = reduce(v);
end

function n = countEntries(v)
    if ischar(v) || isstring(v)
        n = double(any(strlength(string(v)) > 0));
    else
        n = numel(v);
    end
end

function n = periodSamples(v)
    n = 0;
    if isnumeric(v) && size(v, 2) == 2, n = sum(v(:, 2) - v(:, 1) + 1); end
end

function x = scalarOf(v)
    x = NaN;
    if isnumeric(v) || islogical(v)
        if ~isempty(v), x = double(v(1)); end
    end
end

function c = metricColumns(s, field, n)
    c = repmat({NaN}, n, 1);
    if ~isstruct(s) || ~isfield(s, field), return, end
    m = s.(field);
    for i = 1:min(n, size(m, 2)), c{i} = m(:, i); end
end

function x = metricScalars(s, field, n)
    x = nan(n, 1);
    if ~isstruct(s) || ~isfield(s, field), return, end
    m = s.(field);
    k = min(n, numel(m));
    x(1:k) = m(1:k);
end
