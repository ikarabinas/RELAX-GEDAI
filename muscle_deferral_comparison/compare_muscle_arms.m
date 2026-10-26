% compare_muscle_arms() - Compare the cleaned files of the two muscle-screen placements.
%
% Usage:
%   >> compare_muscle_arms
%
% Reads <myPath>/RELAXProcessed_muscleBefore and _muscleAfter (written by run_arm.sh) and,
% for every file cleaned by both arms, measures the final saved data the same way in both:
% per-channel muscle slope (RELAX_metrics_muscle, 1 s epochs, 7-75 Hz log-log slope) and
% mean log power by band. Channels are grouped so both arms are compared on the same set:
%   baselineMuscle - electrodes the Before arm rejected for muscle (interpolated there; in
%                    the After arm cleaned in place unless the deferred screen rejected them)
%   keptInBoth     - electrodes neither arm interpolated
% Writes results/muscle_arms_by_file.csv and results/muscle_arms.mat next to this script,
% and runs the residual screen's score stage on each arm into its own Residual_screen/.

function compare_muscle_arms()

    myPath = '/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Oct26/rest_treat_joined/';
    arms = ["Before" "After"];
    muscleCfg = struct('MuscleSlopeThreshold', -0.31, 'LowPassFilter', 200); % as in RELAX_SET_PARAMETERS_AND_RUN_slurm.m
    bands = struct('name', {'alpha', 'beta', 'gamma', 'highGamma'}, 'hz', {[8 13], [15 30], [30 50], [65 75]});
    holeDepthLimit = 1.4; % as in RELAX_residual_screen.m
    resultsDir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~isfolder(resultsDir), mkdir(resultsDir); end

    cleanedDirs = arrayfun(@(a) fullfile(myPath, "RELAXProcessed_muscle" + a, "Cleaned_Data"), arms);
    names = intersect(cleanedNames(cleanedDirs(1)), cleanedNames(cleanedDirs(2)));
    fprintf('compare_muscle_arms: %d files cleaned by both arms\n', numel(names));

    rows = cell(numel(names), 1);
    spectra = struct('file', {}, 'arm', {}, 'freqs', {}, 'baselineMuscle', {}, 'keptInBoth', {});
    for f = 1:numel(names)
        measured = struct();
        for a = 1:2
            EEG = pop_loadset('filename', char(names(f)), 'filepath', char(cleanedDirs(a)));
            measured.(arms(a)) = measureArm(EEG, muscleCfg, holeDepthLimit);
        end
        baselineMuscle = measured.Before.muscleRejected;
        keptInBoth = setdiff(measured.Before.labels, [measured.Before.interpolated, measured.After.interpolated], 'stable');

        row = table(string(regexprep(names(f), '_RELAX.*\.set$', '')), numel(baselineMuscle), numel(keptInBoth), ...
            'VariableNames', {'File', 'NBaselineMuscleChannels', 'NKeptInBoth'});
        for a = 1:2
            m = measured.(arms(a));
            p = arms(a) + "_";
            row.(p + "NInterpolated") = numel(m.interpolated);
            row.(p + "NMuscleRejected") = numel(m.muscleRejected);
            row.(p + "MuscleOverThresholdBeforeCleaning") = m.overThresholdBeforeCleaning;
            row.(p + "MuscleOverThresholdAfterCleaning") = m.overThresholdAfterCleaning;
            row.(p + "HoleDepth") = m.holeDepth;
            row.(p + "Seconds") = m.seconds;
            row.(p + "ProportionExcludedExtreme") = m.proportionExcludedExtreme;
            row.(p + "ICsMuscle") = m.icsMuscle;
            row.(p + "VarianceMuscleICs") = m.varianceMuscleICs;
            row.(p + "SENSAI") = m.sensai;
            groups = struct('baselineMuscle', baselineMuscle, 'keptInBoth', keptInBoth, 'all', m.labels);
            for group = string(fieldnames(groups))'
                inGroup = ismember(m.labels, groups.(group));
                row.(p + group + "_MuscleProportion") = mean(m.proportion(inGroup), 'omitnan');
                row.(p + group + "_MuscleExceedance") = mean(m.exceedance(inGroup), 'omitnan');
                for b = bands
                    inBand = m.freqs >= b.hz(1) & m.freqs <= b.hz(2);
                    row.(p + group + "_" + b.name + "LogPower") = mean(m.logPower(inGroup, inBand), 'all', 'omitnan');
                end
            end
            spectra(end + 1) = struct('file', row.File, 'arm', arms(a), 'freqs', m.freqs, ... %#ok<AGROW>
                'baselineMuscle', mean(m.logPower(ismember(m.labels, baselineMuscle), :), 1, 'omitnan'), ...
                'keptInBoth', mean(m.logPower(ismember(m.labels, keptInBoth), :), 1, 'omitnan'));
        end
        rows{f} = row;
        fprintf('compare_muscle_arms %d/%d %s\n', f, numel(names), row.File);
    end
    results = vertcat(rows{:});
    writetable(results, fullfile(resultsDir, 'muscle_arms_by_file.csv'));
    save(fullfile(resultsDir, 'muscle_arms.mat'), 'results', 'spectra');

    for a = arms
        outputFolderName = "RELAXProcessed_muscle" + a;
        RELAX_screen_run_score(struct('myPath', myPath, 'outputFolderName', char(outputFolderName), ...
            'outDir', fullfile(myPath, outputFolderName, 'Residual_screen'), 'cutoff', 3, ...
            'holeDepthLimit', holeDepthLimit, 'sampleSize', 10, 'seed', 20261022));
    end
end

function m = measureArm(EEG, muscleCfg, holeDepthLimit)
    m.labels = string({EEG.chanlocs.labels});
    m.interpolated = setdiff(string({EEG.allchan.labels}), flattened(EEG.RELAX.ListOfChannelsAfterRejections), 'stable');
    rejections = EEG.RELAXProcessingExtremeRejections;
    m.muscleRejected = flattened(rejections.MuscleBasedElectrodesToReject);
    m.overThresholdBeforeCleaning = fieldOr(rejections, 'NumberOfMuscleContaminatedChannelsBeforeCleaning', ...
        rejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete);
    m.overThresholdAfterCleaning = NaN;
    if fieldOr(rejections, 'RejectMuscleChannelsAfterCleaning', 0) == 1
        m.overThresholdAfterCleaning = rejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete;
    end
    hole = RELAX_screen_hole_depth(EEG.allchan, m.interpolated, holeDepthLimit);
    m.holeDepth = hole.depth;
    m.seconds = EEG.pnts / EEG.srate;
    m.proportionExcludedExtreme = fieldOr(rejections, 'ProportionExcludedForExtremeOutlier', NaN);
    wICA = fieldOr(EEG, 'RELAXProcessing_wICA', struct());
    m.icsMuscle = fieldOr(wICA, 'ProportionICs_was_Muscle', NaN);
    m.varianceMuscleICs = fieldOr(wICA, 'ProportionVariance_was_MuscleICs', NaN);
    m.sensai = NaN;
    if isfield(EEG.etc, 'GEDAI'), m.sensai = EEG.etc.GEDAI.SENSAI_score; end

    EEG.RELAX.Data_has_been_cleaned = 1;
    [~, epochedEEG] = RELAX_metrics_muscle(EEG, [], muscleCfg);
    slopes = epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels; % channels x epochs
    m.proportion = mean(slopes > muscleCfg.MuscleSlopeThreshold, 2)';
    m.exceedance = mean(max(slopes - muscleCfg.MuscleSlopeThreshold, 0), 2)';
    power = epochedEEG.RELAXProcessing.Details.Muscle_Slopes; % epochs x channels x 7:75 Hz
    m.freqs = 7:75;
    m.logPower = squeeze(mean(log10(power), 1, 'omitnan'));
end

function value = fieldOr(s, name, default)
    value = default;
    if isfield(s, name) && ~isempty(s.(name)), value = s.(name); end
end

function names = cleanedNames(cleanedDir)
    listing = dir(fullfile(cleanedDir, '*_RELAX*.set'));
    names = string({listing.name});
end

function list = flattened(value)
    if iscell(value)
        parts = cellfun(@flattened, value(:)', 'UniformOutput', false);
        list = [strings(1, 0), parts{:}];
    elseif ischar(value) || isstring(value)
        list = reshape(string(value), 1, []);
    else
        list = strings(1, 0);
    end
    list = list(~ismissing(list) & list ~= "");
end
