function test_RELAX_QC_load_tables()
% Joins per-participant metric structs (indexed by FileNumber) and covariate tables
% (keyed by aFileName) onto file base names, without erroring on missing pieces.

    root = tempname; statsDir = fullfile(root, 'stats'); cleanedDir = fullfile(root, 'cleaned');
    mkdir(statsDir); mkdir(cleanedDir);
    cleanup = onCleanup(@() rmdir(root, 's'));

    % m001: everything present, two cleaned files
    RELAX_cfg = struct('files', {{'m001_day1.set', 'm001_day2.set'}}, 'DownSample_to_X_Hz', 500); %#ok<NASGU>
    save(fullfile(statsDir, 'RELAX_cfg_m001.mat'), 'RELAX_cfg');
    RawMetrics = struct('BlinkAmplitudeRatio', [1 2; 3 4]); %#ok<NASGU>
    save(fullfile(statsDir, 'RawMetrics_m001.mat'), 'RawMetrics');
    CleanedMetrics = struct('BlinkAmplitudeRatio', [1.1 1.2; 1.3 1.4], ...
        'ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel', [0.1 0.2], ...
        'MeanMuscleStrengthFromOnlySuperThresholdValues', [0.05 0.07]); %#ok<NASGU>
    save(fullfile(statsDir, 'CleanedMetrics_m001.mat'), 'CleanedMetrics');
    RELAXProcessingExtremeRejectionsAllParticipants = table( ...
        {'m001_day1'; 'm001_day2'}, {{'E1', 'E2'}; {}}, {cell(1,0); {'E9'}}, ...
        {[1 500; 1001 1250]; zeros(0, 2)}, [10; 20], [1e-4; 2e-4], [0; 0], [0; 1], ...
        'VariableNames', {'aFileName', 'ImpedanceBasedChannelToReject', 'PopBasedChannelToReject', ...
        'ExtremelyBadPeriodsForDeletion', 'PopNSteps', 'PopFractionSamplesInterpolated', ...
        'PopExcisedSeconds', 'BlinkElectrodeFallback'}); %#ok<NASGU>
    save(fullfile(statsDir, 'RELAXProcessingExtremeRejectionsAllParticipants_m001.mat'), ...
        'RELAXProcessingExtremeRejectionsAllParticipants');
    RELAXProcessing_wICA_AllParticipants = table({'m001_day1'; 'm001_day2'}, [0.3; 0.4], ...
        'VariableNames', {'aFileName', 'Proportion_artifactICs_reduced_by_wICA'}); %#ok<NASGU>
    save(fullfile(statsDir, 'ProcessingStatistics_wICA_m001.mat'), 'RELAXProcessing_wICA_AllParticipants');
    touch(fullfile(cleanedDir, 'm001_day1_RELAX.set')); touch(fullfile(cleanedDir, 'm001_day2_RELAX.set'));

    % m002: no RawMetrics, second file excluded, third never reached (array task still running)
    RELAX_cfg = struct('files', {{'m002_day1.set', 'm002_day2.set', 'm002_day3.set'}}, 'DownSample_to_X_Hz', 500); %#ok<NASGU>
    save(fullfile(statsDir, 'RELAX_cfg_m002.mat'), 'RELAX_cfg');
    CleanedMetrics = struct('BlinkAmplitudeRatio', [1.5; 1.6], ...
        'ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel', 0.3, ...
        'MeanMuscleStrengthFromOnlySuperThresholdValues', 0.09); %#ok<NASGU>
    save(fullfile(statsDir, 'CleanedMetrics_m002.mat'), 'CleanedMetrics');
    RELAXProcessingExtremeRejectionsAllParticipants = table({'m002_day1'}, {{'E3'}}, ...
        'VariableNames', {'aFileName', 'ImpedanceBasedChannelToReject'}); %#ok<NASGU>
    save(fullfile(statsDir, 'RELAXProcessingExtremeRejectionsAllParticipants_m002.mat'), ...
        'RELAXProcessingExtremeRejectionsAllParticipants');
    RELAXProcessing_wICA_AllParticipants = table({'m002_day1'}, 0.5, ...
        'VariableNames', {'aFileName', 'Proportion_artifactICs_reduced_by_wICA'}); %#ok<NASGU>
    save(fullfile(statsDir, 'ProcessingStatistics_wICA_m002.mat'), 'RELAXProcessing_wICA_AllParticipants');
    RELAX_unprocessed_files = table({'m002_day2.set'}, {'excluded'}, {'RELAX:excludeRecording'}, ...
        {'excessive bridging'}, datetime('now'), ...
        'VariableNames', {'File', 'Outcome', 'Identifier', 'Message', 'When'}); %#ok<NASGU>
    save(fullfile(statsDir, 'RELAX_unprocessed_files_m002.mat'), 'RELAX_unprocessed_files');
    touch(fullfile(cleanedDir, 'm002_day1_RELAX.set'));

    lastwarn('');
    [T, cfgs] = RELAX_QC_load_tables(statsDir, cleanedDir);
    [~, warnId] = lastwarn();
    assert(strcmp(warnId, 'RELAX_QC:missingStats'), 'a missing per-participant file should warn');

    assert(isequal(T.File, ["m001_day1"; "m001_day2"; "m002_day1"; "m002_day2"; "m002_day3"]), 'rows follow cfg.files');
    assert(isequal(T.Outcome, ["cleaned"; "cleaned"; "cleaned"; "excluded"; "notProcessed"]), 'outcomes');
    assert(T.Message(4) == "excessive bridging", 'exclusion message carried');
    assert(isequal(T.NChanImpedance(1:3), [2; 0; 1]), 'impedance channel counts');
    assert(isequal(T.NChanPop(1:2), [0; 1]), 'pop channel counts');
    assert(all(isnan(T.NChanBridge)), 'absent column gives NaN, not 0');
    assert(isequal(T.ExtremeExcisedSeconds(1:2), [1.5; 0]), '(500 + 250) samples at 500 Hz');
    assert(isequal(T.PopNSteps(1:2), [10; 20]), 'pop step counts');
    assert(isequal(T.WicaArtifactICsReduced(1:3), [0.3; 0.4; 0.5]), 'wICA proportion');
    assert(isequal(T.TableRawBlinkRatio{1}, [1; 3]), 'FileNumber 1 is column 1 of the metric array');
    assert(isequal(T.TableCleanedBlinkRatio{2}, [1.2; 1.4]), 'FileNumber 2 is column 2');
    assert(all(isnan(T.TableRawBlinkRatio{3})), 'missing RawMetrics gives NaN');
    assert(T.TableCleanedMuscleProportion(3) == 0.3, 'm002 day1 muscle proportion');
    assert(all(isnan(T.TableCleanedMuscleProportion(4:5))), 'non-cleaned rows carry no cleaned metrics');
    assert(isequal(sort(keys(cfgs)), {'m001', 'm002'}), 'one cfg per participant');

    testWarnsAboutCleanedFilesMissingFromTables();
    testNoConfigsRaisesClearError();
    fprintf('test_RELAX_QC_load_tables: PASS\n');
end

function testWarnsAboutCleanedFilesMissingFromTables()
    root = tempname; statsDir = fullfile(root, 'stats'); cleanedDir = fullfile(root, 'cleaned');
    mkdir(statsDir); mkdir(cleanedDir);
    cleanup = onCleanup(@() rmdir(root, 's'));
    RELAX_cfg = struct('files', {{'m001_day1.set'}}, 'DownSample_to_X_Hz', 500); %#ok<NASGU>
    save(fullfile(statsDir, 'RELAX_cfg_m001.mat'), 'RELAX_cfg');
    touch(fullfile(cleanedDir, 'm001_day1_RELAX.set'));
    touch(fullfile(cleanedDir, 'm009_day1_RELAX.set'));
    touch(fullfile(cleanedDir, 'm009_day2_RELAX.set'));

    previous = warning('off', 'RELAX_QC:missingStats'); restore = onCleanup(@() warning(previous));
    lastwarn('');
    RELAX_QC_load_tables(statsDir, cleanedDir);
    [message, warnId] = lastwarn();
    assert(strcmp(warnId, 'RELAX_QC:filesNotInTables'), 'cleaned files absent from the tables should warn, got %s', warnId);
    assert(contains(message, '2') && contains(message, 'm009_day1'), 'names the count and first files: %s', message);
end

function testNoConfigsRaisesClearError()
    root = tempname; statsDir = fullfile(root, 'stats'); mkdir(statsDir);
    cleanup = onCleanup(@() rmdir(root, 's'));
    try
        RELAX_QC_load_tables(statsDir, fullfile(root, 'cleaned'));
        error('test:noError', 'expected RELAX_QC:noConfigs');
    catch err
        assert(strcmp(err.identifier, 'RELAX_QC:noConfigs') && contains(err.message, statsDir), 'got %s: %s', err.identifier, err.message);
    end
end

function touch(f)
    fclose(fopen(f, 'w'));
end
