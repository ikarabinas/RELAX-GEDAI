function test_RELAX_check_boundary_support()
% Startup guard for stages that cannot handle boundary events.

    EEG_with    = struct('pnts', 200, 'event', struct('type', {'boundary'}, 'latency', {100.5}));
    EEG_without = struct('pnts', 200, 'event', struct('type', {'X'}, 'latency', {100}));

    safe_cfg = struct( ...
        'Perform_PREP_bad_channel_rejection', 0, ...
        'Do_MWF_Once', 0, 'Do_MWF_Twice', 0, 'Do_MWF_Thrice', 0, ...
        'Run_GEDAI', 0, ...
        'NotchFilterType', 'ZaplinePlus', ...
        'Perform_wICA_on_ICLabel', 0, ...
        'ComputeSER_ARR', 0, ...
        'Perform_ICA_subtract', 0);

    % A fully safe config must not throw, with or without boundaries
    RELAX_check_boundary_support(EEG_with, safe_cfg);
    RELAX_check_boundary_support(EEG_without, safe_cfg);

    % Each unsafe stage must throw when boundaries are present
    unsafe_fields = {'Perform_PREP_bad_channel_rejection', 'Do_MWF_Once', 'Do_MWF_Twice', ...
                     'Do_MWF_Thrice', 'Run_GEDAI', 'Perform_wICA_on_ICLabel', 'ComputeSER_ARR'};
    for k = 1:numel(unsafe_fields)
        cfg = safe_cfg;
        cfg.(unsafe_fields{k}) = 1;
        threw = false;
        try
            RELAX_check_boundary_support(EEG_with, cfg);
        catch err
            threw = true;
            assert(strcmp(err.identifier, 'RELAX:boundaryUnsafeStage'), ...
                'wrong error identifier for %s', unsafe_fields{k});
        end
        assert(threw, 'guard must throw when %s is enabled with boundaries present', unsafe_fields{k});

        % ...and must NOT throw on data with no boundaries
        RELAX_check_boundary_support(EEG_without, cfg);
    end

    % The Butterworth notch is unsafe; ZaplinePlus and PMnotch are not checked here
    cfg = safe_cfg;
    cfg.NotchFilterType = 'Butterworth';
    threw = false;
    try
        RELAX_check_boundary_support(EEG_with, cfg);
    catch err
        threw = true;
        assert(strcmp(err.identifier, 'RELAX:boundaryUnsafeStage'), 'wrong identifier for notch');
    end
    assert(threw, 'guard must throw for Butterworth notch with boundaries present');

    % Perform_ICA_subtract is spatial-only and must NEVER be guarded
    cfg = safe_cfg;
    cfg.Perform_ICA_subtract = 1;
    RELAX_check_boundary_support(EEG_with, cfg);

    % Missing fields must be tolerated (treated as disabled), not error
    RELAX_check_boundary_support(EEG_with, struct('Run_GEDAI', 0));

    fprintf('test_RELAX_check_boundary_support: PASS\n');
end
