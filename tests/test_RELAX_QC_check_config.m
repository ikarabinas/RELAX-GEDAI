function test_RELAX_QC_check_config()
% The baseline rebuild is only valid when nothing but extreme-period excision, average
% reference and wICA ran after Extremes_Rejected was saved.

    good = struct('Do_MWF_Once', 0, 'Do_MWF_Twice', 0, 'Do_MWF_Thrice', 0, 'Run_GEDAI', 0, ...
        'LowPassFilterBeforeMWF', 'yes', 'saveextremesrejected', 1, 'DownSample', 'yes', ...
        'Perform_targeted_wICA', 1, 'InterpolateRejectedElectrodesAfterCleaning', 'yes', 'ParticipantTag', 'm001');
    RELAX_QC_check_config(good);
    RELAX_QC_check_config(rmfield(good, 'Do_MWF_Thrice'));   % absent MWF flag means off

    bad = {'Do_MWF_Once', 1; 'Run_GEDAI', 1; 'LowPassFilterBeforeMWF', 'no'; 'saveextremesrejected', 0; ...
           'DownSample', 'no'; 'Perform_targeted_wICA', 0; 'InterpolateRejectedElectrodesAfterCleaning', 'no'};
    for k = 1:size(bad, 1)
        cfg = good; cfg.(bad{k, 1}) = bad{k, 2};
        try
            RELAX_QC_check_config(cfg);
            error('test:noError', '%s = %s should be rejected', bad{k, 1}, string(bad{k, 2}));
        catch err
            assert(strcmp(err.identifier, 'RELAX_QC:unsupportedConfig'), 'wrong error for %s: %s', bad{k, 1}, err.message);
            assert(contains(err.message, 'm001'), 'error names the participant');
        end
    end

    gedai = good; gedai.Run_GEDAI = 1; gedai.GEDAI_BroadbandOnly = true;
    RELAX_QC_check_config(gedai);   % broadband-only GEDAI runs after wICA, so the baseline holds
    gedai.GEDAI_BroadbandOnly = false;
    try
        RELAX_QC_check_config(gedai);
        error('test:noError', 'GEDAI with wavelet bands should be rejected');
    catch err
        assert(strcmp(err.identifier, 'RELAX_QC:unsupportedConfig'), 'wrong error: %s', err.message);
    end
    fprintf('test_RELAX_QC_check_config: PASS\n');
end
