% RELAX_check_boundary_support() - Error if a boundary-blind pipeline stage is enabled
%                                  on data that contains boundary events.
%
% Usage:
%   >> RELAX_check_boundary_support(EEG, RELAX_cfg);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%   RELAX_cfg - RELAX configuration structure
%
% Outputs:
%   none. Throws error 'RELAX:boundaryUnsafeStage' if an unsafe stage is enabled.
%
% The stages listed here all process the record as one unbroken time series and
% would smear across discontinuities. See
% docs/superpowers/specs/2026-09-08-boundary-aware-preprocessing-design.md for the
% analysis of each. RELAX_cfg.Perform_ICA_subtract is deliberately NOT listed: it
% is ICA plus pop_subcomp, a purely spatial projection with no time-domain
% filtering and no wavelets, so it is unaffected by discontinuities.

function RELAX_check_boundary_support(EEG, RELAX_cfg)

    segments = RELAX_get_segments(EEG);
    if size(segments,1) <= 1
        return % no boundaries, nothing to guard against
    end

    unsafe = {};

    % Fields whose ABSENCE means "enabled". RELAX_Wrapper_slurm defaults an absent
    % Perform_PREP_bad_channel_rejection to 1 so that callers predating the flag keep
    % their original behaviour - which means an absent field here implies PREP WILL
    % run, and must be guarded. Every other flag counts as disabled when absent.
    enabledWhenAbsent = {'Perform_PREP_bad_channel_rejection'};

    flagChecks = { ...
        'Perform_PREP_bad_channel_rejection', 'PREP findNoisyChannels tiles its correlation windows flatly from sample 1, so windows straddle boundaries'; ...
        'Do_MWF_Once',                        'MWF delay-stacking builds one lagged covariance across the whole record, so lags span boundaries'; ...
        'Do_MWF_Twice',                       'MWF delay-stacking builds one lagged covariance across the whole record, so lags span boundaries'; ...
        'Do_MWF_Thrice',                      'MWF delay-stacking builds one lagged covariance across the whole record, so lags span boundaries'; ...
        'Run_GEDAI',                          'GEDAI has not been checked for boundary handling'; ...
        'Perform_wICA_on_ICLabel',            'RELAX_wICA_on_ICLabel_artifacts applies wavelet transforms across the whole record'; ...
        'ComputeSER_ARR',                     'the SER/ARR marking block calls RELAX_muscle, RELAX_horizontaleye, RELAX_drift and RELAX_pad_brief_mask_periods, all of which are boundary-blind'};

    for k = 1:size(flagChecks,1)
        thisField = flagChecks{k,1};
        if isfield(RELAX_cfg, thisField)
            thisEnabled = (RELAX_cfg.(thisField) == 1);
        else
            thisEnabled = any(strcmp(thisField, enabledWhenAbsent));
        end
        if thisEnabled
            unsafe{end+1} = sprintf('  RELAX_cfg.%s: %s', thisField, flagChecks{k,2}); %#ok<AGROW>
        end
    end

    if isfield(RELAX_cfg,'NotchFilterType') && strcmp(RELAX_cfg.NotchFilterType,'Butterworth')
        unsafe{end+1} = '  RELAX_cfg.NotchFilterType=''Butterworth'': the Butterworth notch filters across the whole record';
    end

    if ~isempty(unsafe)
        error('RELAX:boundaryUnsafeStage', ...
            ['This data contains %d boundary events (%d segments), but the following enabled ' ...
             'stages do not handle discontinuities:\n%s\n' ...
             'Either disable them, or implement boundary handling for them first. See ' ...
             'docs/superpowers/specs/2026-09-08-boundary-aware-preprocessing-design.md.'], ...
            size(segments,1)-1, size(segments,1), strjoin(unsafe, '\n'));
    end
end
