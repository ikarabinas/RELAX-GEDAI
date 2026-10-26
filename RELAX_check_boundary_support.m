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
% IMPORTANT: this only sees discontinuities present in the INPUT file at the time it
% is called (RELAX_Wrapper_slurm.m calls it once, immediately after pop_loadset). It
% is NOT re-checked later in the run: eeg_eegrej adds more boundaries mid-pipeline
% (rejecting extreme periods), and stages that run after that point are not covered
% by this guard even if this call raised no objection.
%
% The stages listed here all process the record as one unbroken time series and
% would smear across discontinuities. See
% docs/superpowers/specs/2026-09-08-boundary-aware-preprocessing-design.md for the
% analysis of each. RELAX_cfg.Perform_ICA_subtract is deliberately NOT listed: it
% is ICA plus pop_subcomp, a purely spatial projection with no time-domain
% filtering and no wavelets, so it is unaffected by discontinuities. The Butterworth
% notch is also deliberately NOT listed: RELAX_filtbutter's bandstop branch derives
% its reflection-pad width from the stop-band width (a notch rings at its bandwidth,
% not its centre frequency), so it is boundary-aware when RELAX_cfg.RespectBoundaries
% is passed through, as the wrapper does.

function RELAX_check_boundary_support(EEG, RELAX_cfg)

    segments = RELAX_get_segments(EEG);
    if size(segments,1) <= 1
        return % no boundaries, nothing to guard against
    end

    unsafe = {};

    % Fields whose ABSENCE means "enabled". RELAX_Wrapper_slurm defaults an absent
    % Perform_PREP_bad_channel_rejection to 1 so that callers predating the flag keep
    % their original behaviour - which means an absent field here implies PREP WILL
    % run, and must be guarded.
    %
    % This is NOT symmetric with wrapper behaviour for the rest of flagChecks below.
    % Only Perform_PREP_bad_channel_rejection and ComputeSER_ARR are actually
    % defaulted by the wrapper (ComputeSER_ARR to 0, at RELAX_Wrapper_slurm.m:401).
    % The others - Do_MWF_Once/Twice/Thrice, Run_GEDAI, Perform_wICA_on_ICLabel - are
    % dereferenced directly by the wrapper (e.g. `if RELAX_cfg.Do_MWF_Once==1`), so an
    % absent field there means the wrapper CRASHES on an undefined-field error, not
    % that the stage runs "off". Treating them as disabled here is still correct
    % (nothing unsafe can run if the wrapper never gets that far), just not for the
    % reason "absent means off" - it's "absent means the run fails before reaching it".
    enabledWhenAbsent = {'Perform_PREP_bad_channel_rejection'};

    flagChecks = { ...
        'Perform_PREP_bad_channel_rejection', 'PREP findNoisyChannels tiles its correlation windows flatly from sample 1, so windows straddle boundaries'; ...
        'Do_MWF_Once',                        'MWF delay-stacking builds one lagged covariance across the whole record, so lags span boundaries'; ...
        'Do_MWF_Twice',                       'MWF delay-stacking builds one lagged covariance across the whole record, so lags span boundaries'; ...
        'Do_MWF_Thrice',                      'MWF delay-stacking builds one lagged covariance across the whole record, so lags span boundaries'; ...
        'Run_GEDAI',                          'GEDAI has only been checked for boundary handling in broadband-only mode (RELAX_cfg.GEDAI_BroadbandOnly)'; ...
        'Perform_wICA_on_ICLabel',            'RELAX_wICA_on_ICLabel_artifacts applies wavelet transforms across the whole record'; ...
        'ComputeSER_ARR',                     'the SER/ARR marking block calls RELAX_muscle, RELAX_horizontaleye, RELAX_drift and RELAX_pad_brief_mask_periods, all of which are boundary-blind'};

    for k = 1:size(flagChecks,1)
        thisField = flagChecks{k,1};
        if strcmp(thisField, 'Run_GEDAI') && isfield(RELAX_cfg, 'GEDAI_BroadbandOnly') && RELAX_cfg.GEDAI_BroadbandOnly
            thisEnabled = false; % the boundary-aware broadband stage is the one GEDAI_alpha_run already uses on boundary-containing sessions
        elseif isfield(RELAX_cfg, thisField)
            thisEnabled = (RELAX_cfg.(thisField) == 1);
        else
            thisEnabled = any(strcmp(thisField, enabledWhenAbsent));
        end
        if thisEnabled
            unsafe{end+1} = sprintf('  RELAX_cfg.%s: %s', thisField, flagChecks{k,2}); %#ok<AGROW>
        end
    end

    if ~isempty(unsafe)
        error('RELAX:boundaryUnsafeStage', ...
            ['The INPUT file contains %d boundary events (%d segments), but the following enabled ' ...
             'stages do not handle discontinuities:\n%s\n' ...
             'Either disable them, or implement boundary handling for them first. See ' ...
             'docs/superpowers/specs/2026-09-08-boundary-aware-preprocessing-design.md. ' ...
             'Note this check only covers boundaries already present in the input file: ' ...
             'eeg_eegrej adds more boundaries later in the run, which this check cannot see.'], ...
            size(segments,1)-1, size(segments,1), strjoin(unsafe, '\n'));
    end
end
