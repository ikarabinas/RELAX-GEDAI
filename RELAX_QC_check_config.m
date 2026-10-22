% RELAX_QC_check_config() - Stop the QC extractor on any config where rebuilding the
%                           pre-wICA baseline from Extremes_Rejected would be wrong.
%
% Usage:
%   >> RELAX_QC_check_config(RELAX_cfg);
%
% The baseline is valid only when nothing but extreme-period excision, average reference
% and targeted wICA ran after Extremes_Rejected was saved, and the cleaned file was
% interpolated back to the full montage. GEDAI is allowed only in broadband-only mode,
% which runs after targeted wICA and so leaves the baseline valid. Errors RELAX_QC:unsupportedConfig otherwise.

function RELAX_QC_check_config(cfg)
    problems = strings(0, 1);
    if any([field(cfg, 'Do_MWF_Once', 0), field(cfg, 'Do_MWF_Twice', 0), field(cfg, 'Do_MWF_Thrice', 0)] ~= 0)
        problems(end+1) = "MWF is on";
    end
    if field(cfg, 'Run_GEDAI', 0) ~= 0 && ~field(cfg, 'GEDAI_BroadbandOnly', false)
        problems(end+1) = "GEDAI is on with its wavelet bands (only broadband-only GEDAI runs after the baseline's wICA step without altering it)";
    end
    if ~strcmp(field(cfg, 'LowPassFilterBeforeMWF', ''), 'yes')
        problems(end+1) = "the low-pass runs after Extremes_Rejected is saved (LowPassFilterBeforeMWF is not 'yes')";
    end
    if field(cfg, 'saveextremesrejected', 0) ~= 1, problems(end+1) = "Extremes_Rejected was not saved"; end
    if ~strcmp(field(cfg, 'DownSample', ''), 'yes'), problems(end+1) = "DownSample is off"; end
    if field(cfg, 'Perform_targeted_wICA', 0) ~= 1, problems(end+1) = "targeted wICA is off"; end
    if ~strcmp(field(cfg, 'InterpolateRejectedElectrodesAfterCleaning', ''), 'yes')
        problems(end+1) = "cleaned files are not interpolated to the full montage";
    end
    if ~isempty(problems)
        error('RELAX_QC:unsupportedConfig', 'QC figures do not support the config for %s: %s', ...
              char(string(field(cfg, 'ParticipantTag', 'unknown participant'))), strjoin(problems, '; '));
    end
end

function v = field(cfg, name, default)
    v = default;
    if isfield(cfg, name), v = cfg.(name); end
end
