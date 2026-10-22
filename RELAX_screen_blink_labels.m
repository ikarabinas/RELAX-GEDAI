% RELAX_screen_blink_labels() - The electrodes the blink method actually used for one file.
%
% Usage:
%   >> labels = RELAX_screen_blink_labels(cfg, fallback);
%
% Inputs:
%   cfg      - RELAX_cfg; uses BlinkElectrodes
%   fallback - BlinkElectrodeFallback covariate: 0 configured row, 1 backup row added,
%              2 no blink marking. NaN (covariate absent) is treated as 0
%
% Outputs:
%   labels - cell of labels; empty when fallback is 2

function labels = RELAX_screen_blink_labels(cfg, fallback)
    if isnan(fallback), fallback = 0; end
    switch fallback
        case 0
            labels = cfg.BlinkElectrodes(:);
        case 1
            labels = [cfg.BlinkElectrodes(:); RELAX_backup_blink_electrodes()];
        otherwise
            labels = {};
    end
end
