% RELAX_backup_blink_electrodes() - The net row one higher than the default blink electrodes,
%                                   used when fewer than two of RELAX_cfg.BlinkElectrodes survive.
%
% Usage:
%   >> labels = RELAX_backup_blink_electrodes();
%
% Outputs:
%   labels - {5 x 1} electrode labels

function labels = RELAX_backup_blink_electrodes()
    labels = {'E12'; 'E20'; 'E27'; 'E34'; 'E21'};
end
