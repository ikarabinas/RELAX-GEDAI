function RELAX_add_vendored_noisetools()
% RELAX_ADD_VENDORED_NOISETOOLS - put the vendored NoiseTools closure on the path.
%
% RELAX_repair_steps calls nt_detrend. Rather than depend on a NoiseTools checkout
% living next to this repo, vendor/noisetools holds the 10-file runtime closure of
% nt_detrend (see vendor/noisetools/VENDORING.md). Unlike RELAX_add_vendored_zapline
% this shadows nothing - NoiseTools is not an EEGLAB plugin - so the call order
% relative to eeglab() does not matter.

    vendorDir = fullfile(fileparts(mfilename('fullpath')), 'vendor', 'noisetools');
    if ~isfolder(vendorDir)
        error('RELAX_add_vendored_noisetools:missing', ...
            'Vendored NoiseTools not found at %s', vendorDir);
    end

    addpath(vendorDir);

    resolved = which('nt_detrend');
    if ~strncmp(resolved, vendorDir, numel(vendorDir))
        error('RELAX_add_vendored_noisetools:shadowed', ...
            ['nt_detrend resolves to\n  %s\ninstead of the vendored copy in\n  %s\n' ...
             'Another NoiseTools copy is earlier on the path - remove it, or move this ' ...
             'call later so its addpath wins.'], resolved, vendorDir);
    end
end
