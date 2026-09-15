function RELAX_add_vendored_noisetools()
% RELAX_ADD_VENDORED_NOISETOOLS - put the vendored NoiseTools closure on the path.
%
% RELAX_repair_steps calls nt_detrend. Rather than depend on a NoiseTools checkout
% living next to this repo, vendor/noisetools holds the 10-file runtime closure of
% nt_detrend (see vendor/noisetools/VENDORING.md).
%
% ORDERING MATTERS. vendor/zapline-plus vendors its own copies of five of these names
% - nt_demean, nt_fold, nt_unfold, nt_version and nt_greetings - and the first four
% are NOT byte-identical to ours. MATLAB resolves a name to whichever directory sits
% earlier on the path, so whichever of RELAX_add_vendored_noisetools and
% RELAX_add_vendored_zapline is called LAST wins for those shared names. This function
% must therefore be called AFTER RELAX_add_vendored_zapline (which is the order
% RELAX_SET_PARAMETERS_AND_RUN_slurm.m uses). Calling it first leaves nt_detrend
% running against zapline-plus's nt_demean/nt_fold/nt_unfold rather than the closure
% it was vendored with.
%
% The check below verifies the WHOLE closure, not just nt_detrend, so a shadowed
% dependency fails loudly here instead of silently changing nt_detrend's numerics.
% It runs immediately after this function's own addpath, which prepends vendorDir, so
% it can only fail if something is genuinely wrong (a file missing from the vendored
% closure, or a shadow that addpath cannot outrank) - not merely because another
% NoiseTools copy exists somewhere on the path.

    vendorDir = fullfile(fileparts(mfilename('fullpath')), 'vendor', 'noisetools');
    if ~isfolder(vendorDir)
        error('RELAX_add_vendored_noisetools:missing', ...
            'Vendored NoiseTools not found at %s', vendorDir);
    end

    addpath(vendorDir);

    % Derive the list from the directory rather than hardcoding it, so the check cannot
    % drift out of step with the vendored closure when it is re-synced.
    listing = dir(fullfile(vendorDir, '*.m'));
    for k = 1:numel(listing)
        [~, fnName] = fileparts(listing(k).name);
        resolved = which(fnName);
        if ~strncmp(resolved, vendorDir, numel(vendorDir))
            error('RELAX_add_vendored_noisetools:shadowed', ...
                ['%s resolves to\n  %s\ninstead of the vendored copy in\n  %s\n' ...
                 'Another NoiseTools copy is earlier on the path - remove it, or move ' ...
                 'this call later so its addpath wins. Note vendor/zapline-plus ships ' ...
                 'its own non-identical nt_demean, nt_fold, nt_unfold and nt_version.'], ...
                fnName, resolved, vendorDir);
        end
    end
end
