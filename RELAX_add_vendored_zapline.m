function RELAX_add_vendored_zapline()
% RELAX_ADD_VENDORED_ZAPLINE - put the vendored zapline-plus ahead of the EEGLAB plugin.
%
% The EEGLAB plugin copy of zapline-plus is boundary-aware (it passes boundaries through
% to zapline) but still crashes on data whose inter-boundary segments are short, because
% its pwelch frequency grid is not pinned (see
% docs/superpowers/plans/2026-09-09-zapline-nfft-grid-fix.md). vendor/zapline-plus
% holds rebased, fixed copies of the four files that needed changing. Call this AFTER
% eeglab(), because eeglab() prepends its own plugin directories to the path.

    vendorDir = fullfile(fileparts(mfilename('fullpath')), 'vendor', 'zapline-plus');
    if ~isfolder(vendorDir)
        error('RELAX_add_vendored_zapline:missing', 'Vendored zapline-plus not found at %s', vendorDir);
    end

    addpath(vendorDir);  % addpath prepends, so this wins over the plugin directory

    resolved = which('clean_data_with_zapline_plus');
    if ~strncmp(resolved, vendorDir, numel(vendorDir))
        error('RELAX_add_vendored_zapline:shadowed', ...
            ['clean_data_with_zapline_plus resolves to\n  %s\ninstead of the vendored copy in\n  %s\n' ...
             'Something re-prepended the EEGLAB plugin path after this call - move the ' ...
             'RELAX_add_vendored_zapline() call to after the last eeglab() invocation.'], resolved, vendorDir);
    end

    if isempty(which('RELAX_reflect_pad'))
        error('RELAX_add_vendored_zapline:missingDependency', ...
            ['RELAX_reflect_pad does not resolve on the path.\n' ...
             'vendor/zapline-plus/zapline_filter_by_segments.m calls RELAX_reflect_pad, so the ' ...
             'vendored tree depends on the RELAX-GEDAI repo root being on the path, not just on ' ...
             'MATLAB/EEGLAB toolboxes. Add the repo root to the path before calling ' ...
             'RELAX_add_vendored_zapline().']);
    end
end
