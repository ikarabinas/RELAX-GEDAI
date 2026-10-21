function gedaiDir = GEDAI_alpha_paths()
% GEDAI_ALPHA_PATHS  Put the exported GEDAI be4e13a copy first on the path.
%
%   gedaiDir = GEDAI_alpha_paths()
%
%   Removes every path entry under gitRepos/GEDAI-master (the mid-merge checkout), adds the
%   copy at ~/scratch/gedai_review/alpha_figures/GEDAI-be4e13a with its subfolders, and returns
%   the folder holding GEDAI.m in the form `which` reports it.

    copyDir = fullfile(getenv('HOME'), 'scratch', 'gedai_review', 'alpha_figures', 'GEDAI-be4e13a');
    assert(isfile(fullfile(copyDir, 'GEDAI.m')), 'GEDAI_alpha:missingCopy', ...
        'GEDAI copy not found at %s', copyDir);

    entries = strsplit(path, pathsep);
    stale = entries(contains(entries, 'gitRepos/GEDAI-master'));
    if ~isempty(stale)
        rmpath(strjoin(stale, pathsep));
    end

    addpath(genpath(copyDir), '-begin');
    gedaiDir = fileparts(which('GEDAI'));
end
