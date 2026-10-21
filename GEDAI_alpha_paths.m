function gedaiDir = GEDAI_alpha_paths()
% GEDAI_ALPHA_PATHS  Put the exported GEDAI be4e13a copy first on the path.
%
%   gedaiDir = GEDAI_alpha_paths()
%
%   Removes every path entry under gitRepos/GEDAI-master (the mid-merge checkout), adds the
%   copy at ~/scratch/gedai_review/alpha_figures/GEDAI-be4e13a with its subfolders, and returns
%   the canonical copy folder, after asserting that `which('GEDAI')` resolves inside it.

    linkedDir = fullfile(getenv('HOME'), 'scratch', 'gedai_review', 'alpha_figures', 'GEDAI-be4e13a');
    assert(isfile(fullfile(linkedDir, 'GEDAI.m')), 'GEDAI_alpha:missingCopy', ...
        'GEDAI copy not found at %s', linkedDir);
    copyDir = char(java.io.File(linkedDir).getCanonicalPath());

    entries = strsplit(path, pathsep);
    stale = entries(contains(entries, 'gitRepos/GEDAI-master'));
    if ~isempty(stale)
        rmpath(strjoin(stale, pathsep));
    end

    addpath(genpath(copyDir), '-begin');
    if ~startsWith(which('GEDAI'), copyDir)
        error('GEDAI_alpha:wrongGEDAI', 'GEDAI resolves to %s, not the copy in %s', which('GEDAI'), copyDir);
    end
    gedaiDir = copyDir;
end
