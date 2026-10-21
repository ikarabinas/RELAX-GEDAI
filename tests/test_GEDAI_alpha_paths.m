function test_GEDAI_alpha_paths()
    p0 = path; restorePath = onCleanup(@() path(p0));
    expected = char(java.io.File(fullfile(getenv('HOME'), 'scratch', 'gedai_review', 'alpha_figures', 'GEDAI-be4e13a')).getCanonicalPath);
    d = GEDAI_alpha_paths();
    assert(strcmp(d, expected), 'GEDAI_alpha_paths returned %s, expected %s', d, expected);
    assert(startsWith(which('GEDAI'), expected), 'GEDAI resolves to %s, not the be4e13a copy', which('GEDAI'));
    assert(~contains(path, '/gitRepos/GEDAI-master'), 'mid-merge GEDAI checkout is on the path');
    assert(contains(fileread(which('GEDAI')), 'broadband_only'), 'GEDAI source does not mention broadband_only');
    assert(abs(nargin('GEDAI')) >= 16, 'GEDAI lacks the broadband_only argument');
    opts = GEDAI_alpha_figures('opts');
    assert(isequal(opts.gedaiArgs([1 2 3 5 6 11 12 14]), {'auto', 300, 4.0, false, false, 'auto-', 1.0, true}));
    assert(isequal(size(opts.bands), [6 2]) && opts.edgeSec == 0.3 && opts.minSegSec == 2);
end
