function test_GEDAI_alpha_paths()
    p0 = path; restorePath = onCleanup(@() path(p0));
    d = GEDAI_alpha_paths();
    assert(startsWith(which('GEDAI'), d), 'GEDAI resolves to %s, not the be4e13a copy', which('GEDAI'));
    assert(~contains(path, '/gitRepos/GEDAI-master'), 'mid-merge GEDAI checkout is on the path');
    assert(abs(nargin('GEDAI')) >= 16, 'GEDAI lacks the broadband_only argument');
    opts = GEDAI_alpha_figures('opts');
    assert(isequal(opts.gedaiArgs([1 2 3 5 6 11 12 14]), {'auto', 300, 4.0, false, false, 'auto-', 1.0, true}));
    assert(isequal(size(opts.bands), [6 2]) && opts.edgeSec == 0.3 && opts.minSegSec == 2);
end
