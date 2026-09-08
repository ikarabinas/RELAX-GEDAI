function run_all_tests()
% Discover and run every tests/test_*.m, reporting PASS/FAIL per test.
% Errors at the end if any test failed, so `matlab -batch` exits non-zero.

    thisDir = fileparts(mfilename('fullpath'));
    testFiles = dir(fullfile(thisDir, 'test_*.m'));

    if isempty(testFiles)
        error('run_all_tests: no test_*.m files found in %s', thisDir);
    end

    failures = {};
    for k = 1:numel(testFiles)
        [~, testName] = fileparts(testFiles(k).name);
        try
            feval(testName);
        catch err
            failures{end+1} = sprintf('%s: %s', testName, err.message); %#ok<AGROW>
            fprintf('%s: FAIL - %s\n', testName, err.message);
        end
    end

    fprintf('\n%d/%d test files passed.\n', numel(testFiles) - numel(failures), numel(testFiles));

    if ~isempty(failures)
        error('run_all_tests: %d test file(s) failed:\n%s', numel(failures), strjoin(failures, '\n'));
    end
end
