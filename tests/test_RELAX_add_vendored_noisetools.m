function test_RELAX_add_vendored_noisetools()
% The vendored NoiseTools closure is present, self-contained, and resolves on the path.

    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    vendorDir = fullfile(repoRoot, 'vendor', 'noisetools');

    assert(isfolder(vendorDir), 'vendor/noisetools must exist');

    expected = {'nt_detrend','nt_regw','nt_demean','nt_demean2','nt_fold', ...
                'nt_unfold','nt_vecadd','nt_vecmult','nt_greetings','nt_version'};
    for k = 1:numel(expected)
        assert(isfile(fullfile(vendorDir, [expected{k} '.m'])), ...
            'vendored file %s.m is missing', expected{k});
    end

    d = dir(fullfile(vendorDir, '*.m'));
    assert(numel(d) == numel(expected), ...
        'vendor/noisetools should hold exactly %d .m files, found %d', numel(expected), numel(d));

    % The helper puts nt_detrend on the path and points at the vendored copy
    RELAX_add_vendored_noisetools();
    resolved = which('nt_detrend');
    assert(strncmp(resolved, vendorDir, numel(vendorDir)), ...
        'nt_detrend resolves to %s, not the vendored copy in %s', resolved, vendorDir);

    % And it actually runs: a pure linear ramp is removed by an order-1 fit
    x = (1:1000)' * 0.5;
    y = nt_detrend(x, 1);
    assert(max(abs(y)) < 1e-6, 'order-1 detrend should flatten a linear ramp, max residual %g', max(abs(y)));

    fprintf('test_RELAX_add_vendored_noisetools: PASS\n');
end
