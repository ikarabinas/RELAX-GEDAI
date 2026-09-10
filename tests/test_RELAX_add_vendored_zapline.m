function test_RELAX_add_vendored_zapline()
% The vendored zapline-plus copies must shadow the EEGLAB plugin's.

    RELAX_add_vendored_zapline();

    vendorDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'vendor', 'zapline-plus');

    % The mirror is complete, so EVERY zapline-plus function must resolve from it -
    % not just the customized ones. A partial shadow is the failure mode that would
    % pair a patched clean_data_with_zapline_plus with the plugin's unpatched nt_smooth.
    mFiles = dir(fullfile(vendorDir, '*.m'));
    % 26 upstream .m files. Task 2 adds pwelch_segmented.m, making it 27; the
    % floor stays at 26 so this test is correct both before and after that.
    assert(numel(mFiles) >= 26, ...
        'expected the full upstream tree (26 .m files), found %d', numel(mFiles));
    for k = 1:numel(mFiles)
        [~, fname] = fileparts(mFiles(k).name);
        resolved = which(fname);
        assert(strncmp(resolved, vendorDir, numel(vendorDir)), ...
            '%s resolved to %s, expected the vendored copy in %s', fname, resolved, vendorDir);
    end

    % Calling it twice must be harmless
    RELAX_add_vendored_zapline();
    assert(strncmp(which('nt_smooth'), vendorDir, numel(vendorDir)), 'must be idempotent');

    fprintf('test_RELAX_add_vendored_zapline: PASS\n');
end
