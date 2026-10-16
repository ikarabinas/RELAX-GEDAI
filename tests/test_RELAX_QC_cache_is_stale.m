function test_RELAX_QC_cache_is_stale()
% A cleaned file newer than the cache means the pipeline was rerun after extraction.
    d = tempname; mkdir(d); cleanup = onCleanup(@() rmdir(d, 's'));
    fclose(fopen(fullfile(d, 'x_RELAX.set'), 'w'));
    cache.meta.created = datetime('now') - hours(1);
    assert(RELAX_QC_cache_is_stale(cache, d), 'older cache is stale');
    cache.meta.created = datetime('now') + minutes(1);
    assert(~RELAX_QC_cache_is_stale(cache, d), 'newer cache is fresh');
    assert(~RELAX_QC_cache_is_stale(cache, fullfile(d, 'missing')), 'no cleaned folder, nothing to be stale against');
    fprintf('test_RELAX_QC_cache_is_stale: PASS\n');
end
