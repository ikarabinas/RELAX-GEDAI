% RELAX_QC_cache_is_stale() - True when any cleaned file is newer than the QC cache, i.e.
%                             the pipeline was rerun after extraction.

function stale = RELAX_QC_cache_is_stale(cache, cleanedDir)
    files = dir(fullfile(cleanedDir, '*_RELAX*.set'));
    stale = ~isempty(files) && any([files.datenum] > datenum(cache.meta.created));
end
