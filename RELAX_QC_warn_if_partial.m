% RELAX_QC_warn_if_partial() - Warn when the cache is a mid-run checkpoint, not a finished
%                              extract; the draw stage still draws what the cache holds.

function RELAX_QC_warn_if_partial(cache, cachePath)
    if isfield(cache.meta, 'partial') && cache.meta.partial
        warning('RELAX_QC:partialCache', ['%s is a checkpoint from an unfinished extract run; figures show only the files ' ...
                'summarised so far and may lack example windows.'], cachePath);
    end
end
