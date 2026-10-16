% RELAX_QC_captions_footer() - Last line of captions.md: when and from what it was made,
%                              and which files were skipped for alignment.

function footer = RELAX_QC_captions_footer(cache, cachePath)
    isSkipped = [cache.summaries.status]' == "skipped";
    skippedFiles = string(cache.tables.File(isSkipped));
    reasons = [cache.summaries(isSkipped).skipReason]';
    count = nnz(isSkipped);
    if count == 0
        skipped = "No files were skipped for alignment.";
    else
        noun = "files"; if count == 1, noun = "file"; end
        skipped = sprintf('%d %s skipped for alignment: %s.', count, noun, strjoin(skippedFiles + " (" + reasons + ")", '; '));
    end
    footer = sprintf('_Generated %s from %s. %s_\n', string(datetime('now')), cachePath, skipped);
end
