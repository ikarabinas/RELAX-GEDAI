% RELAX_QC_figures() - Artifact-removal QC figure set for a RELAX_Wrapper_slurm.m run.
%                      Design: docs/superpowers/specs/2026-10-15-qc-figures-design.md
%
% Usage:
%   >> RELAX_QC_figures('extract')   % stage 1, slow: submit with sbatch run_RELAX_QC_figures.sh extract
%   >> RELAX_QC_figures('examples')  % rebuild only the Figure 1 example windows in the cache
%   >> RELAX_QC_figures('draw')      % stage 2, fast: redraw every figure from the cache
%   >> RELAX_QC_figures('drawMatched') % redraw from the cache, only for files also cleaned in comparisonFolderName
%   >> RELAX_QC_figures('all')       % both
%   >> RELAX_QC_figures('extract', {'m144_dlpfc_day1_stitched'})   % trial on a few files
%
% Edit myPath and overrides below for a different dataset or hand-picked examples.

function RELAX_QC_figures(stage, onlyFiles)

    myPath = '/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Oct26/rest_treat_joined/';
    outputFolderName = 'RELAXProcessed_broadbandGEDAI'; % the RELAX_cfg.OutputFolderName of the run being plotted
    comparisonFolderName = 'RELAXProcessed'; % the 'drawMatched' stage keeps only files cleaned in this run too
    % Hand-picked Figure 1 examples, replacing the automatic choice when non-empty, e.g.
    %   struct('file', {'m144_dlpfc_day1_stitched'}, 'kind', {'blink'}, 'centreSeconds', {412.3})
    % centreSeconds NaN keeps the automatic centre for that kind.
    % Pinned to the examples and windows of the original RELAXProcessed figure set, so the
    % two sets compare directly. Pop windows stay automatic (the pop stage runs before GEDAI,
    % so it finds the same pops).
    overrides = struct('file', {'m358_dlpfc_day4_stitched', 'm358_dlpfc_day4_stitched', ...
                                'm497_dlpfc_day2_stitched', 'm497_dlpfc_day2_stitched', ...
                                'm159_dlpfc_day2_stitched', 'm159_dlpfc_day2_stitched', ...
                                'm275_dlpfc_day5_stitched', 'm275_dlpfc_day5_stitched'}, ...
                       'kind', {'blink', 'muscle', 'blink', 'muscle', 'blink', 'muscle', 'blink', 'muscle'}, ...
                       'centreSeconds', {1179.222, 612.090, 55.752, 1030.846, 147.016, 676.414, 1001.436, 811.496});

    if nargin < 1, stage = 'all'; end
    if nargin < 2, onlyFiles = {}; end

    outDir = fullfile(myPath, outputFolderName, 'QC_figures');
    if ~isfolder(outDir), mkdir(outDir); end
    opts = struct('myPath', myPath, 'outputFolderName', outputFolderName, 'cachePath', fullfile(outDir, 'qc_cache.mat'), ...
                  'overrides', overrides, 'onlyFiles', {onlyFiles}, ...
                  'windowSeconds', 2, 'maxHz', 80, 'blinkHalfSeconds', 0.5, 'topoHalfSeconds', 0.05, ...
                  'muscleWindowSeconds', 2, 'exampleSeconds', 6, 'minBlinksForGroup', 5);

    if any(strcmp(stage, {'extract', 'all'}))
        RELAX_QC_extract(opts);
    end
    if strcmp(stage, 'examples')
        rebuildExamples(opts);
    end
    if any(strcmp(stage, {'draw', 'all'}))
        drawAll(opts, outDir);
    end
    if strcmp(stage, 'drawMatched')
        drawMatched(opts, fullfile(myPath, outputFolderName, 'QC_figures_matched'), fullfile(myPath, comparisonFolderName, 'Cleaned_Data'));
    end
end

function rebuildExamples(opts)
    loaded = load(opts.cachePath, 'cache');
    cache = loaded.cache;
    processed = fullfile(opts.myPath, opts.outputFolderName);
    [~, cfgs] = RELAX_QC_load_tables(fullfile(processed, 'PerParticipantStats'), fullfile(processed, 'Cleaned_Data'));
    inputs = dir(fullfile(opts.myPath, '**', '*.set'));
    inputs = inputs(~contains({inputs.folder}, 'RELAXProcessed'));  % also excludes RELAXProcessed_* siblings
    cache.examples = RELAX_QC_build_examples(cache.picks, cache.tables, cache.summaries, cfgs, processed, inputs, opts);
    save(opts.cachePath, 'cache', '-v7.3');
    fprintf('RELAX_QC_figures: %d example file(s) rebuilt in %s\n', numel(cache.examples), opts.cachePath);
end

function drawMatched(opts, outDir, comparisonCleanedDir)
    if ~isfolder(outDir), mkdir(outDir); end
    loaded = load(opts.cachePath, 'cache');
    cache = loaded.cache;
    comparisonFiles = regexprep(string({dir(fullfile(comparisonCleanedDir, '*_RELAX*.set')).name}), '_RELAX(_GEDAI)?\.set$', '');
    keep = ismember(cache.tables.File, comparisonFiles);
    fprintf('RELAX_QC_figures: %d of %d cached files were also cleaned in %s\n', sum(keep), numel(keep), comparisonCleanedDir);
    cache.tables = cache.tables(keep, :);
    cache.summaries = cache.summaries(keep);
    drawCache(cache, opts, outDir);
end

function drawAll(opts, outDir)
    loaded = load(opts.cachePath, 'cache');
    drawCache(loaded.cache, opts, outDir);
end

function drawCache(cache, opts, outDir)
    RELAX_QC_warn_if_partial(cache, opts.cachePath);
    if RELAX_QC_cache_is_stale(cache, fullfile(opts.myPath, opts.outputFolderName, 'Cleaned_Data'))
        warning('RELAX_QC:staleCache', 'Cleaned_Data has files newer than %s; rerun the extract stage.', opts.cachePath);
    end
    captions = [RELAX_QC_fig_examples(cache, outDir); RELAX_QC_fig_spectra(cache, outDir); ...
                RELAX_QC_fig_paired_metrics(cache, outDir); RELAX_QC_fig_data_taken(cache, outDir); ...
                RELAX_QC_fig_blink_topo(cache, outDir)];
    writeCaptions(cache, captions, opts.cachePath, outDir);
    fprintf('RELAX_QC_figures: figures and captions.md written to %s\n', outDir);
end

function writeCaptions(cache, captions, cachePath, outDir)
    fid = fopen(fullfile(outDir, 'captions.md'), 'w');
    fprintf(fid, '%s\n', strjoin([captions(:); RELAX_QC_captions_footer(cache, cachePath)], newline));
    fclose(fid);
end
