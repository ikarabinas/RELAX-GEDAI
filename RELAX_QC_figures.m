% RELAX_QC_figures() - Artifact-removal QC figure set for a RELAX_Wrapper_slurm.m run.
%                      Design: docs/superpowers/specs/2026-10-15-qc-figures-design.md
%
% Usage:
%   >> RELAX_QC_figures('extract')   % stage 1, slow: submit with sbatch run_RELAX_QC_figures.sh extract
%   >> RELAX_QC_figures('examples')  % rebuild only the Figure 1 example windows in the cache
%   >> RELAX_QC_figures('draw')      % stage 2, fast: redraw every figure from the cache
%   >> RELAX_QC_figures('all')       % both
%   >> RELAX_QC_figures('extract', {'m144_dlpfc_day1_stitched'})   % trial on a few files
%
% Edit myPath and overrides below for a different dataset or hand-picked examples.

function RELAX_QC_figures(stage, onlyFiles)

    myPath = '/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Oct26/rest_treat_joined/';
    % Hand-picked Figure 1 examples, replacing the automatic choice when non-empty, e.g.
    %   struct('file', {'m144_dlpfc_day1_stitched'}, 'kind', {'blink'}, 'centreSeconds', {412.3})
    % centreSeconds NaN keeps the automatic centre for that kind.
    overrides = struct('file', {}, 'kind', {}, 'centreSeconds', {});

    if nargin < 1, stage = 'all'; end
    if nargin < 2, onlyFiles = {}; end

    outDir = fullfile(myPath, 'RELAXProcessed', 'QC_figures');
    if ~isfolder(outDir), mkdir(outDir); end
    opts = struct('myPath', myPath, 'cachePath', fullfile(outDir, 'qc_cache.mat'), ...
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
end

function rebuildExamples(opts)
    loaded = load(opts.cachePath, 'cache');
    cache = loaded.cache;
    processed = fullfile(opts.myPath, 'RELAXProcessed');
    [~, cfgs] = RELAX_QC_load_tables(fullfile(processed, 'PerParticipantStats'), fullfile(processed, 'Cleaned_Data'));
    inputs = dir(fullfile(opts.myPath, '**', '*.set'));
    inputs = inputs(~contains({inputs.folder}, 'RELAXProcessed'));
    cache.examples = RELAX_QC_build_examples(cache.picks, cache.tables, cache.summaries, cfgs, processed, inputs, opts);
    save(opts.cachePath, 'cache', '-v7.3');
    fprintf('RELAX_QC_figures: %d example file(s) rebuilt in %s\n', numel(cache.examples), opts.cachePath);
end

function drawAll(opts, outDir)
    loaded = load(opts.cachePath, 'cache');
    cache = loaded.cache;
    RELAX_QC_warn_if_partial(cache, opts.cachePath);
    if RELAX_QC_cache_is_stale(cache, fullfile(opts.myPath, 'RELAXProcessed', 'Cleaned_Data'))
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
