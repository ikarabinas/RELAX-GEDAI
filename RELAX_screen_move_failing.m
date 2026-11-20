% RELAX_screen_move_failing() - Move stage of the residual-noise screen: move each file whose
%                               ScreenOutcome is "fail" from Cleaned_Data/ into Cleaned_Data/bad_QC_files/.
%
% Usage:
%   >> moved = RELAX_screen_move_failing(opts);
%
% opts: built by RELAX_residual_screen (.myPath, .outputFolderName, .outDir); needs the score
% stage's residual_screen.mat. Every Cleaned_Data file named <File>_RELAX* moves, so a .fdt goes
% with its .set. Passing, unscreened and notCleaned files stay. Returns the base names moved.
% Downstream loaders that list Cleaned_Data/*.set without recursing no longer see the moved files;
% the pipeline's skip-already-cleaned scan doesn't either, so a pipeline rerun reprocesses them.

function moved = RELAX_screen_move_failing(opts)

    load(fullfile(opts.outDir, 'residual_screen.mat'), 'screen');
    cleanedDir = fullfile(opts.myPath, opts.outputFolderName, 'Cleaned_Data');
    badDir = RELAX_screen_bad_QC_folder(opts);

    moved = strings(0, 1);
    for file = screen.File(screen.ScreenOutcome == "fail")'
        cleanedFiles = dir(fullfile(cleanedDir, file + "_RELAX*"));
        if isempty(cleanedFiles), continue, end
        if ~isfolder(badDir), mkdir(badDir); end
        for k = 1:numel(cleanedFiles)
            movefile(fullfile(cleanedDir, cleanedFiles(k).name), fullfile(badDir, cleanedFiles(k).name));
        end
        moved(end + 1, 1) = file; %#ok<AGROW>
        fprintf('RELAX_screen_move_failing: moved %s to %s\n', file, badDir);
    end
    fprintf('RELAX_screen_move_failing: %d failing file(s) moved\n', numel(moved));
end
