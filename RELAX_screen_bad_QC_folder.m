% RELAX_screen_bad_QC_folder() - Where the move stage puts the files that fail the residual screen.
%
% Usage:
%   >> badDir = RELAX_screen_bad_QC_folder(opts);

function badDir = RELAX_screen_bad_QC_folder(opts)
    badDir = fullfile(opts.myPath, opts.outputFolderName, 'Cleaned_Data', 'bad_QC_files');
end
