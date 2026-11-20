% RELAX_screen_assert_no_moved_files() - Refuse to score while bad_QC_files holds moved files:
%                                        scoring would count them as not cleaned and leave them out of
%                                        the medians and scales the remaining files are judged against.
%
% Usage:
%   >> RELAX_screen_assert_no_moved_files(opts);

function RELAX_screen_assert_no_moved_files(opts)
    badDir = RELAX_screen_bad_QC_folder(opts);
    if ~isempty(dir(fullfile(badDir, '*_RELAX*.set')))
        error('RELAX_screen:movedFilesPresent', ...
              ['%s holds files an earlier move stage took out of Cleaned_Data. Move them back before ' ...
               'scoring again, or the screen judges the remaining files against a different baseline.'], badDir);
    end
end
