% RELAX_QC_cleaned_suffix() - Suffix of the cleaned .set a run saved for one file:
%                            "_RELAX_GEDAI.set" when GEDAI ran, else "_RELAX.set".
%
% Usage:
%   >> suffix = RELAX_QC_cleaned_suffix(cleanedDir, file);

function suffix = RELAX_QC_cleaned_suffix(cleanedDir, file)
    suffix = "_RELAX.set";
    if isfile(fullfile(cleanedDir, file + "_RELAX_GEDAI.set"))
        suffix = "_RELAX_GEDAI.set";
    end
end
