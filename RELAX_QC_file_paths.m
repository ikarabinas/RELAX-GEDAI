% RELAX_QC_file_paths() - Input, Extremes_Rejected and cleaned .set paths of one file.
%
% Usage:
%   >> paths = RELAX_QC_file_paths(processed, inputs, file);
%
% inputs is dir(fullfile(myPath, '**', '*.set')) without the RELAXProcessed folder.

function paths = RELAX_QC_file_paths(processed, inputs, file)
    hit = find(strcmp({inputs.name}, char(file + ".set")), 1);
    if isempty(hit), error('RELAX_QC:missingInput', 'No input file %s.set under myPath', file); end
    paths.input    = fullfile(inputs(hit).folder, inputs(hit).name);
    paths.extremes = char(fullfile(processed, 'Extremes_Rejected', file + "_Extremes_Rejected.set"));
    paths.cleaned  = char(fullfile(processed, 'Cleaned_Data', file + "_RELAX.set"));
end
