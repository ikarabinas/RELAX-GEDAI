% RELAX_QC_saved_pop_steps() - The pop steps the wrapper saved on an Extremes_Rejected
%                              EEG, as a struct of columns; [] when none were saved.
%
% Usage:
%   >> steps = RELAX_QC_saved_pop_steps(extremesEEG);

function steps = RELAX_QC_saved_pop_steps(extremesEEG)
    steps = [];
    if ~isfield(extremesEEG.RELAX, 'PopScreen') || ~isfield(extremesEEG.RELAX.PopScreen, 'steps'), return, end
    steps = extremesEEG.RELAX.PopScreen.steps;
    if istable(steps), steps = table2struct(steps, 'ToScalar', true); end
end
