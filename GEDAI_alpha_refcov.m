function [refCOV, Lref] = GEDAI_alpha_refcov(labels, leadfieldDir, leadfieldKey)
% GEDAI_ALPHA_REFCOV  Referenced leadfield covariance for the surviving channels, as the wrapper builds it.
%
%   [refCOV, Lref] = GEDAI_alpha_refcov(labels, leadfieldDir, leadfieldKey)
%
%   Selects the rows of the participant leadfield by channel label (in the order of labels),
%   re-references over the nch+1 channel montage (L - sum(L,1)/(nch+1)), and returns
%   refCOV = Lref*Lref'. Lref is nch x 3*nsrc, columns grouped by 3 per source.
%   Errors with GEDAI_alpha:leadfieldChannels when a label is missing from the leadfield.

    p = inputParser;
    addRequired(p, 'labels', @iscellstr);
    addRequired(p, 'leadfieldDir', @ischar);
    addRequired(p, 'leadfieldKey', @ischar);
    parse(p, labels, leadfieldDir, leadfieldKey);

    [L, leadfieldLabels] = RELAX_load_ppt_leadfield(leadfieldDir, [leadfieldKey '_x.set']);
    [found, rows] = ismember(labels(:)', leadfieldLabels);
    if ~all(found)
        error('GEDAI_alpha:leadfieldChannels', 'Channels missing from the %s leadfield: %s', ...
            leadfieldKey, strjoin(labels(~found), ', '));
    end

    nch = numel(labels);
    Lref = L(rows, :);
    Lref = Lref - sum(Lref, 1) / (nch + 1);
    refCOV = Lref * Lref';
end
