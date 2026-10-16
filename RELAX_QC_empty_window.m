% RELAX_QC_empty_window() - Template for one Figure 1 example window.

function w = RELAX_QC_empty_window(kind)
    w = struct('kind', string(kind), 'available', false, 'note', "", 'channels', {{}}, ...
               't', [], 'before', [], 'after', [], 'removed', [], 'shade', [], 'boundaries', [], 'excised', zeros(0, 2));
end
