function gof = GEDAI_alpha_dipfit(topo, Lref)
% GEDAI_ALPHA_DIPFIT  Goodness of fit of the best single free-orientation dipole to a topography.
%
%   gof = GEDAI_alpha_dipfit(topo, Lref)
%
%   Lref is nch x 3*nsrc with columns grouped by 3. gof = 1 - residual/total after removing
%   the mean from topo and from each leadfield column.

    p = inputParser;
    addRequired(p, 'topo', @isnumeric);
    addRequired(p, 'Lref', @isnumeric);
    parse(p, topo, Lref);

    topo = topo(:) - mean(topo(:));
    Lc = Lref - mean(Lref, 1);
    bestResidual = inf;
    for src = 1:size(Lref, 2) / 3
        G = Lc(:, 3*src-2:3*src);
        r = topo - G * (G \ topo);
        bestResidual = min(bestResidual, r' * r);
    end
    gof = 1 - bestResidual / (topo' * topo);
end
