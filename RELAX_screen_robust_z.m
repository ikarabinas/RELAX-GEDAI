% RELAX_screen_robust_z() - Median/MAD z-scores of one screening metric.
%
% Usage:
%   >> [z, centre, scale] = RELAX_screen_robust_z(x);
%
% Inputs:
%   x      - [n x 1] metric values, NaN where missing
%
% Outputs:
%   z      - [n x 1] (x - centre) / scale; NaN wherever x is NaN, and everywhere when
%            scale is 0
%   centre - median of x
%   scale  - 1.4826 * MAD; when MAD is 0, 1.2533 * mean absolute deviation instead.
%            NaN with centre when every x is NaN

function [z, centre, scale] = RELAX_screen_robust_z(x)

    centre = median(x, 'omitnan');
    absoluteDeviation = abs(x - centre);
    medianAbsoluteDeviation = median(absoluteDeviation, 'omitnan');
    if medianAbsoluteDeviation == 0
        scale = 1.2533 * mean(absoluteDeviation, 'omitnan');
    else
        scale = 1.4826 * medianAbsoluteDeviation;
    end

    if scale == 0
        z = NaN(size(x));
    else
        z = (x - centre) / scale;
    end
end
