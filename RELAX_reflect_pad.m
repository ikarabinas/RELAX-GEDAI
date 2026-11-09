% RELAX_reflect_pad() - Even (mirror) reflection padding about the edge samples.
%
% Usage:
%   >> padded = RELAX_reflect_pad(x, p);
%
% Inputs:
%   x - [nChannels x nSamples] data
%   p - number of samples to add at each end. Clamped to nSamples-1.
%
% Outputs:
%   padded - [nChannels x (nSamples + 2*p)] data
%
% Mirrors the interior about each edge sample, flip(interior), without repeating it.
% Not the odd reflection filtfilt applies internally, 2*x(edge) - flip(interior): that
% sets the pad's level from the single edge sample, so when the edge sample is far from
% the local level (mains at an arbitrary phase, or a residual TMS spike at a segment
% start) the pad is offset by twice that distance, and a high-pass rings on the step for
% hundreds of ms into the segment. Even reflection keeps the pad at the local level at
% the cost of a slope cusp at the edge. Measurement (m468, 0.5 Hz high-pass):
% TMS_response_variability implementation_notes/decay/line-noise-nuisance.md.

function padded = RELAX_reflect_pad(x, p)

    nSamples = size(x,2);
    p = min(p, nSamples-1);

    if p < 1
        padded = x;
        return
    end

    leftPad  = x(:, p+1:-1:2);
    rightPad = x(:, end-1:-1:end-p);
    padded = [leftPad, x, rightPad];
end
