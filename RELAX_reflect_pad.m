% RELAX_reflect_pad() - Odd (antisymmetric) reflection padding about the edge samples.
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
% Uses odd reflection, 2*x(edge) - flip(interior), which is the convention
% filtfilt applies internally. Even (mirror) reflection would place a cusp at the
% segment edge that a high pass filter then rings on; odd reflection continues the
% local trend across the edge instead, which matters here because the segments
% between boundary events are only ~8 s long.

function padded = RELAX_reflect_pad(x, p)

    nSamples = size(x,2);
    p = min(p, nSamples-1);

    if p < 1
        padded = x;
        return
    end

    leftPad  = 2*x(:,1)   - x(:, p+1:-1:2);
    rightPad = 2*x(:,end) - x(:, end-1:-1:end-p);
    padded = [leftPad, x, rightPad];
end
