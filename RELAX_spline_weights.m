% RELAX_spline_weights() - Spherical-spline interpolation weights, as an explicit linear map.
%
% Usage:
%   >> w = RELAX_spline_weights(basisXYZ, targetXYZ);
%   >> estimate = w * data(basisChannels, samples);
%
% Reproduces EEGLAB's eeg_interp(EEG, bad, 'spherical') - Perrin et al. (1989) spherical
% splines with EEGLAB's own conventions (coordinates normalised to the unit sphere, the
% 1 - Euclidean-distance argument to the Legendre series, m = 4, 7 terms, no
% regularisation, and the per-time-point mean over the basis channels removed before the
% fit and added back after) - but returns the weights instead of applying them, so the same
% interpolation can be applied to a few samples of one channel without building an EEG
% struct, and one solve serves every interval that shares a basis set.
%
% Inputs:
%   basisXYZ  - [nBasis x 3] Cartesian coordinates (chanlocs X, Y, Z) of the channels the
%               estimate is built from
%   targetXYZ - [nTarget x 3] coordinates of the channel(s) to estimate
%
% Outputs:
%   w - [nTarget x nBasis]; each row sums to 1, so a constant field is reproduced exactly.
%       For basis values V ([nBasis x T]), eeg_interp would return exactly w * V.

function w = RELAX_spline_weights(basisXYZ, targetXYZ)

    m = 4; nTerms = 7;                          % eeg_interp 'spherical': params = [0 4 7]

    B = basisXYZ ./ sqrt(sum(basisXYZ.^2, 2));  % unit sphere
    T = targetXYZ ./ sqrt(sum(targetXYZ.^2, 2));
    nB = size(B, 1);

    Gbb = computeg(B, B, m, nTerms);
    Gtb = computeg(T, B, m, nTerms);

    % eeg_interp: C = pinv([G; 1']) * [values - mean; 0]; result = Gtb * C + mean.
    % With P = pinv([Gbb; ones]) that is linear in the values:
    %   result = (Gtb * P(:,1:nB)) * (values - mean) + mean
    P  = pinv([Gbb; ones(1, nB)]);
    w0 = Gtb * P(:, 1:nB);
    % fold the mean removal / restoration into the weights: mean = (1/nB) * sum(values)
    w  = w0 + (1 - sum(w0, 2)) / nB;
end

function g = computeg(A, B, m, nTerms)
% EEGLAB's computeg: Legendre series in (1 - distance) between every A row and B row.
    nA = size(A, 1); nB = size(B, 1);
    dist = sqrt((A(:,1) - B(:,1)').^2 + (A(:,2) - B(:,2)').^2 + (A(:,3) - B(:,3)').^2);
    EI = 1 - dist;
    g = zeros(nA, nB);
    for n = 1:nTerms
        L = legendre(n, EI);                    % [n+1 x nA x nB]; order-0 term is L(1,:,:)
        g = g + ((2*n + 1) / (n^m * (n + 1)^m)) * reshape(L(1, :, :), nA, nB);
    end
    g = g / (4 * pi);
end
