function test_RELAX_spline_weights()
% Spherical-spline weight row: reproduces EEGLAB eeg_interp('spherical') as an explicit
% linear map, so a short stretch of one channel can be rebuilt from its neighbours.

    rng(11);
    % 32 electrodes on the upper hemisphere (Fibonacci lattice), unit sphere
    N = 32; k = (0:N-1)'; phi = pi * (3 - sqrt(5)) * k;
    zc = 1 - k / (N - 1) * 0.9;                    % z from 1 down to 0.1: no electrodes on the rim
    r = sqrt(1 - zc.^2);
    XYZ = [r .* cos(phi), r .* sin(phi), zc] * 85;  % 85 mm head radius, like real chanlocs

    % 1. weights are a proper interpolation: they sum to one
    w = RELAX_spline_weights(XYZ(2:end, :), XYZ(1, :));
    assert(isequal(size(w), [1 N-1]), 'weights must be a 1 x nBasis row');
    assert(abs(sum(w) - 1) < 1e-10, 'weights must sum to 1 (constant fields interpolate exactly)');

    % 2. a smooth field (dipole-like first-order harmonic) is recovered leave-one-out
    d = [0.3 0.5 0.8]; d = d / norm(d);
    field = (XYZ / 85) * d';                         % N x 1 value per electrode
    err = zeros(N, 1);
    for c = 1:N
        others = setdiff(1:N, c);
        w = RELAX_spline_weights(XYZ(others, :), XYZ(c, :));
        err(c) = w * field(others) - field(c);
    end
    inner = 1:N/2;                                   % the lattice's upper half; the lowest ring extrapolates
    assert(rms(err(inner)) / rms(field(inner)) < 0.05, ...
        'interior leave-one-out error %.3f of field rms', rms(err(inner)) / rms(field(inner)));
    assert(rms(err) / rms(field) < 0.35, 'whole-lattice error %.3f (rim extrapolation dominates)', rms(err) / rms(field));

    % 3. a multi-target call returns one row per target
    W = RELAX_spline_weights(XYZ(4:end, :), XYZ(1:3, :));
    assert(isequal(size(W), [3 N-3]));
    assert(all(abs(sum(W, 2) - 1) < 1e-10));

    % 4. matches EEGLAB's eeg_interp bit-for-bit (when EEGLAB is on the path)
    if exist('eeg_interp', 'file') == 2
        T = 50;
        data = randn(N, T) .* (1 + field) + field * linspace(-1, 1, T) * 10;
        EEG = eeg_emptyset;
        EEG.data = data; EEG.nbchan = N; EEG.pnts = T; EEG.trials = 1; EEG.srate = 100;
        EEG.xmin = 0; EEG.xmax = (T-1)/100;
        cl = struct('labels', {}, 'X', {}, 'Y', {}, 'Z', {}, 'theta', {}, 'radius', {}, ...
                    'sph_theta', {}, 'sph_phi', {}, 'sph_radius', {}, 'type', {}, 'urchan', {}, 'ref', {});
        for c = 1:N
            cl(c) = struct('labels', sprintf('E%d', c), 'X', XYZ(c,1), 'Y', XYZ(c,2), 'Z', XYZ(c,3), ...
                           'theta', 1, 'radius', 0.5, 'sph_theta', 1, 'sph_phi', 1, 'sph_radius', 85, ...
                           'type', 'EEG', 'urchan', c, 'ref', '');
        end
        EEG.chanlocs = cl;
        bad = 7;
        ref = evalc_quiet(@() eeg_interp(EEG, bad, 'spherical'));
        others = setdiff(1:N, bad);
        w = RELAX_spline_weights(XYZ(others, :), XYZ(bad, :));
        mine = w * data(others, :);
        assert(max(abs(mine - ref.data(bad, :))) < 1e-6, ...
            'weights must reproduce eeg_interp (max diff %.2e)', max(abs(mine - ref.data(bad, :))));
    else
        fprintf('  (eeg_interp not on path - skipping the EEGLAB cross-check)\n');
    end

    fprintf('test_RELAX_spline_weights: PASS\n');
end

function out = evalc_quiet(f)
% eeg_interp prints progress lines; keep the test output clean
    [~, out] = evalc('f()');
end
