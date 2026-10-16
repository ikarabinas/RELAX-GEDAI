function chanlocs = make_qc_chanlocs(nChan)
% nChan electrodes spread over the upper hemisphere (Fibonacci spiral), labelled E1..En,
% with every coordinate system EEGLAB's interpolation and topoplot need.
    k = (0:nChan-1)';
    z = 0.05 + 0.95 * (1 - k / (nChan - 1));
    r = sqrt(1 - z.^2);
    phi = k * pi * (3 - sqrt(5));
    labels = arrayfun(@(i) sprintf('E%d', i), 1:nChan, 'UniformOutput', false);
    chanlocs = struct('labels', labels, 'X', num2cell((r .* cos(phi))'), ...
                      'Y', num2cell((r .* sin(phi))'), 'Z', num2cell(z'));
    chanlocs = convertlocs(chanlocs, 'cart2all');
end
