function test_zapline_boundary_crash()
% Regression for the reported m208_dlpfc_day1 failure: clean_data_with_zapline_plus
% died in its fine noise-peak search on data whose inter-boundary segments are ~8 s,
% because the frequency grid derived from those segments had no bin inside the
% +/-0.05 Hz search window.
%
%   Unable to perform assignment because the left and right sides have a
%   different number of elements.
%   Error in clean_data_with_zapline_plus (line 492)

    RELAX_add_vendored_zapline();

    fs     = 250;
    segLen = round(7.93 * fs);   % same 0.12207 Hz grid as the 500 Hz production case
    nSeg   = 20;                 % ~159 s, enough for several >=30 s adaptive chunks
    nChan  = 8;

    rng(7);
    t      = (0:segLen-1)' / fs;
    line   = 5 * sin(2*pi*60.02*t);
    data   = zeros(nSeg*segLen, nChan);
    bounds = zeros(nSeg-1, 1);
    for k = 1:nSeg
        rows = (k-1)*segLen + (1:segLen);
        data(rows, :) = repmat(line, 1, nChan) + randn(segLen, nChan);
        if k > 1
            bounds(k-1) = rows(1);
        end
    end

    cfg = struct('plotResults', 0, 'boundaries', bounds);
    [cleaned, cfgOut, analytics] = clean_data_with_zapline_plus(data, fs, cfg);

    assert(isequal(size(cleaned), size(data)), 'output size must match input');
    assert(all(isfinite(cleaned(:))), 'cleaned data must be finite');

    % The 60 Hz line must actually be attenuated
    [pRaw,   f] = pwelch(data(1:segLen,:),    hanning(segLen), [], 8192, fs);
    [pClean, ~] = pwelch(cleaned(1:segLen,:), hanning(segLen), [], 8192, fs);
    band = f > 59.5 & f < 60.5;
    attenDb = 10*log10(mean(mean(pRaw(band,:))) / mean(mean(pClean(band,:))));
    assert(attenDb > 6, 'expected >6 dB attenuation at 60 Hz, got %.2f dB', attenDb);

    assert(cfgOut.detailedFreqBoundsUpper(2) == 0.05, 'guard: default fine bounds assumed');
    assert(all(isfinite(analytics.proportionRemovedNoise)), ...
        'proportionRemovedNoise must be finite - it was NaN when the +/-0.05 Hz index was empty');

    fprintf('test_zapline_boundary_crash: PASS\n');
end
