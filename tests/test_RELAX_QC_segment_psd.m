function test_RELAX_QC_segment_psd()
% Welch windows must never span a boundary: a DC jump at a boundary, or a short noisy
% segment that cannot hold a whole window, must leave no trace in the spectrum.

    srate = 500;
    t = (0:1999) / srate;
    sine = 3 * sin(2 * pi * 10 * t);
    data = [sine, sine + 1000, 500 * randn(1, 200)];
    segments = [1 2000; 2001 4000; 4001 4200];

    [freqs, power] = RELAX_QC_segment_psd(data, srate, segments, 2, 80);
    assert(abs(freqs(2) - freqs(1) - 0.5) < 1e-12, 'resolution is 1/windowSeconds');
    assert(freqs(1) == 0 && freqs(end) == 80, 'grid runs 0..maxHz');
    [~, peak] = max(power);
    assert(freqs(peak) == 10, '10 Hz peak');
    assert(max(power(freqs >= 30)) < 1e-6 * max(power), 'no broadband leak from the jump or the short segment');
    total = sum(power) * (freqs(2) - freqs(1));
    assert(abs(total - 4.5) / 4.5 < 0.02, 'density integrates to the sine variance (3^2/2), got %g', total);

    [~, twoChannels] = RELAX_QC_segment_psd([data; 2 * data], srate, segments, 2, 80);
    assert(max(abs(twoChannels(2, :) - 4 * twoChannels(1, :))) < 1e-9 * max(twoChannels(2, :)), 'per-channel rows');

    [~, none] = RELAX_QC_segment_psd(randn(2, 300), srate, [1 300], 2, 80);
    assert(all(isnan(none), 'all'), 'no whole window gives NaN');

    fprintf('test_RELAX_QC_segment_psd: PASS\n');
end
