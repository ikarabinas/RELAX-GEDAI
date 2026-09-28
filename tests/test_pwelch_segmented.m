function test_pwelch_segmented()
% Boundary-aware Welch PSD used by the vendored zapline-plus.

    RELAX_add_vendored_zapline();

    fs     = 250;
    segLen = 1983;   % 7.93 s - the m208_dlpfc_day1 inter-train segment length
    nChan  = 4;

    % ---------------------------------------------------------------------
    % 1. Identical segments must produce the PSD of any one of them.
    % The pre-fix code assigned pwelch's frequency-vector output back into the
    % variable it passed as pwelch's nfft argument, so segments 2..N ran through
    % pwelch's Goertzel path, which omits the one-sided power doubling and
    % returns exactly half the PSD. That made the average 2/3 of the truth.
    % ---------------------------------------------------------------------
    rng(0);
    one    = randn(segLen, nChan);
    data   = [one; one; one];
    bounds = [segLen + 1; 2*segLen + 1];

    [pxxN, ~] = pwelch_segmented(data, [], [], 8192, fs, bounds);
    [pxx1, ~] = pwelch(one, hanning(segLen), [], 8192, fs);
    relErr = max(abs(pxxN(:) - pxx1(:))) / max(pxx1(:));
    assert(relErr < 1e-10, ...
        'identical segments must average to one segment''s PSD (rel err %.4g)', relErr);

    % ---------------------------------------------------------------------
    % 2. Segment PSDs are weighted by the number of samples they represent.
    % A longer segment must have proportionally more influence than a shorter
    % segment when their spectra differ.
    % ---------------------------------------------------------------------
    shortLen = 1000;
    longLen  = 3000;
    short   = randn(shortLen, nChan);
    long    = 2 * randn(longLen, nChan);
    unequal = [short; long];
    unequalBounds = shortLen + 1;
    [pxxWeighted, ~] = pwelch_segmented(unequal, [], [], 8192, fs, unequalBounds);
    [pxxShort, ~] = pwelch(short, hanning(shortLen), [], 8192, fs);
    [pxxLong, ~]  = pwelch(long, hanning(shortLen), [], 8192, fs);
    expected = (shortLen * pxxShort + longLen * pxxLong) / (shortLen + longLen);
    weightedErr = max(abs(pxxWeighted(:) - expected(:))) / max(expected(:));
    assert(weightedErr < 1e-10, ...
        'segment PSDs must be weighted by segment length (rel err %.4g)', weightedErr);

    % ---------------------------------------------------------------------
    % 3. An explicit nfft must give the +/-0.05 Hz fine-peak window real bins.
    % ---------------------------------------------------------------------
    nfft = 2^nextpow2(round(fs / (0.05/1.5)));   % 8192 at fs = 250
    [~, f] = pwelch_segmented(data, [], [], nfft, fs, bounds);
    assert(sum(f > 59.95 & f < 60.05) >= 1, ...
        'the +/-0.05 Hz fine-peak window must contain at least one bin');

    % ---------------------------------------------------------------------
    % 4. Regression witness: omitting nfft reintroduces the crash condition.
    % The window-derived grid on 7.93 s segments has bins at 59.937 and 60.059
    % and none in between. Kept so nobody "simplifies" the call sites back.
    % ---------------------------------------------------------------------
    [~, fDefault] = pwelch_segmented(data, [], [], [], fs, bounds);
    assert(sum(fDefault > 59.95 & fDefault < 60.05) == 0, ...
        'expected the window-derived grid to be too coarse - this is what passing nfft fixes');

    % ---------------------------------------------------------------------
    % 5. nfft below the window length must be raised, never wrap the data.
    % ---------------------------------------------------------------------
    [~, fSmall] = pwelch_segmented(data, [], [], 256, fs, bounds);
    dfSmall = fSmall(2) - fSmall(1);
    assert(dfSmall <= fs / 2^nextpow2(segLen) + eps, ...
        'nfft smaller than the window must be raised so pwelch zero-pads instead of wrapping');

    % ---------------------------------------------------------------------
    % 6. Continuous data must behave exactly like a plain full-window pwelch.
    % ---------------------------------------------------------------------
    cont     = randn(4000, nChan);
    [pc, ~]  = pwelch_segmented(cont, [], [], [], fs, []);
    [pr, ~]  = pwelch(cont, hanning(4000), [], 2^nextpow2(4000), fs);
    assert(isequal(size(pc), size(pr)) && max(abs(pc(:) - pr(:))) < 1e-12, ...
        'boundary-free data must match a plain full-window pwelch');

    % ---------------------------------------------------------------------
    % 7. A short fragment at a chunk edge no longer collapses the window: the
    % default window is sized from the segments of at least 1 s (see section 5),
    % so the grid stays as fine as the full segments allow. (Until 2026-09-18
    % this section documented the opposite as a known, deferred limitation.)
    % ---------------------------------------------------------------------
    ragged  = [randn(50, nChan); one; one];
    rbounds = [51; 51 + segLen];
    [~, fRagged] = pwelch_segmented(ragged, [], [], [], fs, rbounds);
    [~, fFull]   = pwelch_segmented(one, [], [], [], fs, []);
    assert(abs((fRagged(2) - fRagged(1)) - (fFull(2) - fFull(1))) < 1e-12, ...
        'a 50-sample fragment must not set the window (grid %.4g Hz vs %.4g Hz)', ...
        fRagged(2) - fRagged(1), fFull(2) - fFull(1));

    % ---------------------------------------------------------------------
    % 5. A sub-window sliver between two boundaries (an excision that landed a
    % sample away from a TMS boundary, m292 day 3) must neither size the window
    % nor error: the window comes from the segments that are actually usable.
    % ---------------------------------------------------------------------
    sliver = [one; randn(1, nChan); one];
    sliverBounds = [segLen + 1; segLen + 2];
    [pxxS, ~] = pwelch_segmented(sliver, [], [], 8192, fs, sliverBounds);
    relErr = max(abs(pxxS(:) - pxx1(:))) / max(pxx1(:));
    assert(relErr < 1e-10, 'a 1-sample segment must be ignored (rel err %.4g)', relErr);
    short = [one; randn(30, nChan); one];               % 0.12 s: too short to be useful
    [pxxT, ~] = pwelch_segmented(short, [], [], 8192, fs, [segLen + 1; segLen + 31]);
    relErr = max(abs(pxxT(:) - pxx1(:))) / max(pxx1(:));
    assert(relErr < 1e-10, 'a segment shorter than 1 s must not set the window (rel err %.4g)', relErr);

    fprintf('test_pwelch_segmented: PASS\n');
end
