function test_pwelch_segmented()
% Boundary-aware Welch PSD used by the vendored zapline-plus.

    RELAX_add_vendored_zapline();

    fs     = 250;
    segLen = 1983;   % 7.93 s - the m208_dlpfc_day1 inter-train segment length
    nChan  = 4;

    % ---------------------------------------------------------------------
    % 1. Every segment must contribute with equal weight.
    % Three identical segments must average to the PSD of any one of them.
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
    % 2. An explicit nfft must give the +/-0.05 Hz fine-peak window real bins.
    % ---------------------------------------------------------------------
    nfft = 2^nextpow2(round(fs / (0.05/1.5)));   % 8192 at fs = 250
    [~, f] = pwelch_segmented(data, [], [], nfft, fs, bounds);
    assert(sum(f > 59.95 & f < 60.05) >= 1, ...
        'the +/-0.05 Hz fine-peak window must contain at least one bin');

    % ---------------------------------------------------------------------
    % 3. Regression witness: omitting nfft reintroduces the crash condition.
    % The window-derived grid on 7.93 s segments has bins at 59.937 and 60.059
    % and none in between. Kept so nobody "simplifies" the call sites back.
    % ---------------------------------------------------------------------
    [~, fDefault] = pwelch_segmented(data, [], [], [], fs, bounds);
    assert(sum(fDefault > 59.95 & fDefault < 60.05) == 0, ...
        'expected the window-derived grid to be too coarse - this is what passing nfft fixes');

    % ---------------------------------------------------------------------
    % 4. nfft below the window length must be raised, never wrap the data.
    % ---------------------------------------------------------------------
    [~, fSmall] = pwelch_segmented(data, [], [], 256, fs, bounds);
    dfSmall = fSmall(2) - fSmall(1);
    assert(dfSmall <= fs / 2^nextpow2(segLen) + eps, ...
        'nfft smaller than the window must be raised so pwelch zero-pads instead of wrapping');

    % ---------------------------------------------------------------------
    % 5. Continuous data must behave exactly like a plain full-window pwelch.
    % ---------------------------------------------------------------------
    cont     = randn(4000, nChan);
    [pc, ~]  = pwelch_segmented(cont, [], [], [], fs, []);
    [pr, ~]  = pwelch(cont, hanning(4000), [], 2^nextpow2(4000), fs);
    assert(isequal(size(pc), size(pr)) && max(abs(pc(:) - pr(:))) < 1e-12, ...
        'boundary-free data must match a plain full-window pwelch');

    % ---------------------------------------------------------------------
    % 6. Known limitation, characterised deliberately: the window still
    % collapses to the SHORTEST segment, so a short fragment at a chunk edge
    % degrades the whole chunk's spectrum. Deferred - see the plan's "Deferred"
    % section. Update this test when that is addressed.
    % ---------------------------------------------------------------------
    ragged  = [randn(50, nChan); one; one];
    rbounds = [51; 51 + segLen];
    [~, fRagged] = pwelch_segmented(ragged, [], [], [], fs, rbounds);
    assert((fRagged(2) - fRagged(1)) > 0.5, ...
        'documenting that a 50-sample fragment still collapses the window');

    fprintf('test_pwelch_segmented: PASS\n');
end
