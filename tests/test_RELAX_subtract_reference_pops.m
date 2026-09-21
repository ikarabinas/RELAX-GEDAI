function test_RELAX_subtract_reference_pops()
% Subtracting a detected reference pop restores every channel's own signal: exact before
% the onset, the reference waveform removed through the transition, the amplitude removed
% from the offset to the end of the segment. Other segments, flat channels and a channel's
% own pop are untouched. Applying the same correction to an auxiliary EEG gives the same
% rows. After the repair, neither detector sees the reference pop any more.

    srate = 1000; n = 20000; rng(9); N = 64; noise = 10;

    function EEG = mkEEG(data, boundaries)
        if nargin < 2; boundaries = []; end
        EEG = struct();
        EEG.data = data; EEG.nbchan = size(data,1); EEG.pnts = size(data,2);
        EEG.srate = srate; EEG.trials = 1; EEG.xmin = 0; EEG.xmax = (EEG.pnts-1)/srate;
        EEG.event = struct('type', {}, 'latency', {});
        for b = boundaries(:)'
            EEG.event(end+1) = struct('type', 'boundary', 'latency', b);
        end
        EEG.chanlocs = struct('labels', arrayfun(@(k) sprintf('E%d', k), 1:size(data,1), 'UniformOutput', false));
    end
    % Exactly 0 before onset, exactly `amplitude` from offset on, a smooth cubic in between:
    % a sigmoid's tails never reach their plateaus, and section 1 asserts exact recovery.
    function r = smoothStep(onset, offset, amplitude)
        r = zeros(1, n); r(offset:end) = amplitude;
        u = ((onset+1:offset-1) - onset) / (offset - onset);
        r(onset+1:offset-1) = amplitude * (3*u.^2 - 2*u.^3);
    end

    S = noise * randn(N, n) + 30 * sin(2*pi*0.05*(1:n)/srate + 2*pi*rand(N,1));
    r = smoothStep(7960, 8040, 80) + smoothStep(11940, 12060, -120);   % two reference pops, one segment
    X = S - r;                                                          % channels see -r: -80 then +120
    X(5, 15000:end) = X(5, 15000:end) + 300;                          % E5's own pop
    X(1, :) = 0;                                                        % E1 is flat

    %% 1. hand-built steps with an exact reference trace: near-perfect recovery
    W = 150;
    refTrace = -r;
    refSteps = struct('onset', [7960; 11940], 'offset', [8040; 12060], 'amplitude', [-80; 120]);
    EEG = mkEEG(X);
    [out, info] = RELAX_subtract_reference_pops(EEG, refSteps, refTrace, 'WindowMs', W);
    assert(info.available && info.nSteps == 2 && info.nChannelsCorrected == N - 1);
    assert(isequal(info.skippedLabels, {'E1'}));
    truth = S; truth(5, 15000:end) = truth(5, 15000:end) + 300;
    live = setdiff(1:N, 1);
    resid = out.data(live, :) - truth(live, :);
    assert(max(abs(resid(:))) < 1e-6, 'exact trace + exact steps restores s_i (max %.2g)', max(abs(resid(:))));
    assert(isequal(out.data(1, :), X(1, :)), 'flat channel untouched');
    assert(isequal(out.data(:, 1:7960), X(:, 1:7960)), 'bitwise identical before the first onset');
    assert(abs(info.transitionSeconds - ((8040-7960-1) + (12060-11940-1)) / srate) < 1e-12);
    assert(abs(info.medianAmplitude - 100) < 1e-12, 'median |amplitude|');
    assert(isa(info.correction, 'single') && isequal(size(info.correction), [1 n]));
    assert(isequal(size(out.data), size(X)) && strcmp(class(out.data), class(X)));

    %% 2. a boundary after the first pop: only that segment's remainder is corrected
    EEG2 = mkEEG(X, 10000.5);
    out2 = RELAX_subtract_reference_pops(EEG2, struct('onset', 7960, 'offset', 8040, 'amplitude', -80), refTrace, 'WindowMs', W);
    assert(isequal(out2.data(:, 10001:end), X(:, 10001:end)), 'the next segment is untouched');
    resid2 = out2.data(live, 8040:10000) - (X(live, 8040:10000) + 80);
    assert(max(abs(resid2(:))) < 1e-6);

    %% 3. the same correction on an auxiliary EEG equals the matching rows
    aux = mkEEG(X([7 8], :)); aux.chanlocs = EEG.chanlocs([7 8]);
    outAux = RELAX_subtract_reference_pops(aux, refSteps, refTrace, 'WindowMs', W);
    assert(isequal(outAux.data, out.data([7 8], :)));

    %% 4. empty steps: nothing changes
    e = struct('onset', zeros(0,1), 'offset', zeros(0,1), 'amplitude', zeros(0,1));
    [out4, info4] = RELAX_subtract_reference_pops(EEG, e, refTrace);
    assert(isequal(out4.data, X) && info4.available && info4.nSteps == 0 && isnan(info4.medianAmplitude));

    %% 5. single-precision data stays single and is still corrected
    EEGs = mkEEG(single(X));
    outS = RELAX_subtract_reference_pops(EEGs, refSteps, refTrace, 'WindowMs', W);
    assert(isa(outS.data, 'single'));
    assert(max(abs(double(outS.data(live, :)) - truth(live, :)), [], 'all') < 1e-2);

    %% 6. end to end with the detector: detect -> subtract -> nothing left, E5's pop still there
    [refSteps6, refInfo6] = RELAX_detect_reference_pops(EEG);
    assert(refInfo6.nSteps == 2, 'both reference pops detected (got %d)', refInfo6.nSteps);
    [out6, info6] = RELAX_subtract_reference_pops(EEG, refSteps6, refInfo6.refChannel.data, ...
        'WindowMs', refInfo6.detectInfo.params.WindowMs);
    assert(info6.nSteps == 2);
    resid6 = out6.data(live, :) - truth(live, :);
    assert(max(abs(resid6(:))) < 3 * noise, 'detected steps restore s_i to within the noise (max %.1f)', max(abs(resid6(:))));
    after = out6; after.data(1, :) = 0;
    [~, refInfoAfter] = RELAX_detect_reference_pops(after);
    assert(refInfoAfter.nSteps == 0, 'no reference pop left (found %d)', refInfoAfter.nSteps);
    [steps, pinfo] = RELAX_detect_steps(after, 'CommonModeChannelFraction', 0.15);
    assert(pinfo.nCommonModeEvents == 0, 'no common-mode event left (found %d)', pinfo.nCommonModeEvents);
    assert(any(steps.channel == 5 & abs(steps.sample - 15000) < 100), 'E5''s own pop is still there for the per-channel stage');

    fprintf('test_RELAX_subtract_reference_pops: PASS\n');
end