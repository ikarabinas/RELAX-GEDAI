function test_RELAX_detect_reference_pops()
% Reference pops are steps on the cross-channel median. A common sigmoid step subtracted
% from every channel is found once, with its interval and amplitude, and passes the
% uniformity gate; a single-channel pop, a mixed-sign movement event and a same-sign
% variable-amplitude event are not reported as reference pops.

    srate = 1000; n = 20000; rng(5); N = 64; noise = 10;

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
    function r = sigmoidStep(at, amplitude, riseSamples)
        r = zeros(1, n);
        t = (1:n) - at;
        r = r + amplitude ./ (1 + exp(-t / (riseSamples / 6)));
    end
    cols = {'sample','onset','offset','amplitude','riseMs','z','channelZ', ...
            'signAgreement','amplitudeSpread','fractionNear'};

    % drift: a slow sine of a different phase per channel, so channels reorder over time
    drift = 30 * sin(2*pi*0.05*(1:n)/srate + 2*pi*rand(N,1));
    S = noise * randn(N, n) + drift;
    r = sigmoidStep(8000, 80, 20);                           % reference pops +80 uV at 8 s
    X = S - r;
    X(5, 14000:end) = X(5, 14000:end) + 300;                 % one contact pop on E5 at 14 s
    move = zeros(N, n);                                      % movement at 17 s: 20 channels,
    mvCh = 10:29; mvAmp = (100 + 300*rand(20,1)) .* sign(randn(20,1));   % mixed sign
    move(mvCh, 17000:end) = repmat(mvAmp, 1, n - 17000 + 1);
    X = X + move;

    %% 1. the reference pop is found, once, with the right interval and amplitude
    [refSteps, refInfo] = RELAX_detect_reference_pops(mkEEG(X));
    assert(refInfo.available, 'available: %s', refInfo.status);
    for c = cols
        assert(isnumeric(refSteps.(c{1})) && size(refSteps.(c{1}), 2) == 1, 'refSteps.%s must be a column', c{1});
    end
    assert(refInfo.nSteps == 1 && numel(refSteps.sample) == 1, 'expected 1 reference pop, got %d (candidates %d, rejected %d)', ...
        refInfo.nSteps, refInfo.nCandidates, refInfo.nRejectedUniformity);
    assert(abs(refSteps.sample - 8000) <= 30, 'mid-crossing at %d', refSteps.sample);
    assert(refSteps.onset < 8000 && refSteps.offset > 8000 && refSteps.offset - refSteps.onset < 200);
    assert(abs(refSteps.amplitude - (-80)) < 12, 'amplitude %.1f (channels see -r, so -80)', refSteps.amplitude);
    assert(refSteps.signAgreement >= 0.95 && refSteps.amplitudeSpread <= 0.25 && refSteps.fractionNear >= 0.75, ...
        'uniformity: agreement %.2f spread %.2f near %.2f', refSteps.signAgreement, refSteps.amplitudeSpread, refSteps.fractionNear);
    assert(isequal(size(refInfo.perChannelAmplitude), [N 1]));
    assert(all(abs(refInfo.perChannelAmplitude(setdiff(1:N, 5)) - (-80)) < 4 * noise));
    assert(abs(refInfo.stepsPerMinute - 1 / (n / srate / 60)) < 1e-9);
    assert(numel(refInfo.candidates.accepted) == refInfo.nCandidates);

    %% 2. nothing near the single-channel pop or the mixed-sign movement
    assert(~any(abs(refSteps.sample - 14000) < 300), 'a single-channel pop is not a reference pop');
    assert(~any(abs(refSteps.sample - 17000) < 300), 'a mixed-sign movement event is not a reference pop');

    %% 3. same-sign, variable-amplitude event: a candidate on the median, rejected by the gate
    X3 = S - r;
    gnd = zeros(N, n); gCh = 1:40; gAmp = 50 + 350*rand(40,1);
    gnd(gCh, 12000:end) = repmat(gAmp, 1, n - 12000 + 1);
    [refSteps3, refInfo3] = RELAX_detect_reference_pops(mkEEG(X3 + gnd));
    assert(refInfo3.nSteps == 1 && abs(refSteps3.sample - 8000) <= 30, 'only the true reference pop survives');
    assert(refInfo3.nRejectedUniformity >= 1, 'the variable-amplitude event was a candidate and was rejected');
    assert(~any(abs(refSteps3.sample - 12000) < 300));
    rejected = refInfo3.candidates;
    k = find(abs(rejected.sample - 12000) < 300, 1);
    assert(~isempty(k) && ~rejected.accepted(k) && rejected.amplitudeSpread(k) > 0.5);

    %% 4. gates off: the same event is then reported
    [refSteps4, ~] = RELAX_detect_reference_pops(mkEEG(X3 + gnd), 'MinSignAgreement', 0, 'MaxAmplitudeSpread', Inf);
    assert(any(abs(refSteps4.sample - 12000) < 300));

    %% 5. Detect options are forwarded: a 200 uV floor hides the 80 uV pop
    [refSteps5, refInfo5] = RELAX_detect_reference_pops(mkEEG(X), 'Detect', {'MinAmplitudeMicrovolts', 200});
    assert(refInfo5.available && refInfo5.nSteps == 0 && isequal(size(refSteps5.sample), [0 1]));
    assert(refInfo5.detectInfo.params.MinAmplitudeMicrovolts == 200);

    %% 6. boundaries: the pop is still found when the record is split at 10 s
    [refSteps6, ~] = RELAX_detect_reference_pops(mkEEG(X, 10000.5));
    assert(numel(refSteps6.sample) == 1 && abs(refSteps6.sample - 8000) <= 30);

    %% 7. quiet data: available, empty, 0x1 columns
    [refSteps7, refInfo7] = RELAX_detect_reference_pops(mkEEG(S));
    assert(refInfo7.available && refInfo7.nSteps == 0);
    for c = cols
        assert(isequal(size(refSteps7.(c{1})), [0 1]), 'refSteps.%s must be 0x1 when empty', c{1});
    end
    assert(isequal(size(refInfo7.perChannelAmplitude), [N 0]));

    %% 8. too few channels: unavailable, empty
    [refSteps8, refInfo8] = RELAX_detect_reference_pops(mkEEG(X(1:4, :)));
    assert(~refInfo8.available && ~isempty(refInfo8.status) && isequal(size(refSteps8.sample), [0 1]));

    %% 9. excluded labels are not used for the median or the uniformity vector
    [~, refInfo9] = RELAX_detect_reference_pops(mkEEG(X), 'ExcludeLabels', {'E5'});
    assert(refInfo9.nSteps == 1 && isnan(refInfo9.perChannelAmplitude(5)));

    %% 10. a staircase reference pop is returned as non-overlapping steps whose subtraction
    %      restores the channels; every candidate records the step it was merged into
    r10 = sigmoidStep(9000, 60, 20) + sigmoidStep(9040, 60, 20);
    X10 = S - r10;
    [refSteps10, refInfo10] = RELAX_detect_reference_pops(mkEEG(X10));
    assert(refInfo10.nSteps >= 1 && all(abs(refSteps10.sample - 9020) < 100), 'staircase found: %s', mat2str(refSteps10.sample'));
    [on, order] = sort(refSteps10.onset); off = refSteps10.offset(order);
    assert(all(on(2:end) > off(1:end-1)), 'returned steps never overlap: %s', mat2str([on off]));
    assert(abs(sum(refSteps10.amplitude) - (-120)) < 15, 'total amplitude %.1f', sum(refSteps10.amplitude));
    assert(isnumeric(refInfo10.nMerged) && refInfo10.nMerged == refInfo10.nCandidates - refInfo10.nRejectedUniformity - refInfo10.nSteps);
    assert(isequal(size(refInfo10.candidates.mergedInto), [refInfo10.nCandidates 1]));
    assert(all(refInfo10.candidates.mergedInto(refInfo10.candidates.accepted) >= 1) && ...
           all(refInfo10.candidates.mergedInto(~refInfo10.candidates.accepted) == 0));
    assert(isequal(size(refInfo10.perChannelAmplitude), [N refInfo10.nSteps]));
    assert(isequal(size(refInfo10.candidates.perChannelAmplitude), [N refInfo10.nCandidates]), ...
        'candidates.perChannelAmplitude must be per-candidate, not per-step');
    out10 = RELAX_subtract_reference_pops(mkEEG(X10), refSteps10, refInfo10.refChannel.data, ...
        'WindowMs', refInfo10.detectInfo.params.WindowMs);
    resid10 = out10.data - S;
    assert(max(abs(resid10(:))) < 3 * noise, 'subtracting the returned steps restores s_i (max %.1f)', max(abs(resid10(:))));
    assert(isequal(size(refInfo7.candidates.mergedInto), [0 1]) && refInfo7.nMerged == 0);

    %% 11. a rejected candidate is merged nowhere (the merge itself is pinned by
    %      test_RELAX_merge_overlapping_steps; the detector may or may not double-bracket the staircase)
    assert(isfield(refInfo3.candidates, 'mergedInto') && refInfo3.candidates.mergedInto(k) == 0, 'a rejected candidate is merged nowhere');

    fprintf('test_RELAX_detect_reference_pops: PASS\n');
end
