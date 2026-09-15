function test_RELAX_repair_steps()
% Step repair: removes steps and drift, preserves signal above the cutoff, degrades safely.

    RELAX_add_vendored_noisetools();

    srate = 1000;
    n     = 10000;                       % 10 s
    rng(11);

    function EEG = mkEEG(data)
        EEG = struct();
        EEG.data     = data;
        EEG.nbchan   = size(data,1);
        EEG.pnts     = size(data,2);
        EEG.srate    = srate;
        EEG.trials   = 1;
        EEG.event    = struct('type', {}, 'latency', {});
        EEG.chanlocs = struct('labels', ...
            arrayfun(@(k) sprintf('E%d', k), 1:size(data,1), 'UniformOutput', false));
    end

    function s = mkSteps(ch, samp)
        % ch/sample are parallel arrays, one entry per step (RELAX_detect_steps's
        % contract). A scalar ch with multiple samples broadcasts to "same channel,
        % several steps" rather than silently truncating to one step.
        if isscalar(ch) && numel(samp) > 1
            ch = repmat(ch, numel(samp), 1);
        end
        s = struct('channel', ch(:), 'label', {cell(numel(ch),1)}, 'sample', samp(:), ...
                   'amplitude', zeros(numel(ch),1), 'z', zeros(numel(ch),1));
    end

    % ---- 1. a 500 uV step is removed ----
    x = randn(1, n) * 5;
    x(4001:end) = x(4001:end) + 500;
    EEG = mkEEG(x);
    [out, info] = RELAX_repair_steps(EEG, mkSteps(1, 4001));
    assert(info.available, 'repair should be available: %s', info.status);
    pre  = median(out.data(1, 3001:4000));
    post = median(out.data(1, 4001:5000));
    assert(abs(post - pre) < 10, ...
        'the 500 uV level difference should be gone, still %g uV', abs(post - pre));

    % ---- 2. THE DESIGN-CLAIM TEST: 2 Hz survives, 0.3 Hz does not ----
    % This is the test that certifies the ladder, not just the implementation. If the
    % measured c(n) constants were wrong, or the thresholds drifted out of sync with
    % MaxCutoffHz, nothing else in this suite notices.
    t = (0:n-1) / srate;
    for f = [2.0, 0.3]
        sig = sin(2*pi*f*t) * 100;
        [outF, ~] = RELAX_repair_steps(mkEEG(sig), mkSteps([], []));
        retained = std(outF.data(1,:)) / std(sig);
        if f == 2.0
            assert(retained > 0.95, ...
                '2 Hz must survive with >95%% amplitude, retained %.3f', retained);
        else
            assert(retained < 0.5, ...
                '0.3 Hz must be strongly attenuated, retained %.3f', retained);
        end
    end

    % ---- 3. slow drift is removed even with no steps at all ----
    drift = linspace(0, 300, n) + randn(1, n) * 5;
    [outD, infoD] = RELAX_repair_steps(mkEEG(drift), mkSteps([], []));
    assert(abs(median(outD.data(1, 9001:end)) - median(outD.data(1, 1:1000))) < 20, ...
        'drift must be removed on a channel with no steps');
    assert(infoD.subSegmentsPerChannel(1) >= 1, 'a stepless channel is still partitioned');

    % ---- 4. ladder selection at each duration threshold ----
    % Sub-segment durations are set by placing steps; check the chosen rung.
    % thresholds: >=2.000 s -> order3oa, >=1.635 -> order3, >=1.240 -> order2,
    %             >=0.845 -> order1, else meanOnly
    cases = { 3.000, 'order3oa'; 1.800, 'order3'; 1.400, 'order2'; ...
              1.000, 'order1';   0.500, 'meanOnly' };
    for k = 1:size(cases,1)
        segLen = round(cases{k,1} * srate);
        nSteps = floor(n / segLen) - 1;
        stepSamples = (1:nSteps) * segLen + 1;
        [~, infoK] = RELAX_repair_steps(mkEEG(randn(1,n)*5), mkSteps(ones(1,nSteps), stepSamples));
        assert(infoK.treatmentCounts.(cases{k,2}) > 0, ...
            'a %.3f s sub-segment should use rung %s, counts were %s', ...
            cases{k,1}, cases{k,2}, jsonencode(infoK.treatmentCounts));
    end

    % ---- 4a. overlap-add survives a MaxPolyOrder above the response table ----
    % chooseOrder saturates at numel(C)=3, so gating overlap-add on the RAW configured
    % MaxPolyOrder makes the equality unsatisfiable the moment MaxPolyOrder exceeds 3:
    % every long sub-segment then silently collapses to a single global fit whose cutoff
    % sits far below the 0.5 Hz floor the ladder guarantees. The gate must compare against
    % the ACHIEVABLE ceiling, min(MaxPolyOrder, numel(C)). This test fails if that is
    % reverted to `order == opt.MaxPolyOrder`.
    for mpo = [3 4 5]
        wOff = warning('off', 'RELAX_repair_steps:maxPolyOrderClamped');
        cleanupW = onCleanup(@() warning(wOff));
        [~, infoM] = RELAX_repair_steps(mkEEG(randn(1,n)*5), mkSteps([], []), ...
            'MaxPolyOrder', mpo);
        clear cleanupW
        assert(infoM.treatmentCounts.order3oa > 0, ...
            ['a 10 s stepless channel must use order-3 overlap-add at MaxPolyOrder=%d ' ...
             '(the achievable ceiling is 3), counts were %s'], ...
            mpo, jsonencode(infoM.treatmentCounts));
        assert(infoM.treatmentCounts.order3 == 0, ...
            ['MaxPolyOrder=%d must not fall back to a single global order-3 fit, ' ...
             'counts were %s'], mpo, jsonencode(infoM.treatmentCounts));
    end

    % Lowering MaxPolyOrder must report the order ACTUALLY used, not a hardcoded order3oa.
    wOff2 = warning('off', 'RELAX_repair_steps:maxPolyOrderClamped');
    [~, infoM2] = RELAX_repair_steps(mkEEG(randn(1,n)*5), mkSteps([], []), 'MaxPolyOrder', 2);
    warning(wOff2);
    assert(infoM2.treatmentCounts.order2oa > 0, ...
        'MaxPolyOrder=2 overlap-add must be counted as order2oa, counts were %s', ...
        jsonencode(infoM2.treatmentCounts));
    assert(infoM2.treatmentCounts.order3oa == 0, ...
        'an order-2 overlap-add must NOT be reported as order3oa, counts were %s', ...
        jsonencode(infoM2.treatmentCounts));

    % ---- 4b. ladder selection AT the exact documented rung boundaries ----
    % The shipped case table above (3.000/1.800/1.400/1.000/0.500 s) sits safely inside
    % each rung, so it never exercises an edge. C(n)*n/(2*T) and MaxCutoffHz are
    % independently-rounded floating-point values that need not cancel exactly even when
    % T is chosen to land exactly on a documented boundary - at T = 1.635 s, order 3, the
    % raw comparison evaluates to 1.0 + 1 ULP and fails a strict <=, which would silently
    % select order 2 at a duration the docstring documents as order 3. This block pins
    % every rung's exact edge from the docstring table:
    %   T = 2.000 s  -> order 3, single window (not overlap-add: wsize is 2.000 s too,
    %                   and overlap-add requires numel(x) > wsize strictly)
    %   T = 1.635 s  -> order 3, single window (THE case that was failing)
    %   T = 1.240 s  -> order 2, single window
    %   T = 0.845 s  -> order 1, single window
    %
    % The record length is sized to an EXACT multiple of segLen (rather than reusing the
    % fixed n=10000 s record and tiling steps into it), so every sub-segment - including
    % the last one - is exactly segLen samples long. A ragged tail sub-segment left over
    % from an inexact tiling would land in whatever rung ITS OWN duration implies, and a
    % mere ">0" count on the target rung cannot tell that apart from the boundary segments
    % actually landing there - it only requires assert()ing an exact count, and that every
    % OTHER rung's count is exactly zero, that the test can no longer be fooled by a tail.
    edgeCases = { 2.000, 'order3'; 1.635, 'order3'; 1.240, 'order2'; 0.845, 'order1' };
    tcFieldsAll = {'order3oa','order2oa','order1oa','order3','order2','order1','meanOnly','fallback'};
    nSubSegsEdge = 5;
    for k = 1:size(edgeCases,1)
        segLen       = round(edgeCases{k,1} * srate);
        nCaseSamples = nSubSegsEdge * segLen;              % exact multiple -> no tail
        stepSamplesE = (1:nSubSegsEdge-1) * segLen + 1;
        [~, infoE] = RELAX_repair_steps(mkEEG(randn(1, nCaseSamples) * 5), ...
            mkSteps(ones(1, numel(stepSamplesE)), stepSamplesE));
        assert(infoE.treatmentCounts.(edgeCases{k,2}) == nSubSegsEdge, ...
            'a %.3f s sub-segment (exact rung boundary, no ragged tail) should use rung %s exactly %d times, counts were %s', ...
            edgeCases{k,1}, edgeCases{k,2}, nSubSegsEdge, jsonencode(infoE.treatmentCounts));
        for f = 1:numel(tcFieldsAll)
            if ~strcmp(tcFieldsAll{f}, edgeCases{k,2})
                assert(infoE.treatmentCounts.(tcFieldsAll{f}) == 0, ...
                    'a %.3f s sub-segment should land ONLY on rung %s, but %s was %d (counts were %s)', ...
                    edgeCases{k,1}, edgeCases{k,2}, tcFieldsAll{f}, ...
                    infoE.treatmentCounts.(tcFieldsAll{f}), jsonencode(infoE.treatmentCounts));
            end
        end
    end

    % ---- 5. NaN in a sub-segment takes the mean-only path without propagating ----
    xn = randn(1, n) * 5;
    xn(5000:5010) = NaN;
    [outN, infoN] = RELAX_repair_steps(mkEEG(xn), mkSteps([], []));
    assert(~all(isnan(outN.data(1,:))), 'NaN must not propagate across the whole channel');
    assert(sum(isnan(outN.data(1,:))) == 11, 'only the original 11 NaN samples remain NaN');
    assert(infoN.treatmentCounts.meanOnly > 0 || infoN.treatmentCounts.fallback > 0, ...
        'a NaN sub-segment must be counted on a degraded rung');

    % ---- 6. a constant sub-segment cannot produce NaN or an unhandled error ----
    xc = randn(1, n) * 5;
    xc(2001:4000) = 42;                          % perfectly flat stretch
    [outC, ~] = RELAX_repair_steps(mkEEG(xc), mkSteps(1, [2001 4001]));
    assert(~any(isnan(outC.data(1,:))), 'a constant sub-segment must not yield NaN');

    % ---- 7. steps on one channel do not partition another ----
    x2 = [randn(1,n)*5; randn(1,n)*5];
    [~, info2] = RELAX_repair_steps(mkEEG(x2), mkSteps(1, [3001 6001]));
    assert(info2.subSegmentsPerChannel(1) == 3, 'channel 1 has 2 steps -> 3 sub-segments');
    assert(info2.subSegmentsPerChannel(2) == 1, 'channel 2 has no steps -> 1 sub-segment');
    % Orientation is load-bearing: per-channel info fields must be columns, pinned
    % explicitly rather than via isequal (which can pass for a wrong shape when both
    % sides happen to be empty).
    assert(size(info2.subSegmentsPerChannel,1) == 2 && size(info2.subSegmentsPerChannel,2) == 1, ...
        'subSegmentsPerChannel must be nbchan x 1, got %dx%d', size(info2.subSegmentsPerChannel));
    assert(size(info2.lowFreqVarianceRemoved,1) == 2 && size(info2.lowFreqVarianceRemoved,2) == 1, ...
        'lowFreqVarianceRemoved must be nbchan x 1, got %dx%d', size(info2.lowFreqVarianceRemoved));

    % ---- 8. covariates are populated ----
    [~, info8] = RELAX_repair_steps(mkEEG(drift), mkSteps([], []));
    assert(info8.lowFreqVarianceRemoved(1) > 0, 'variance removed should be positive on drift');
    assert(info8.lowFreqVarianceRemoved(1) <= 1, 'variance removed is a proportion');
    assert(info8.fractionSamplesMeanOnly >= 0 && info8.fractionSamplesMeanOnly <= 1, ...
        'fractionSamplesMeanOnly is a proportion');

    % ---- 9. the data array keeps its class and shape ----
    EEGs = mkEEG(single(randn(2, n) * 5));
    [outS, ~] = RELAX_repair_steps(EEGs, mkSteps([], []));
    assert(isequal(size(outS.data), [2 n]), 'shape must be preserved');

    % ---- 10. degrade-safely: nt_detrend missing must not error, and must still return
    % a fully-formed info struct (every documented field present, correctly shaped) ----
    vendorDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'vendor', 'noisetools');
    rmpath(vendorDir);
    cleanupPath = onCleanup(@() RELAX_add_vendored_noisetools());
    assert(isempty(which('nt_detrend')), 'test setup: nt_detrend must be off the path here');
    x10 = randn(2, n) * 5;
    [out10, info10] = RELAX_repair_steps(mkEEG(x10), mkSteps([], []));
    assert(~info10.available, 'repair must report unavailable when nt_detrend is missing');
    assert(ischar(info10.status) && ~isempty(info10.status), 'status must explain why');
    assert(isequal(out10.data, x10), 'data must be returned unmodified on the failure path');
    assert(size(info10.subSegmentsPerChannel,1) == 2 && size(info10.subSegmentsPerChannel,2) == 1, ...
        'subSegmentsPerChannel must stay nbchan x 1 even on the failure path');
    assert(size(info10.lowFreqVarianceRemoved,1) == 2 && size(info10.lowFreqVarianceRemoved,2) == 1, ...
        'lowFreqVarianceRemoved must stay nbchan x 1 even on the failure path');
    tcFields = {'order3oa','order2oa','order1oa','order3','order2','order1','meanOnly','fallback'};
    for k = 1:numel(tcFields)
        assert(isfield(info10.treatmentCounts, tcFields{k}) && info10.treatmentCounts.(tcFields{k}) == 0, ...
            'treatmentCounts.%s must be present and zero on the failure path', tcFields{k});
    end
    clear cleanupPath   % restores the vendored path before later tests / other test files run

    fprintf('test_RELAX_repair_steps: PASS\n');
end
