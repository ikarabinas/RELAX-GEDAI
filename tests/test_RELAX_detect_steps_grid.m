function test_RELAX_detect_steps_grid()
% 'ThresholdPairs' reuses one pass of the per-channel statistics for several (ZThreshold,
% ChannelZThreshold) pairs. Each pair's result must be exactly what a separate single-pair call
% returns, and running the channels on an open parallel pool must not change any result.

    srate = 1000;
    n     = 20000;          % 20 s
    rng(11);

    function EEG = mkEEG(data, boundaries)
        if nargin < 2; boundaries = []; end
        EEG = struct();
        EEG.data     = data;
        EEG.nbchan   = size(data,1);
        EEG.pnts     = size(data,2);
        EEG.srate    = srate;
        EEG.trials   = 1;
        EEG.event    = struct('type', {}, 'latency', {});
        for b = boundaries(:)'
            EEG.event(end+1) = struct('type', 'boundary', 'latency', b);
        end
        EEG.chanlocs = struct('labels', ...
            arrayfun(@(k) sprintf('E%d', k), 1:size(data,1), 'UniformOutput', false));
    end

    % six channels: pops of several sizes and rise times, a wave, a random walk, a common-mode
    % step on three channels, one constant (degenerate) channel, and two boundaries
    t = 0:n-1;
    X = randn(6, n) * 5;
    X(1,:) = X(1,:) + 200 * (t >= 3000) + 40 * (t >= 12000);
    X(2,:) = X(2,:) + 100 ./ (1 + exp(-(t - 5000) / 10)) + 60 * sin(2*pi*2*t/srate) .* (t > 14000);
    X(3,:) = X(3,:) + cumsum(randn(1, n)) + 30 * (t >= 9000);
    X(1:3,:) = X(1:3,:) + 80 * (t >= 16000);
    X(4,:) = X(4,:) + 25 * (t >= 7000) - 25 * (t >= 7300);
    X(5,:) = 3;
    X(6,:) = X(6,:) + 150 * (t >= 2500) .* (1 - 0.6 * exp(-max(t - 2500, 0) / 40));
    EEG = mkEEG(X, [8000 15000]);

    pairs = [4 4; 12 10; 20 35; 50 10];
    optionSets = {{}, ...
                  {'TransitionMs', 1, 'WindowMs', 51}, ...
                  {'TransitionMs', 10, 'RefineTransition', false}};

    %% ---- 1. every pair equals its own single-pair call ----
    for o = 1:numel(optionSets)
        opts = optionSets{o};
        [gSteps, gInfo] = RELAX_detect_steps(EEG, opts{:}, 'ThresholdPairs', pairs);
        assert(numel(gSteps) == size(pairs, 1) && numel(gInfo) == size(pairs, 1), ...
            'options %d: one steps table and one info per threshold pair', o);
        for k = 1:size(pairs, 1)
            [sSteps, sInfo] = RELAX_detect_steps(EEG, opts{:}, ...
                'ZThreshold', pairs(k,1), 'ChannelZThreshold', pairs(k,2));
            assert(isequaln(gSteps(k), sSteps), ...
                'options %d, pair %d/%d: steps differ from the single-pair call', o, pairs(k,:));
            assertInfoEqual(gInfo(k), sInfo, sprintf('options %d, pair %d/%d', o, pairs(k,:)));
            assert(gInfo(k).params.ZThreshold == pairs(k,1) && ...
                   gInfo(k).params.ChannelZThreshold == pairs(k,2), ...
                'options %d, pair %d/%d: info.params must carry the pair''s own thresholds', ...
                o, pairs(k,:));
        end
        nFound = arrayfun(@(s) numel(s.sample), gSteps);
        assert(nFound(1) > 0, 'options %d: the loosest pair must find steps (fixture check)', o);
    end

    %% ---- 2. a pair order is kept, and one pair is the single-pair call ----
    [gSteps, ~] = RELAX_detect_steps(EEG, 'ThresholdPairs', flipud(pairs));
    [sSteps, ~] = RELAX_detect_steps(EEG, 'ZThreshold', pairs(1,1), 'ChannelZThreshold', pairs(1,2));
    assert(isequaln(gSteps(end), sSteps), 'results follow the order of ThresholdPairs rows');
    [defSteps, defInfo] = RELAX_detect_steps(EEG);
    defaultPair = [defInfo.params.ZThreshold defInfo.params.ChannelZThreshold];
    [oneSteps, oneInfo] = RELAX_detect_steps(EEG, 'ThresholdPairs', defaultPair);
    assert(isequaln(oneSteps, defSteps), 'a single pair equal to the defaults gives the default result');
    assertInfoEqual(oneInfo, defInfo, 'single pair');

    %% ---- 3. a segment too short for any window: every pair reports unavailable ----
    [gSteps, gInfo] = RELAX_detect_steps(mkEEG(randn(2, 300)), 'ThresholdPairs', pairs);
    assert(numel(gInfo) == size(pairs, 1) && ~any([gInfo.available]), ...
        'too-short data: every pair unavailable');
    assert(all(arrayfun(@(s) isempty(s.sample), gSteps)), 'too-short data: no steps for any pair');

    %% ---- 4. malformed ThresholdPairs is an error, not a silent default ----
    for bad = {[12 10 5], [12; 10], zeros(0, 2), [12 NaN]}
        threw = false;
        try
            RELAX_detect_steps(EEG, 'ThresholdPairs', bad{1});
        catch
            threw = true;
        end
        assert(threw, 'ThresholdPairs of size %s must be rejected', mat2str(size(bad{1})));
    end

    %% ---- 5. an open parallel pool changes nothing ----
    % The detector never opens a pool itself; with one open it spreads channels over it.
    hadPool = ~isempty(gcp('nocreate'));
    [refSteps, refInfo] = RELAX_detect_steps(EEG, 'ThresholdPairs', pairs);
    if hadPool
        fprintf(['test_RELAX_detect_steps_grid: a pool was already open, so the serial ' ...
                 'reference could not be taken; section 5 checks run-to-run equality only\n']);
    else
        assert(isempty(gcp('nocreate')), 'RELAX_detect_steps must not open a parallel pool');
        c = parcluster('local');
        c.JobStorageLocation = tempname;
        mkdir(c.JobStorageLocation);
        pool = parpool(c, 2);
        cleanup = onCleanup(@() delete(pool)); %#ok<NASGU>
    end
    [parSteps, parInfo] = RELAX_detect_steps(EEG, 'ThresholdPairs', pairs);
    for k = 1:size(pairs, 1)
        assert(isequaln(parSteps(k), refSteps(k)), 'pair %d/%d: steps differ with a pool open', pairs(k,:));
        assertInfoEqual(parInfo(k), refInfo(k), sprintf('pool, pair %d/%d', pairs(k,:)));
    end
    [parOne, parOneInfo] = RELAX_detect_steps(EEG, 'TransitionMs', 1);
    [serOne, serOneInfo] = RELAX_detect_steps_serial_reference(EEG, hadPool);
    if ~hadPool
        assert(isequaln(parOne, serOne), 'single-pair call: steps differ with a pool open');
        assertInfoEqual(parOneInfo, serOneInfo, 'pool, single-pair call');
    end
end

function [steps, info] = RELAX_detect_steps_serial_reference(EEG, hadPool)
% The single-pair call with no pool open: close the test's pool, run, and leave it closed
% (onCleanup in the caller tolerates an already-deleted pool).
    if hadPool; steps = []; info = []; return; end
    delete(gcp('nocreate'));
    [steps, info] = RELAX_detect_steps(EEG, 'TransitionMs', 1);
end

function assertInfoEqual(a, b, where)
% Every info field equal, except params.ThresholdPairs (the only option the two calls differ in).
    fa = sort(fieldnames(a)); fb = sort(fieldnames(b));
    assert(isequal(fa, fb), '%s: info field names differ', where);
    for f = fa'
        va = a.(f{1}); vb = b.(f{1});
        if strcmp(f{1}, 'params')
            if isfield(va, 'ThresholdPairs'); va = rmfield(va, 'ThresholdPairs'); end
            if isfield(vb, 'ThresholdPairs'); vb = rmfield(vb, 'ThresholdPairs'); end
            va = rmfield(va, {'ZThreshold', 'ChannelZThreshold'});
            vb = rmfield(vb, {'ZThreshold', 'ChannelZThreshold'});
        end
        assert(isequaln(va, vb), '%s: info.%s differs', where, f{1});
    end
end
