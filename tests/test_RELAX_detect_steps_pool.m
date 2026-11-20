function test_RELAX_detect_steps_pool()
% RELAX_detect_steps spreads channels over an already-open parallel pool (the slurm wrapper
% has one) and never opens one itself. With or without a pool the result must be identical.

    srate = 1000;
    n     = 20000;          % 20 s
    rng(11);

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
    EEG = makeEEG(X, [8000 15000], srate);

    optionSets = {{'ZThreshold', 4, 'ChannelZThreshold', 4}, ...
                  {'TransitionMs', 1, 'WindowMs', 51}, ...
                  {'TransitionMs', 10, 'RefineTransition', false}};

    % Earlier tests can leave an auto-created pool open (a parfor elsewhere in the suite);
    % the serial reference needs none.
    delete(gcp('nocreate'));

    serial = cell(size(optionSets));
    for o = 1:numel(optionSets)
        [serial{o}.steps, serial{o}.info] = RELAX_detect_steps(EEG, optionSets{o}{:});
    end
    assert(isempty(gcp('nocreate')), 'RELAX_detect_steps must not open a parallel pool');
    assert(~isempty(serial{1}.steps.sample), 'the loose thresholds must find steps (fixture check)');

    c = parcluster('local');
    c.JobStorageLocation = tempname;
    mkdir(c.JobStorageLocation);
    pool = parpool(c, 2);
    cleanup = onCleanup(@() delete(pool)); %#ok<NASGU>

    for o = 1:numel(optionSets)
        [parSteps, parInfo] = RELAX_detect_steps(EEG, optionSets{o}{:});
        assert(isequaln(parSteps, serial{o}.steps), 'options %d: steps differ with a pool open', o);
        assert(isequaln(parInfo, serial{o}.info), 'options %d: info differs with a pool open', o);
    end

    fprintf('test_RELAX_detect_steps_pool: PASS\n');
end

function EEG = makeEEG(data, boundaries, srate)
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
