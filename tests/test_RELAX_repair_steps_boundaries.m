function test_RELAX_repair_steps_boundaries()
% Repair never fits across a boundary: each segment is processed as if standalone.

    RELAX_add_vendored_noisetools();

    srate = 1000;
    n     = 12000;
    rng(23);

    function EEG = mkEEG(data, boundaryLatencies)
        EEG = struct();
        EEG.data     = data;
        EEG.nbchan   = size(data,1);
        EEG.pnts     = size(data,2);
        EEG.srate    = srate;
        EEG.trials   = 1;
        if isempty(boundaryLatencies)
            EEG.event = struct('type', {}, 'latency', {});
        else
            EEG.event = struct('type', repmat({'boundary'}, 1, numel(boundaryLatencies)), ...
                               'latency', num2cell(boundaryLatencies));
        end
        EEG.chanlocs = struct('labels', ...
            arrayfun(@(k) sprintf('E%d', k), 1:size(data,1), 'UniformOutput', false));
    end

    noSteps = struct('channel', zeros(0,1), 'label', {cell(0,1)}, 'sample', zeros(0,1), ...
                     'amplitude', zeros(0,1), 'z', zeros(0,1));

    % Two segments with very different drift slopes. If the fit ran across the boundary,
    % neither segment would come out flat.
    x = randn(1, n) * 5;
    x(1:6000)     = x(1:6000)     + linspace(0, 400, 6000);
    x(6001:12000) = x(6001:12000) + linspace(0, -900, 6000);

    EEG = mkEEG(x, 6000.5);
    [out, info] = RELAX_repair_steps(EEG, noSteps);

    seg1 = out.data(1, 1:6000);
    seg2 = out.data(1, 6001:12000);
    assert(abs(median(seg1(5001:6000)) - median(seg1(1:1000))) < 20, ...
        'segment 1 must be flattened on its own');
    assert(abs(median(seg2(5001:6000)) - median(seg2(1:1000))) < 20, ...
        'segment 2 must be flattened on its own');

    % Processing each half as a standalone record must give a bit-identical answer.
    [outA, ~] = RELAX_repair_steps(mkEEG(x(1:6000), []), noSteps);
    [outB, ~] = RELAX_repair_steps(mkEEG(x(6001:end), []), noSteps);
    assert(isequal(out.data(1,1:6000),   outA.data(1,:)), ...
        'segment 1 must be bit-identical to processing it alone');
    assert(isequal(out.data(1,6001:end), outB.data(1,:)), ...
        'segment 2 must be bit-identical to processing it alone');

    % Sub-segment count: 2 boundary segments, 1 step in the second -> 3 sub-segments
    stepsOne = struct('channel', 1, 'label', {{'E1'}}, 'sample', 9000, ...
                      'amplitude', 0, 'z', 0);
    [~, info3] = RELAX_repair_steps(mkEEG(x, 6000.5), stepsOne);
    assert(info3.subSegmentsPerChannel(1) == 3, ...
        'boundary + 1 step -> 3 sub-segments, got %d', info3.subSegmentsPerChannel(1));

    % A step landing on the boundary itself must not create a zero-length sub-segment
    stepsAtBoundary = struct('channel', 1, 'label', {{'E1'}}, 'sample', 6001, ...
                             'amplitude', 0, 'z', 0);
    [~, infoB] = RELAX_repair_steps(mkEEG(x, 6000.5), stepsAtBoundary);
    assert(infoB.subSegmentsPerChannel(1) == 2, ...
        'a step coinciding with a boundary must merge into it, got %d sub-segments', ...
        infoB.subSegmentsPerChannel(1));
    assert(info.available, 'sanity: repair was available throughout');

    % Orientation is load-bearing: pin the shape explicitly rather than via isequal.
    assert(size(info.subSegmentsPerChannel,1) == 1 && size(info.subSegmentsPerChannel,2) == 1, ...
        'subSegmentsPerChannel must be nbchan x 1, got %dx%d', size(info.subSegmentsPerChannel));

    fprintf('test_RELAX_repair_steps_boundaries: PASS\n');
end
