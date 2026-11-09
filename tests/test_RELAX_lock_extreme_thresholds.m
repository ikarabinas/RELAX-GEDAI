function test_RELAX_lock_extreme_thresholds()
% The second extreme-marking pass judges cleaned epochs by the first pass's per-channel
% statistics. Fitted and applied back to the same epochs, they must reproduce what
% pop_jointprob and pop_rejkurt computed there, and applying them to other epochs must
% never re-fit them.

    [EEG, labels] = epoched_fixture();
    cfg = struct('ExtremeVoltageShiftThreshold', 25);
    EEG = pop_jointprob(EEG, 1, 1:8, 10, 1000, 0, 0, 0);
    EEG = pop_rejkurt(EEG, 1, 1:8, 10, 1e6, 0, 0);
    locked = RELAX_lock_extreme_thresholds(EEG.data, labels, cfg);
    s = RELAX_score_locked_extremes(locked, EEG.data, labels);

    % Reproduction
    assert(max(abs(s.jpZ - EEG.stats.jpE), [], 'all') < 1e-9, 'jp z does not reproduce pop_jointprob');
    assert(max(abs(s.kurtZ - EEG.stats.kurtE), [], 'all') < 1e-9, 'kurtosis z does not reproduce pop_rejkurt');
    r = squeeze(range(double(EEG.data), 2));
    assert(isequal(locked.shiftUpper, median(r, 2) + 25 * mad(r, 1, 2)), 'shift bound is not median + 25 MAD');
    assert(isequal(s.range, r), 'range is not the per-epoch voltage range');
    assert(isequal(s.maxAbs, squeeze(max(abs(double(EEG.data)), [], 2))), 'maxAbs is not the per-epoch max |v|');
    assert(s.kurtZ(3, 7) > 10, 'the one-sample spike should score extreme kurtosis');
    assert(s.jpZ(6, 20) > 10, 'the 30-sample plateau should score improbable');

    % label_subset: rows are matched by label, not position
    subset = {'E2', 'E5', 'E3'};
    s2 = RELAX_score_locked_extremes(locked, EEG.data([2 5 3], :, :), subset);
    assert(isequal(s2.jpZ, s.jpZ([2 5 3], :)) && isequal(s2.kurtZ, s.kurtZ([2 5 3], :)), ...
        'scores of a channel subset differ from the same channels in the full set');
    assert(isequal(s2.shiftUpper, locked.shiftUpper([2 5 3])), 'shiftUpper not matched by label');

    % Unknown label
    threw = false;
    try
        RELAX_score_locked_extremes(locked, EEG.data(1, :, :), {'E99'});
    catch err
        threw = strcmp(err.identifier, 'RELAX:secondPassUnknownChannel');
    end
    assert(threw, 'a channel the first pass never had did not raise RELAX:secondPassUnknownChannel');

    % Not re-fitted: one epoch alone scores as it did among all of them
    s1 = RELAX_score_locked_extremes(locked, EEG.data(:, :, 7), labels);
    assert(max(abs(s1.jpZ - s.jpZ(:, 7))) < 1e-9 && max(abs(s1.kurtZ - s.kurtZ(:, 7))) < 1e-9, ...
        'scoring a single epoch re-fitted the statistics');

    % out_of_range: values outside the fitted range take the probability floor
    above = repmat(locked.jpHi + 50, 1, 100);
    below = repmat(locked.jpLo - 50, 1, 100);
    sOut = RELAX_score_locked_extremes(locked, cat(3, above, below), labels);
    assert(all(isfinite(sOut.jpZ), 'all'), 'values outside the first pass''s range gave non-finite jp z');
    assert(all(sOut.jpZ > 10, 'all'), 'values outside the first pass''s range should be improbable');

    % single_input
    sSingle = RELAX_score_locked_extremes(locked, single(EEG.data), labels);
    singleDiff = max(abs(sSingle.jpZ - s.jpZ), [], 'all');
    assert(singleDiff < 0.25, 'single-precision data score differently: max |z diff| %.3g', singleDiff);
    assert(isa(sSingle.jpZ, 'double'), 'scores should be double');

    % flat_channel: a channel with no spread in the first pass cannot be judged by joint
    % probability or kurtosis (ShiftLow judges it), so its z is 0, never NaN
    flatData = EEG.data;
    flatData(4, :, :) = 5;
    flatLocked = RELAX_lock_extreme_thresholds(flatData, labels, cfg);
    for later = {flatData, EEG.data}
        sFlat = RELAX_score_locked_extremes(flatLocked, later{1}, labels);
        assert(all(sFlat.jpZ(4, :) == 0) && all(sFlat.kurtZ(4, :) == 0), ...
            'a channel flat across the first pass did not score z = 0');
        assert(all(isfinite(sFlat.jpZ([1:3 5:8], :)), 'all'), 'the flat channel disturbed the others');
    end

    fprintf('test_RELAX_lock_extreme_thresholds: PASS\n');
end

function [EEG, labels] = epoched_fixture()
    rng(1);
    data = 10 * randn(8, 100, 300);
    data(3, 50, 7) = 300;
    data(6, 31:60, 20) = 200 + data(6, 31:60, 20);
    labels = arrayfun(@(c) sprintf('E%d', c), 1:8, 'UniformOutput', false);
    EEG = eeg_emptyset();
    EEG.data = data;
    EEG.nbchan = 8; EEG.pnts = 100; EEG.trials = 300; EEG.srate = 100;
    EEG.xmin = 0; EEG.xmax = 0.99;
    EEG.chanlocs = struct('labels', labels);
    EEG = eeg_checkset(EEG);
end
