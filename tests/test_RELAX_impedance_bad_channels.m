function test_RELAX_impedance_bad_channels()
% Impedance-based channel screening.

    mk = @(labels, values) makeChanlocs(labels, values);

    % --- Basic threshold behaviour, and that the comparison is STRICTLY greater-than ---
    chanlocs = mk({'E1','E2','E3','E4'}, [10, 500, 500.1, 900]);
    [labels, info] = RELAX_impedance_bad_channels(chanlocs, 500, 0.5);
    assert(isequal(labels, {'E3';'E4'}), 'only channels strictly over the limit are rejected');
    assert(info.available, 'impedances were present so the screen is available');
    assert(info.nChannels == 4 && info.nWithImpedance == 4, 'channel counts');
    assert(info.nOverThreshold == 2 && info.nRejected == 2, 'rejection counts');
    assert(~info.capApplied, 'cap must not fire below the limit');
    assert(isequal(info.rejectImpedances, [500.1, 900]), 'rejected impedances reported');

    % Labels come back as a column, so pop_select and the reporting table agree
    assert(iscell(labels) && size(labels,2) == 1, 'labels must be a column cell array');

    % --- Nothing over the limit -> empty, but still available ---
    chanlocs = mk({'E1','E2'}, [10, 20]);
    [labels, info] = RELAX_impedance_bad_channels(chanlocs, 500, 0.5);
    assert(isempty(labels), 'nothing over the limit rejects nothing');
    assert(info.available && info.nOverThreshold == 0 && info.nRejected == 0, 'clean file counts');
    assert(info.medianImpedance == 15 && info.maxImpedance == 20, 'summary stats over known values');

    % --- NaN (unknown impedance) is NOT bad ---
    chanlocs = mk({'E1','E2','E3'}, [NaN, 10, 900]);
    [labels, info] = RELAX_impedance_bad_channels(chanlocs, 500, 0.9);
    assert(isequal(labels, {'E3'}), 'NaN impedance must never be rejected');
    assert(info.nWithImpedance == 2, 'NaN channels are not counted as known');
    assert(info.medianImpedance == 455, 'summary stats ignore NaN');

    % --- All-NaN impedances (real files in this cohort look like this) ---
    chanlocs = mk({'E1','E2','E3'}, [NaN, NaN, NaN]);
    warnState = warning('off', 'RELAX_impedance_bad_channels:noImpedanceValues');
    cleanupWarn = onCleanup(@() warning(warnState));
    [labels, info] = RELAX_impedance_bad_channels(chanlocs, 500, 0.5);
    assert(isempty(labels), 'all-NaN impedances reject nothing');
    assert(~info.available, 'all-NaN must report itself unavailable, not as a clean file');
    assert(info.nWithImpedance == 0, 'no known impedances');
    clear cleanupWarn;

    % --- No impedance field at all ---
    chanlocs = struct('labels', {'E1','E2'});
    warnState = warning('off', 'RELAX_impedance_bad_channels:noImpedanceField');
    cleanupWarn = onCleanup(@() warning(warnState));
    [labels, info] = RELAX_impedance_bad_channels(chanlocs, 500, 0.5);
    assert(isempty(labels), 'a missing impedance field rejects nothing');
    assert(~info.available, 'a missing impedance field is unavailable');
    clear cleanupWarn;

    % --- Empty chanlocs must not error ---
    [labels, info] = RELAX_impedance_bad_channels(struct('labels', {}, 'impedance', {}), 500, 0.5);
    assert(isempty(labels) && ~info.available && info.nChannels == 0, 'empty chanlocs tolerated');

    % --- The cap keeps the WORST channels and reports itself ---
    % 10 channels, 5 over the limit, cap of 0.25 -> floor(2.5) = 2 rejections allowed.
    chanlocs = mk({'E1','E2','E3','E4','E5','E6','E7','E8','E9','E10'}, ...
                  [10, 600, 20, 5000, 30, 700, 40, 900, 50, 4000]);
    warnState = warning('off', 'RELAX_impedance_bad_channels:tooManyOverThreshold');
    cleanupWarn = onCleanup(@() warning(warnState));
    [labels, info] = RELAX_impedance_bad_channels(chanlocs, 500, 0.25);
    assert(info.capApplied, 'cap must fire when more channels are over the limit than allowed');
    assert(info.nOverThreshold == 5, 'the over-threshold count reports the true total');
    assert(info.nRejected == 2, 'only the capped number is rejected');
    % E4 (5000) and E10 (4000) are the two worst; labels stay in channel order
    assert(isequal(labels, {'E4';'E10'}), 'the cap keeps the highest-impedance channels');
    clear cleanupWarn;

    % --- A cap that rounds to zero rejects nothing rather than erroring ---
    chanlocs = mk({'E1','E2'}, [900, 900]);
    warnState = warning('off', 'RELAX_impedance_bad_channels:tooManyOverThreshold');
    cleanupWarn = onCleanup(@() warning(warnState));
    [labels, info] = RELAX_impedance_bad_channels(chanlocs, 500, 0.25);
    assert(isempty(labels) && info.nRejected == 0, 'a cap of floor(0.5)=0 rejects nothing');
    assert(info.capApplied && info.nOverThreshold == 2, 'the cap is still reported');
    clear cleanupWarn;

    % --- Malformed impedance values degrade to NaN rather than erroring ---
    chanlocs = mk({'E1','E2','E3'}, [10, 10, 10]);
    chanlocs(1).impedance = [];          % empty
    chanlocs(2).impedance = 'high';      % char
    chanlocs(3).impedance = 900;
    [labels, info] = RELAX_impedance_bad_channels(chanlocs, 500, 0.9);
    assert(isequal(labels, {'E3'}), 'malformed impedances are treated as unknown, not bad');
    assert(info.nWithImpedance == 1, 'only the well-formed value counts as known');

    fprintf('test_RELAX_impedance_bad_channels: PASS\n');
end

function chanlocs = makeChanlocs(labels, values)
    chanlocs = struct('labels', labels, 'impedance', num2cell(values));
end
