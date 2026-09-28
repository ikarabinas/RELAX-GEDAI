function test_RELAX_bridged_channels()
% Bridged-channel detection orchestration (filtering, boundary-respecting epoching,
% eBridge invocation, and QC reporting), using a test-double eBridge so the tests do
% not depend on the real EEGLAB ebridge plugin being installed.

    test_eBridge_not_found_returns_unavailable();
    test_too_few_epochs_returns_unavailable();
    test_no_ed_cutoff_returns_available_but_clean();
    test_bridged_count_zero_returns_available_but_clean();
    test_simple_pair_is_rejected_no_qc_flags();
    test_oversized_and_excessive_cluster_flags_both();

    fprintf('test_RELAX_bridged_channels: PASS\n');
end

function test_eBridge_not_found_returns_unavailable()
    % No test double on the path here - this exercises the real eBridge-missing branch.
    % If the real EEGLAB ebridge plugin happens to be installed and on the path in this
    % environment, that branch is unreachable this way; skip rather than fail, since a
    % plugin being installed is not itself a defect.
    if ~isempty(which('eBridge'))
        fprintf('  (skipping eBridge-not-found case: a real eBridge is on the path)\n');
        return;
    end
    EEG = makeEEG(4, 2, 10);
    warnState = warning('off', 'RELAX_bridged_channels:eBridgeNotFound');
    cleanupWarn = onCleanup(@() warning(warnState));
    [rejectLabels, info] = RELAX_bridged_channels(EEG);
    assert(isempty(rejectLabels), 'eBridge missing must reject nothing');
    assert(~info.available, 'eBridge missing must report unavailable');
    assert(contains(info.status, 'not found'), 'status must explain why');
end

function test_too_few_epochs_returns_unavailable()
    stubsDir = fullfile(fileparts(mfilename('fullpath')), 'stubs');
    addpath(stubsDir);
    cleanupPath = onCleanup(@() rmpath(stubsDir));

    EEG = makeEEG(4, 2, 10); % 20 samples / 10-sample (1s) epochs = 2 epochs, needs 30
    warnState = warning('off', 'RELAX_bridged_channels:tooFewEpochs');
    cleanupWarn = onCleanup(@() warning(warnState));
    [rejectLabels, info] = RELAX_bridged_channels(EEG, 'EpochSec', 1);
    assert(isempty(rejectLabels), 'too few epochs must reject nothing');
    assert(~info.available, 'too few epochs must report unavailable');
    assert(contains(info.status, 'too few epochs'), 'status must explain why');
end

function test_no_ed_cutoff_returns_available_but_clean()
    stubsDir = fullfile(fileparts(mfilename('fullpath')), 'stubs');
    addpath(stubsDir);
    cleanupPath = onCleanup(@() rmpath(stubsDir));

    EEG = makeEEG(40, 35, 10); % 350 samples / 10-sample epochs = 35 epochs >= 30
    [rejectLabels, info] = RELAX_bridged_channels(EEG, 'EpochSec', 1);
    assert(isempty(rejectLabels), 'no ED cutoff must reject nothing');
    assert(info.available, 'a found-but-empty cutoff is still an available, clean result');
    assert(info.numEpochs == 35, 'numEpochs is reported even when nothing is rejected');
end

function test_bridged_count_zero_returns_available_but_clean()
    stubsDir = fullfile(fileparts(mfilename('fullpath')), 'stubs');
    addpath(stubsDir);
    cleanupPath = onCleanup(@() rmpath(stubsDir));

    EEG = makeEEG(41, 35, 10);
    [rejectLabels, info] = RELAX_bridged_channels(EEG, 'EpochSec', 1);
    assert(isempty(rejectLabels), 'Bridged.Count==0 must reject nothing');
    assert(info.available && info.nBridgedChannels == 0, 'available and clean');
end

function test_simple_pair_is_rejected_no_qc_flags()
    stubsDir = fullfile(fileparts(mfilename('fullpath')), 'stubs');
    addpath(stubsDir);
    cleanupPath = onCleanup(@() rmpath(stubsDir));

    EEG = makeEEG(30, 35, 10); % E1,E2 bridged; 2/30 = 6.7% is under the 10% QC cap
    [rejectLabels, info] = RELAX_bridged_channels(EEG, 'EpochSec', 1);
    assert(isequal(rejectLabels, {'E1';'E2'}), 'both members of the flagged pair are rejected, in channel order');
    assert(info.available && info.nBridgedChannels == 2 && info.nBridgedPairs == 1, 'counts reported');
    assert(~info.oversizedCluster && ~info.excessiveBridging, 'a small, proportionate cluster raises no QC flag');
end

function test_oversized_and_excessive_cluster_flags_both()
    stubsDir = fullfile(fileparts(mfilename('fullpath')), 'stubs');
    addpath(stubsDir);
    cleanupPath = onCleanup(@() rmpath(stubsDir));

    EEG = makeEEG(12, 35, 10); % E1-E2-E3-E4-E5 chain: 5 > MaxClusterSize(4), 5/12=41.7% > 10%
    warnState1 = warning('off', 'RELAX_bridged_channels:oversizedCluster');
    warnState2 = warning('off', 'RELAX_bridged_channels:excessiveBridging');
    cleanupWarn = onCleanup(@() restoreWarnings(warnState1, warnState2));
    [rejectLabels, info] = RELAX_bridged_channels(EEG, 'EpochSec', 1);
    assert(isequal(rejectLabels, {'E1';'E2';'E3';'E4';'E5'}), 'every chained member is rejected');
    assert(info.oversizedCluster, 'a 5-member cluster must exceed the default MaxClusterSize of 4');
    assert(info.excessiveBridging, '5 of 12 channels must exceed the default MaxProportionChannels of 0.10');
end

function restoreWarnings(warnState1, warnState2)
    warning(warnState1);
    warning(warnState2);
end

function EEG = makeEEG(nChan, nSeconds, srate)
    nSamples = round(nSeconds * srate);
    EEG.data = randn(nChan, nSamples);
    EEG.srate = srate;
    EEG.nbchan = nChan;
    EEG.pnts = nSamples;
    EEG.event = [];
    labels = arrayfun(@(k) sprintf('E%d', k), 1:nChan, 'UniformOutput', false);
    EEG.chanlocs = struct('labels', labels);
end
