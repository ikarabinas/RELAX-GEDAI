# Bridged Channel Rejection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a bridged-electrode-cluster rejection stage to `RELAX_Wrapper_slurm.m`, running immediately after the existing impedance-based channel rejection and before the PREP `findNoisyChannels` block, ported from `detect_bridged_channels`/`remove_bridged_channels` in `~/TMS_response_variability/iTBS_pipeline.m` (lines 470-478 + their called functions).

**Architecture:** Two new standalone functions at repo root, following the `RELAX_impedance_bad_channels.m` pattern (pure function, no side effects on `EEG`, caller does the `pop_select`):
- `RELAX_group_bridged_channels.m` — pure connected-components grouping of eBridge's flagged pairs into clusters. Split out from the detection function specifically so its logic is unit-testable without eBridge installed.
- `RELAX_bridged_channels.m` — orchestrates eBridge: optional pre-filter (via the existing `RELAX_filtbutter`), boundary-respecting epoching (via the existing `RELAX_get_segments`), calls `eBridge`, calls `RELAX_group_bridged_channels`, and reports QC info.

`RELAX_Wrapper_slurm.m` gets a new config-gated section (default off) mirroring the impedance section immediately above it: default the new `RELAX_cfg` fields, call `RELAX_bridged_channels`, `pop_select` out anything flagged, and fold the result into `EEG.RELAX_issues_to_check`.

**Tech Stack:** MATLAB, EEGLAB. External dependency: the EEGLAB `ebridge` plugin's `eBridge()` (NOT vendored — same status as PREP/MWF/FieldTrip). Uses base-MATLAB `graph`/`conncomp` (available since R2015b, confirmed present in this environment).

**Spec:** No separate spec doc; source of truth is `~/TMS_response_variability/iTBS_pipeline.m:470-478` (call site) and its local functions `detect_bridged_channels` (line 1861) and `remove_bridged_channels` (line 2065), reworked to match `RELAX_impedance_bad_channels.m`'s calling convention and this repo's boundary-handling helpers (`RELAX_get_segments.m`, `RELAX_filtbutter.m`).

## Global Constraints

- Match `RELAX_impedance_bad_channels.m`'s style exactly: header comment block (Usage/Inputs/Outputs/rationale), a single exported function, `inputParser`-free where the impedance function is (it isn't parameterized) but `RELAX_bridged_channels` DOES need `inputParser` like `detect_bridged_channels` did — keep the same Name-Value convention.
- New `RELAX_cfg` fields default OFF (absent = 0), like `Perform_impedance_bad_channel_rejection` — this is new behaviour and must not silently turn on for users without the `ebridge` plugin installed.
- Rejection labels returned as a column cell array (`{}` if none), matching `RELAX_impedance_bad_channels`'s output contract, so the wrapper's `pop_select(EEG,'nochannel',...)` call is identical in shape to the impedance one.
- Every member of a bridged cluster is removed (never "keep one") — a bridge is an unrecoverable physical short; this is `remove_bridged_channels`'s explicit design decision and must not be diluted.
- No enforced cap on bridge rejections (unlike the impedance screen's `maxProportion` cap): bridging is a direct measurement from the data, not metadata that could itself be wrong, so `MaxProportionChannels`/`MaxClusterSize` stay QC-only warning flags, exactly as in `detect_bridged_channels`.
- Boundary-safe: epoching must go through `RELAX_get_segments` (not a duplicated re-implementation of `compute_segment_windows`), and any optional pre-filter must go through `RELAX_filtbutter`'s existing `respect_boundaries` argument, not a bespoke `apply_filter_by_segments`.
- Tests go in `tests/test_*.m`, auto-discovered by `tests/run_all_tests.m` — no registration step needed.
- The wrapper section itself is not unit-tested (per `CLAUDE.md`: "there is no build step... correctness is verified by running the pipeline on real EEG data"); verify it with `check_matlab_code` (static analysis) and a careful manual trace instead.

---

## File Structure

- Create: `RELAX_group_bridged_channels.m` — pure clustering helper (no eBridge dependency).
- Create: `tests/test_RELAX_group_bridged_channels.m` — unit tests for the clustering helper.
- Create: `RELAX_bridged_channels.m` — orchestrating detection function (mirrors `RELAX_impedance_bad_channels.m`).
- Create: `tests/stubs/eBridge.m` — test double for the EEGLAB `ebridge` plugin's `eBridge()`, used only by the next test file.
- Create: `tests/test_RELAX_bridged_channels.m` — unit tests for the orchestrating function, using the fixture above.
- Modify: `RELAX_Wrapper_slurm.m` — new config-gated section between the impedance block (ends ~line 386) and the PREP comment block (~line 388), plus a small addition to the `EEG.RELAX_issues_to_check` block (~lines 914-937).
- Modify: `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` — new `RELAX_cfg` defaults, documented inline next to the impedance block (~line 338).
- Modify: `CLAUDE.md` — mention the new stage in the "Pipeline architecture" stage-4 bullet and add `ebridge` to "External dependencies".

---

## Task 1: `RELAX_group_bridged_channels.m` — clustering helper

**Files:**
- Create: `RELAX_group_bridged_channels.m`
- Test: `tests/test_RELAX_group_bridged_channels.m`

**Interfaces:**
- Produces: `[clusters, pairLabels] = RELAX_group_bridged_channels(labels, pairPairs)` where `labels` is a 1×n cell array of channel labels (order = channel order) and `pairPairs` is eBridge's `EB.Bridged.Pairs` shape: a 2×nPairs cell array of label pairs. Returns `clusters` (1×nClusters cell array, each element a cell array of ≥2 member labels) and `pairLabels` (1×nMatchedPairs cell array of `'A-B'` strings, pairs referencing an unknown label dropped).

- [ ] **Step 1: Write the failing tests**

```matlab
function test_RELAX_group_bridged_channels()
% Connected-components grouping of eBridge-flagged pairs.

    % --- A single simple pair becomes a 2-member cluster ---
    labels = {'E1','E2','E3','E4'};
    pairs = {'E1'; 'E2'}; % 2x1: one pair
    [clusters, pairLabels] = RELAX_group_bridged_channels(labels, pairs);
    assert(numel(clusters) == 1 && isequal(sort(clusters{1}), {'E1','E2'}), ...
        'a single flagged pair must become one 2-member cluster');
    assert(isequal(pairLabels, {'E1-E2'}), 'pair label string uses channel order');

    % --- A chain of pairs merges into one cluster, not several ---
    labels = {'E1','E2','E3','E4','E5'};
    pairs = [{'E1','E2','E3','E4'}; {'E2','E3','E4','E5'}]; % chain E1-E2-E3-E4-E5
    [clusters, ~] = RELAX_group_bridged_channels(labels, pairs);
    assert(numel(clusters) == 1, 'a chain of pairs must merge into a single cluster');
    assert(isequal(sort(clusters{1}), {'E1','E2','E3','E4','E5'}), ...
        'the merged cluster must contain every chained channel');

    % --- Two disjoint pairs stay as two separate clusters ---
    labels = {'E1','E2','E3','E4','E5','E6'};
    pairs = [{'E1','E4'}; {'E2','E5'}]; % E1-E2 and E4-E5, unrelated
    [clusters, ~] = RELAX_group_bridged_channels(labels, pairs);
    assert(numel(clusters) == 2, 'unrelated pairs must not be merged');
    sizes = cellfun(@numel, clusters);
    assert(all(sizes == 2), 'each disjoint pair is its own 2-member cluster');

    % --- A channel with no flagged pair never appears in any cluster ---
    labels = {'E1','E2','E3'};
    pairs = {'E1'; 'E2'};
    [clusters, ~] = RELAX_group_bridged_channels(labels, pairs);
    allMembers = [clusters{:}];
    assert(~ismember('E3', allMembers), 'an unpaired channel must not appear in any cluster');

    % --- A pair naming an unknown label is silently skipped ---
    labels = {'E1','E2'};
    pairs = [{'E1','E1'}; {'E2','E99'}]; % second pair references E99, not in labels
    [clusters, pairLabels] = RELAX_group_bridged_channels(labels, pairs);
    assert(numel(clusters) == 1 && isequal(sort(clusters{1}), {'E1','E2'}), ...
        'the valid pair must still be grouped');
    assert(isequal(pairLabels, {'E1-E2'}), 'the pair naming an unknown label must be dropped from pairLabels');

    % --- No pairs -> no clusters ---
    [clusters, pairLabels] = RELAX_group_bridged_channels({'E1','E2'}, cell(2,0));
    assert(isempty(clusters) && isempty(pairLabels), 'empty pair list produces no clusters');

    % --- Empty label list -> no clusters, no error ---
    [clusters, pairLabels] = RELAX_group_bridged_channels({}, cell(2,0));
    assert(isempty(clusters) && isempty(pairLabels), 'empty channel list tolerated');

    fprintf('test_RELAX_group_bridged_channels: PASS\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_group_bridged_channels"`
Expected: FAIL with "Undefined function 'RELAX_group_bridged_channels'"

- [ ] **Step 3: Write the implementation**

```matlab
% RELAX_group_bridged_channels() - Group eBridge-flagged channel pairs into
%                                   connected clusters of electrically bridged channels.
%
% Usage:
%   >> [clusters, pairLabels] = RELAX_group_bridged_channels(labels, pairPairs);
%
% Inputs:
%   labels     - cell array of channel labels, in EEG.chanlocs order
%   pairPairs  - eBridge's EB.Bridged.Pairs: 2 x nPairs cell array of channel label
%                pairs flagged as bridged. A pair naming a label not present in
%                `labels` (e.g. a channel already removed upstream) is silently
%                dropped, matching the source pipeline's behaviour.
%
% Outputs:
%   clusters   - 1 x nClusters cell array; each element is a cell array of channel
%                labels electrically shorted together (connected components of the
%                flagged-pair graph, so bridge chains of more than 2 channels are
%                grouped rather than split into separate pairs)
%   pairLabels - 1 x nMatchedPairs cell array of 'A-B' strings, one per matched pair
%                (pairs with an unmatched label excluded), for QC reporting
%
% Bridges chain: electrolyte spreading across a patch commonly shorts 3-5 electrodes
% into one cluster, and handling flagged pairs independently would miss that they are
% one measurement, not several unrelated ones. Split out of RELAX_bridged_channels so
% this grouping logic can be unit-tested without the eBridge plugin installed.

function [clusters, pairLabels] = RELAX_group_bridged_channels(labels, pairPairs)

    clusters = {};
    pairLabels = {};

    numChans = numel(labels);
    numPairs = size(pairPairs, 2);
    if numChans == 0 || numPairs == 0
        return;
    end

    adjacency = false(numChans);
    pairLabels = cell(1, numPairs);
    for k = 1:numPairs
        i = find(strcmp(labels, pairPairs{1,k}), 1);
        j = find(strcmp(labels, pairPairs{2,k}), 1);
        if isempty(i) || isempty(j)
            continue;
        end
        adjacency(i,j) = true;
        adjacency(j,i) = true;
        pairLabels{k} = sprintf('%s-%s', labels{i}, labels{j});
    end
    pairLabels = pairLabels(~cellfun('isempty', pairLabels));

    comp = conncomp(graph(adjacency));
    for c = 1:max(comp)
        members = find(comp == c);
        if numel(members) < 2
            continue; % an unbridged channel sitting in its own component
        end
        clusters{end+1} = labels(members); %#ok<AGROW>
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_group_bridged_channels"`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add RELAX_group_bridged_channels.m tests/test_RELAX_group_bridged_channels.m
git commit -m "Add bridged-channel-pair connected-components grouping helper"
```

---

## Task 2: `RELAX_bridged_channels.m` — detection orchestration

**Files:**
- Create: `RELAX_bridged_channels.m`
- Create: `tests/stubs/eBridge.m`
- Test: `tests/test_RELAX_bridged_channels.m`

**Interfaces:**
- Consumes: `RELAX_group_bridged_channels(labels, pairPairs)` from Task 1 (exact signature above); `RELAX_get_segments(EEG)` (existing, returns `[nSegments x 2]` sample ranges, `[1 EEG.pnts]` when no boundaries); `RELAX_filtbutter(EEG, high, low, ord, type, causal_or_acausal_filter, respect_boundaries)` (existing).
- Produces: `[rejectLabels, info] = RELAX_bridged_channels(EEG, varargin)`. `rejectLabels` is a column cell array of channel labels (channel order), `{}` if none/unavailable. `info` is a struct with fields: `available` (logical), `status` (char), `filtered` (logical), `filtBand` (`[low high]`), `epochSec`, `numEpochs`, `EDcutoff`, `EDscale`, `clusters` (cell array of cell arrays), `pairs` (cell array of `'A-B'` strings), `nBridgedChannels`, `nBridgedPairs`, `oversizedCluster` (logical), `excessiveBridging` (logical), `rejectLabels` (same as first output). Consumed by Task 3's wrapper call.

- [ ] **Step 1: Write the test double eBridge.m used only by the test file**

```matlab
function EB = eBridge(EPO, varargin)
% Test double for the EEGLAB ebridge plugin's eBridge(), used only by
% tests/test_RELAX_bridged_channels.m so that RELAX_bridged_channels can be
% exercised end-to-end without the real ebridge plugin installed. Branches on the
% channel count of EPO.data to return a canned result per test scenario - see
% tests/test_RELAX_bridged_channels.m for what each channel count means.

    nChan = size(EPO.data, 1);
    labels = {EPO.chanlocs.labels};

    switch nChan
        case 40 % "no bridging": no near-zero peak in the ED distribution
            EB.Info.EDcutoff = -1;
            EB.Info.EDscale = NaN;
            EB.Bridged.Count = 0;
            EB.Bridged.Pairs = cell(2,0);
        case 41 % "clean file": a cutoff was found but nothing crossed it
            EB.Info.EDcutoff = 0.05;
            EB.Info.EDscale = 1.0;
            EB.Bridged.Count = 0;
            EB.Bridged.Pairs = cell(2,0);
        case 30 % a single 2-channel bridge, well under any QC cap
            EB.Info.EDcutoff = 0.05;
            EB.Info.EDscale = 1.2;
            EB.Bridged.Pairs = {labels{1}; labels{2}};
            EB.Bridged.Count = 1;
        case 12 % a 5-channel chain: triggers both oversizedCluster and excessiveBridging
            EB.Info.EDcutoff = 0.05;
            EB.Info.EDscale = 1.2;
            EB.Bridged.Pairs = [labels(1:4); labels(2:5)];
            EB.Bridged.Count = 4;
        otherwise
            error('eBridge test double: unexpected channel count %d - add a case for it', nChan);
    end
end
```

- [ ] **Step 2: Write the failing tests**

```matlab
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
    % No fixtures on the path here - eBridge really is absent in this environment,
    % which is the state a user who has not installed the ebridge plugin is in.
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
    cleanupWarn = onCleanup(@() [warning(warnState1), warning(warnState2)]);
    [rejectLabels, info] = RELAX_bridged_channels(EEG, 'EpochSec', 1);
    assert(isequal(rejectLabels, {'E1';'E2';'E3';'E4';'E5'}), 'every chained member is rejected');
    assert(info.oversizedCluster, 'a 5-member cluster must exceed the default MaxClusterSize of 4');
    assert(info.excessiveBridging, '5 of 12 channels must exceed the default MaxProportionChannels of 0.10');
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
```

- [ ] **Step 3: Run test to verify it fails**

Run: `matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_bridged_channels"`
Expected: FAIL with "Undefined function 'RELAX_bridged_channels'"

- [ ] **Step 4: Write the implementation**

```matlab
% RELAX_bridged_channels() - Detect electrically bridged electrode clusters via eBridge.
%
% Usage:
%   >> [rejectLabels, info] = RELAX_bridged_channels(EEG);
%   >> [rejectLabels, info] = RELAX_bridged_channels(EEG, 'EpochSec', 1, 'BCT', 0.5);
%
% Inputs:
%   EEG - EEGLAB EEG structure (uses .data, .chanlocs, .srate, .nbchan, .pnts, .event)
%   Name-value pairs (all optional):
%     EpochSec              - epoch length in seconds for the electrical distance
%                              computation (default 1; eBridge's own auto-range is
%                              0.5-1.0 s)
%     Filter                 - 1 = bandpass a working copy before computing ED, via
%                              RELAX_filtbutter; 0 = use the data as-is (default 0).
%                              Off by default: band limiting rescales every ED value,
%                              and this screen normally runs on data already notched.
%     FiltLow, FiltHigh       - bandpass edges in Hz, used only if Filter=1 (default 0.5, 30)
%     FiltOrder               - Butterworth order for the pre-ED bandpass (default 2)
%     BCT                     - eBridge BCT: fraction of epochs below the ED cutoff
%                              required to call a pair bridged (default 0.5)
%     BinSize                 - eBridge BinSize: bin width of the ED frequency
%                              distribution (default 0.1)
%     MaxEpochs                - cap on epochs handed to eBridge, matching eBridge's own
%                              internal cap (default 1500); subsampled deterministically
%     MaxProportionChannels    - QC only, does not change what is rejected: proportion of
%                              channels above which bridging is flagged as excessive
%                              (default 0.10)
%     MaxClusterSize          - QC only: cluster size above which interpolation over the
%                              resulting gap is considered untrustworthy (default 4)
%     RespectBoundaries        - passed through to RELAX_filtbutter when Filter=1 (default 0)
%
% Outputs:
%   rejectLabels - Column cell array of channel labels to reject: every channel
%                  belonging to a bridged cluster (all members, not all-but-one - a
%                  bridge is a physical short and the two original signals are not
%                  recoverable from one measurement; a retained member would still hold
%                  the smeared mixture while its chanlocs entry claims a precise
%                  location). {} if none / eBridge unavailable / too few epochs.
%   info          - struct of QC fields (store on EEG.RELAX for cohort-level review):
%                     .available          false when eBridge is missing, there were too
%                                         few usable epochs, or eBridge itself errored;
%                                         nothing is rejected
%                     .status             human-readable reason when not available
%                     .filtered, .filtBand, .epochSec, .numEpochs
%                     .EDcutoff, .EDscale
%                     .clusters, .pairs, .nBridgedChannels, .nBridgedPairs
%                     .oversizedCluster, .excessiveBridging   QC flags, see above -
%                                         these do not change what is rejected
%                     .rejectLabels
%
% Computes the pairwise electrical distance ED(i,j) = var(x_i - x_j) over short epochs
% and flags pairs whose ED sits in the near-zero peak of the ED distribution (Alschuler
% et al. 2014, Clin Neurophysiol 125:484-490), via the EEGLAB ebridge plugin's eBridge()
% (not vendored in this repo - install the ebridge plugin to use this screen). ED is
% reference-invariant, so it is valid regardless of the current reference. It is also
% amplitude-sensitive, which matters on high-density nets where adjacent electrodes
% routinely correlate at r > 0.95 for purely neurophysiological reasons - a correlation
% criterion would flag them, ED does not.
%
% Epochs are built from RELAX_get_segments so they respect boundary events instead of
% straddling excised (e.g. TMS) artifact windows; a plain reshape would silently corrupt
% ED on any file with boundaries in it. eBridge's own FiltMode is not used regardless (it
% is only honoured for continuous input and silently ignored for epoched input), so
% filtering, if wanted, happens here via RELAX_filtbutter before epoching.
%
% Unlike RELAX_impedance_bad_channels, there is no enforced cap on how many channels this
% screen can reject: bridging is measured directly from the data (a robust statistical
% test already ran), not metadata that could itself be wrong, so there is no argument for
% under-removing a channel this test flags. MaxProportionChannels/MaxClusterSize are
% QC-only: they mark a recording whose reconstruction should not be trusted (nothing
% derived from a bridged channel is fed back in, so a large contiguous cluster leaves a
% gap with no interior information, and spherical interpolation will happily produce a
% smooth, confident, meaningless surface across it) without changing what gets rejected.
%
% Ported from detect_bridged_channels / remove_bridged_channels in
% TMS_response_variability/iTBS_pipeline.m; grouping of flagged pairs into clusters is
% split out to RELAX_group_bridged_channels for independent unit testing.

function [rejectLabels, info] = RELAX_bridged_channels(EEG, varargin)

    p = inputParser();
    p.KeepUnmatched = true;
    addParameter(p, 'EpochSec', 1);
    addParameter(p, 'Filter', 0);
    addParameter(p, 'FiltLow', 0.5);
    addParameter(p, 'FiltHigh', 30);
    addParameter(p, 'FiltOrder', 2);
    addParameter(p, 'BCT', 0.5);
    addParameter(p, 'BinSize', 0.1);
    addParameter(p, 'MaxEpochs', 1500);
    addParameter(p, 'MaxProportionChannels', 0.10);
    addParameter(p, 'MaxClusterSize', 4);
    addParameter(p, 'RespectBoundaries', 0);
    parse(p, varargin{:});
    epochSec          = p.Results.EpochSec;
    doFilter          = p.Results.Filter;
    filtLow           = p.Results.FiltLow;
    filtHigh          = p.Results.FiltHigh;
    filtOrder         = p.Results.FiltOrder;
    bct               = p.Results.BCT;
    binSize           = p.Results.BinSize;
    maxEpochs         = p.Results.MaxEpochs;
    maxPropChans      = p.Results.MaxProportionChannels;
    maxClusterSz      = p.Results.MaxClusterSize;
    respectBoundaries = p.Results.RespectBoundaries;

    minEpochs = 30; % eBridge asserts on this internally; check first for a clearer message

    rejectLabels = {};
    info = struct('available', false, 'status', 'ok', ...
                  'filtered', logical(doFilter), 'filtBand', [filtLow, filtHigh], ...
                  'epochSec', epochSec, 'numEpochs', 0, ...
                  'EDcutoff', NaN, 'EDscale', NaN, ...
                  'clusters', {{}}, 'pairs', {{}}, ...
                  'nBridgedChannels', 0, 'nBridgedPairs', 0, ...
                  'oversizedCluster', false, 'excessiveBridging', false, ...
                  'rejectLabels', {{}});

    if isempty(which('eBridge'))
        warning('RELAX_bridged_channels:eBridgeNotFound', ...
            ['eBridge.m not found on the MATLAB path. Skipping bridged channel ' ...
             'detection. Install the EEGLAB ebridge plugin to use this screen.']);
        info.status = 'eBridge not found';
        return;
    end

    if doFilter
        workingEEG = RELAX_filtbutter(EEG, filtLow, filtHigh, filtOrder, 'bandpass', 'acausal', respectBoundaries);
    else
        workingEEG = EEG;
    end

    segments = RELAX_get_segments(workingEEG);
    windows = local_epoch_windows(segments, workingEEG.srate, epochSec);
    numEpochs = size(windows, 1);

    if numEpochs < minEpochs
        warning('RELAX_bridged_channels:tooFewEpochs', ...
            ['Only %d epoch(s) of length %.2f s available (eBridge requires at ' ...
             'least %d). Skipping bridged channel detection.'], numEpochs, epochSec, minEpochs);
        info.status = sprintf('too few epochs (%d < %d)', numEpochs, minEpochs);
        return;
    end

    if numEpochs > maxEpochs
        % Deterministic even subsampling (not random) to keep runs reproducible:
        keep = unique(round(linspace(1, numEpochs, maxEpochs)));
        windows = windows(keep, :);
        numEpochs = size(windows, 1);
    end
    info.numEpochs = numEpochs;

    numPts = windows(1,2) - windows(1,1) + 1;
    epochData = zeros(workingEEG.nbchan, numPts, numEpochs);
    for k = 1:numEpochs
        epochData(:,:,k) = double(workingEEG.data(:, windows(k,1):windows(k,2)));
    end

    % eBridge only reads .data, .chanlocs and .srate, so a minimal struct suffices:
    EPO = struct('data', epochData, 'chanlocs', EEG.chanlocs, 'srate', EEG.srate);
    try
        EB = eBridge(EPO, 'BCT', bct, 'BinSize', binSize, 'PlotMode', 0, 'Verbose', 1);
    catch ebException
        warning('RELAX_bridged_channels:eBridgeFailed', ...
            'eBridge failed (%s). Skipping bridged channel detection.', ebException.message);
        info.status = sprintf('eBridge error: %s', ebException.message);
        return;
    end
    clear epochData EPO

    info.available = true;
    if isfield(EB, 'Info')
        if isfield(EB.Info, 'EDcutoff'), info.EDcutoff = EB.Info.EDcutoff; end
        if isfield(EB.Info, 'EDscale'),  info.EDscale  = EB.Info.EDscale;  end
    end

    % EDcutoff == -1 means no near-zero peak was found in the ED distribution, i.e. no
    % bridging - distinct from an error:
    if info.EDcutoff == -1
        info.status = 'no ED cutoff (no bridging)';
        return;
    end
    if ~isfield(EB, 'Bridged') || EB.Bridged.Count == 0
        return;
    end

    labels = {EEG.chanlocs.labels};
    [clusters, pairLabels] = RELAX_group_bridged_channels(labels, EB.Bridged.Pairs);
    if isempty(clusters)
        return;
    end

    clusterSizes = cellfun(@numel, clusters);
    info.clusters = clusters;
    info.pairs = pairLabels;
    info.nBridgedPairs = numel(pairLabels);
    info.nBridgedChannels = sum(clusterSizes);

    memberIdx = [];
    for c = 1:numel(clusters)
        memberIdx = [memberIdx, find(ismember(labels, clusters{c}))]; %#ok<AGROW>
    end
    memberIdx = unique(memberIdx);
    rejectLabels = labels(memberIdx)';
    info.rejectLabels = rejectLabels;

    if any(clusterSizes > maxClusterSz)
        info.oversizedCluster = true;
        warning('RELAX_bridged_channels:oversizedCluster', ...
            ['%d bridged cluster(s) exceed MaxClusterSize (%d); largest is %d ' ...
             'channels. Interpolation over a patch this large is unreliable - ' ...
             'consider excluding this recording.'], ...
            sum(clusterSizes > maxClusterSz), maxClusterSz, max(clusterSizes));
    end
    if info.nBridgedChannels > maxPropChans * numel(labels)
        info.excessiveBridging = true;
        warning('RELAX_bridged_channels:excessiveBridging', ...
            ['%d of %d channels (%.1f%%) are bridged, above the %.0f%% threshold. ' ...
             'Consider excluding this recording.'], ...
            info.nBridgedChannels, numel(labels), 100*info.nBridgedChannels/numel(labels), 100*maxPropChans);
    end
end

function windows = local_epoch_windows(segments, fs, epochSec)
    % Non-overlapping windows of length epochSec within each boundary-respecting
    % segment, never crossing a segment edge.
    windowSamples = round(epochSec * fs);
    windows = zeros(0, 2);
    for s = 1:size(segments,1)
        segStart = segments(s,1);
        segEnd = segments(s,2);
        segLen = segEnd - segStart + 1;
        if segLen < windowSamples
            continue; % segment too short to fit even one window
        end
        winStarts = (segStart:windowSamples:(segEnd - windowSamples + 1))';
        if isempty(winStarts)
            continue;
        end
        winEnds = winStarts + windowSamples - 1;
        windows = [windows; winStarts, winEnds]; %#ok<AGROW>
    end
end
```

- [ ] **Step 5: Run test to verify it passes**

Run: `matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_bridged_channels"`
Expected: PASS

- [ ] **Step 6: Run the full test suite to check for regressions**

Run: `matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"`
Expected: all test files pass, including the two new ones

- [ ] **Step 7: Commit**

```bash
git add RELAX_bridged_channels.m tests/stubs/eBridge.m tests/test_RELAX_bridged_channels.m
git commit -m "Add bridged-channel detection via eBridge, boundary- and RELAX-style"
```

---

## Task 3: Wire bridged-channel rejection into `RELAX_Wrapper_slurm.m`

**Files:**
- Modify: `RELAX_Wrapper_slurm.m:386-388` (insert new section between the impedance block and the PREP comment block)
- Modify: `RELAX_Wrapper_slurm.m:914-937` (fold bridge counts into `EEG.RELAX_issues_to_check`)

**Interfaces:**
- Consumes: `RELAX_bridged_channels(EEG, varargin)` from Task 2 (exact signature above).
- Produces: `EEG.RELAXProcessingExtremeRejections.BridgedChannelToReject` (column cell array, always present, `{}` when off/nothing rejected — same contract as the existing `ImpedanceBasedChannelToReject`/`PREPBasedChannelToReject` fields other wrapper code already reads); `EEG.RELAX.BridgeScreen` (the `info` struct, only set when the stage runs); new `RELAX_cfg` fields `Perform_bridged_channel_rejection` (default 0), `BridgeMaxProportionChannels` (default 0.10), `BridgeMaxClusterSize` (default 4).

- [ ] **Step 1: Insert the new section after the impedance block**

In `RELAX_Wrapper_slurm.m`, immediately after the impedance block's closing (the line reading `end` that closes the `if RELAX_cfg.Perform_impedance_bad_channel_rejection == 1 ... else ... end` at line 386) and before the `%% Clean flat channels...` comment at line 388, insert:

```matlab
    %% Reject bridged electrode clusters detected via eBridge:
    % Bridged electrodes (electrolyte spreading between adjacent contacts) report a
    % near-identical, spatially smeared voltage; RELAX_bridged_channels flags them by an
    % electrical-distance criterion computed directly from the data (Alschuler et al.
    % 2014). Unlike the impedance screen above, this is a data-driven measurement, not
    % metadata that could itself be wrong, so nothing here is capped the way the
    % impedance screen caps its rejections - every flagged channel is removed. A bridge
    % cannot be undone (both electrodes measure one shorted voltage), so the removed
    % channels are reconstructed later by spherical interpolation from the surviving
    % channels, the same treatment any other bad channel receives.
    %
    % Runs before PREP for the same reason the impedance screen does: findNoisyChannels
    % builds a robust average reference and channel-to-channel correlations across the
    % whole montage, so a near-perfectly-correlated bridged pair would skew its verdict
    % on every other channel too.
    %
    % Absent = off: eBridge is an external EEGLAB plugin not bundled with this repo, and
    % this is new behaviour, so it must not switch itself on for a user who has not
    % installed it or opted in.
    if ~isfield(RELAX_cfg,'Perform_bridged_channel_rejection')
        RELAX_cfg.Perform_bridged_channel_rejection = 0;
    end
    if ~isfield(RELAX_cfg,'BridgeMaxProportionChannels')
        RELAX_cfg.BridgeMaxProportionChannels = 0.10; % QC flag only, see RELAX_bridged_channels
    end
    if ~isfield(RELAX_cfg,'BridgeMaxClusterSize')
        RELAX_cfg.BridgeMaxClusterSize = 4; % QC flag only, see RELAX_bridged_channels
    end
    EEG.RELAXProcessingExtremeRejections.BridgedChannelToReject={};
    if RELAX_cfg.Perform_bridged_channel_rejection == 1
        [BridgedChannelToReject, BridgeScreenInfo] = RELAX_bridged_channels(EEG, ...
            'RespectBoundaries', RELAX_cfg.RespectBoundaries, ...
            'MaxProportionChannels', RELAX_cfg.BridgeMaxProportionChannels, ...
            'MaxClusterSize', RELAX_cfg.BridgeMaxClusterSize);
        EEG.RELAXProcessingExtremeRejections.BridgedChannelToReject = BridgedChannelToReject;
        EEG.RELAX.BridgeScreen = BridgeScreenInfo;
        if ~BridgeScreenInfo.available
            warning('Bridge screen unavailable for %s (%s) - no channels rejected on bridging.', FileName, BridgeScreenInfo.status);
        elseif isempty(BridgedChannelToReject)
            fprintf('Bridge screen: no bridged channels detected (%d epochs).\n', BridgeScreenInfo.numEpochs);
        else
            fprintf('Bridge screen: removing %d bridged channel(s) in %d cluster(s): %s\n', ...
                numel(BridgedChannelToReject), numel(BridgeScreenInfo.clusters), strjoin(BridgedChannelToReject', ', '));
            EEG=pop_select(EEG,'nochannel',BridgedChannelToReject);
            EEG = eeg_checkset( EEG );
        end
    else
        fprintf('Skipping bridged-channel rejection (RELAX_cfg.Perform_bridged_channel_rejection=0).\n');
    end

```

- [ ] **Step 2: Fold bridge counts into `EEG.RELAX_issues_to_check`**

Find this existing block (around line 920-921 of `RELAX_Wrapper_slurm.m`):

```matlab
    EEG.RELAX_issues_to_check.Impedance_screen_unavailable=0;
    EEG.RELAX_issues_to_check.Impedance_rejected_too_many_electrodes=0;
    if isfield(EEG,'RELAX') && isfield(EEG.RELAX,'ImpedanceScreen')
        if ~EEG.RELAX.ImpedanceScreen.available
            EEG.RELAX_issues_to_check.Impedance_screen_unavailable=1;
        elseif EEG.RELAX.ImpedanceScreen.capApplied
            EEG.RELAX_issues_to_check.Impedance_rejected_too_many_electrodes=EEG.RELAX.ImpedanceScreen.nOverThreshold;
        end
    end
```

Insert immediately after it:

```matlab
    % The bridge screen enforces no cap on itself (see RELAX_bridged_channels), so
    % there is nothing analogous to "rejected too many" here - just the two QC flags it
    % raises about the recording's own reconstructability:
    EEG.RELAX_issues_to_check.Bridge_screen_unavailable=0;
    EEG.RELAX_issues_to_check.Bridge_oversized_cluster=0;
    EEG.RELAX_issues_to_check.Bridge_excessive_bridging=0;
    if isfield(EEG,'RELAX') && isfield(EEG.RELAX,'BridgeScreen')
        if ~EEG.RELAX.BridgeScreen.available
            EEG.RELAX_issues_to_check.Bridge_screen_unavailable=1;
        else
            if EEG.RELAX.BridgeScreen.oversizedCluster
                EEG.RELAX_issues_to_check.Bridge_oversized_cluster=EEG.RELAX.BridgeScreen.nBridgedChannels;
            end
            if EEG.RELAX.BridgeScreen.excessiveBridging
                EEG.RELAX_issues_to_check.Bridge_excessive_bridging=EEG.RELAX.BridgeScreen.nBridgedChannels;
            end
        end
    end
```

Then find the existing electrode-budget sum a few lines below:

```matlab
    if (EEG.RELAXProcessingExtremeRejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete...
            +EEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete...
            +size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1))...
            >=RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted*size(EEG.allchan,2)
        EEG.RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold=...
            (EEG.RELAXProcessingExtremeRejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete...
            +EEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete...
            +size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1));
    else
        EEG.RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold=0;
    end
```

Replace it with (adding the bridged count into the same budget sum, the same way PREP's count already is):

```matlab
    if (EEG.RELAXProcessingExtremeRejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete...
            +EEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete...
            +size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1)...
            +size(EEG.RELAXProcessingExtremeRejections.BridgedChannelToReject,1))...
            >=RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted*size(EEG.allchan,2)
        EEG.RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold=...
            (EEG.RELAXProcessingExtremeRejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete...
            +EEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete...
            +size(EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject,1)...
            +size(EEG.RELAXProcessingExtremeRejections.BridgedChannelToReject,1));
    else
        EEG.RELAX_issues_to_check.ElectrodeRejectionRecommendationsMetOrExceededThreshold=0;
    end
```

- [ ] **Step 3: Static-analyze the edited file**

Run the MATLAB Code Analyzer on the file (via the `check_matlab_code` MCP tool, or `checkcode('RELAX_Wrapper_slurm.m')` in MATLAB) and confirm no new warnings were introduced by this change (pre-existing warnings elsewhere in the file are out of scope).

- [ ] **Step 4: Manual trace check**

With `RELAX_cfg.Perform_bridged_channel_rejection` left unset (or 0), confirm by reading the diff that:
- `EEG.RELAXProcessingExtremeRejections.BridgedChannelToReject` is still always initialized to `{}` (so downstream code reading `size(...,1)` never errors on a missing field), and
- no other line in the file references `RELAX_cfg.Perform_bridged_channel_rejection`, `BridgeMaxProportionChannels`, or `BridgeMaxClusterSize` before they are defaulted (i.e. the new block's `isfield` guards are the first read of each).

- [ ] **Step 5: Commit**

```bash
git add RELAX_Wrapper_slurm.m
git commit -m "Wire bridged-channel rejection into RELAX_Wrapper_slurm before PREP"
```

---

## Task 4: Config defaults and documentation

**Files:**
- Modify: `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` (~line 338, right after the impedance config block)
- Modify: `CLAUDE.md` (Pipeline architecture stage-4 bullet; External dependencies section)

- [ ] **Step 1: Add default `RELAX_cfg` fields to `RELAX_SET_PARAMETERS_AND_RUN_slurm.m`**

Immediately after this existing line (~338):

```matlab
RELAX_cfg.MaxImpedanceKOhm=500; % kOhm. Strictly greater-than; an unknown (NaN) impedance is never rejected
```

Insert:

```matlab

% Reject clusters of electrically bridged electrodes (electrolyte spreading between
% adjacent contacts), detected via the EEGLAB ebridge plugin's eBridge() from the
% electrical-distance criterion (Alschuler et al. 2014). Off by default: eBridge is an
% external plugin not bundled with this repo (install it to use this screen), and this
% is new behaviour. Every channel in a flagged cluster is removed (a bridge cannot be
% undone) and reconstructed later by spherical interpolation, same as any other bad
% channel. MaxProportionChannels/MaxClusterSize are QC-only warnings on the recording's
% own reconstructability, not caps on what gets rejected - see RELAX_bridged_channels.
RELAX_cfg.Perform_bridged_channel_rejection=0; % 1 = reject bridged channel clusters, 0 = skip
RELAX_cfg.BridgeMaxProportionChannels=0.10; % QC flag only
RELAX_cfg.BridgeMaxClusterSize=4; % QC flag only
```

- [ ] **Step 2: Update `CLAUDE.md`'s pipeline stage-4 bullet**

Find:

```
4. **Bad-channel rejection** — PREP's `findNoisyChannels` (off by default in the SLURM path via
   `RELAX_cfg.Perform_PREP_bad_channel_rejection=0`; see "Boundary/discontinuity handling" below).
```

Replace with:

```
4. **Bad-channel rejection** — three screens run in order, all optional and off/on per
   `RELAX_cfg` flags: recorded-impedance rejection (`RELAX_cfg.Perform_impedance_bad_channel_rejection`,
   `RELAX_impedance_bad_channels.m`), electrically-bridged-cluster rejection
   (`RELAX_cfg.Perform_bridged_channel_rejection`, `RELAX_bridged_channels.m`, needs the
   EEGLAB `ebridge` plugin), then PREP's `findNoisyChannels` (off by default in the SLURM
   path via `RELAX_cfg.Perform_PREP_bad_channel_rejection=0`; see "Boundary/discontinuity
   handling" below). Impedance and bridge screens run first so PREP's whole-montage
   correlation/reference statistics are never skewed by a channel already known to be bad.
```

- [ ] **Step 3: Update `CLAUDE.md`'s External dependencies list**

Find:

```
- ERPLAB, if `RELAX_cfg.NotchFilterType='PMnotch'`
```

Add a line after it:

```
- The EEGLAB `ebridge` plugin, if `RELAX_cfg.Perform_bridged_channel_rejection=1`
```

- [ ] **Step 4: Commit**

```bash
git add RELAX_SET_PARAMETERS_AND_RUN_slurm.m CLAUDE.md
git commit -m "Document bridged-channel rejection defaults and dependency"
```

---

## Self-Review Notes

- **Spec coverage:** placement (after impedance, before PREP) — Task 3 Step 1; ported detection logic — Task 2; ported removal-of-every-cluster-member logic — Task 2 (`RELAX_bridged_channels`, no all-but-one behaviour); style/doc/QC matched to `RELAX_impedance_bad_channels.m` — headers, `info` struct, `EEG.RELAX_issues_to_check` folding all present in Tasks 2-3; config flag defaulting to off — Task 3 Step 1 and Task 4 Step 1.
- **Deliberate deviations from the source** (called out inline in each task, not left implicit): grouping logic split into its own file for testability; no enforced rejection cap (QC-only flags instead), because bridging is data-derived rather than metadata; epoching goes through `RELAX_get_segments` instead of a duplicated `compute_segment_windows`; pre-filtering goes through the existing `RELAX_filtbutter` instead of a duplicated `apply_filter_by_segments`.
- **Out of scope:** `RELAX_Wrapper.m` / `RELAX_Wrapper_bart.m` (the user's request named `RELAX_Wrapper_slurm.m` specifically); vendoring `eBridge.m` itself (it stays an external plugin dependency like PREP/MWF/FieldTrip).
