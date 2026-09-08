# Boundary-Aware Preprocessing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the `RELAX_Wrapper_slurm.m` preprocessing path handle `boundary` events correctly, by processing each stretch of data between boundaries independently with reflection-padded edges.

**Architecture:** Three new shared helper functions turn "a stretch of data between boundary events" into a first-class concept: one derives the segment list from `EEG.event`, one reflection-pads a data block, and one applies an arbitrary function to each segment independently. The functions that currently smear across discontinuities (`RELAX_filtbutter`, `RELAX_blinks_IQR_method`, `RELAX_targeted_wICA`) route their time-domain operations through those helpers, gated behind a new `RELAX_cfg.RespectBoundaries` flag that defaults off so sibling wrappers are unaffected. Two further stages that are boundary-blind but not needed (the dead SER/ARR marking block, and PREP's `findNoisyChannels`) are switched off by flag rather than fixed, and a startup guard errors if any known-unsafe stage is enabled on data that contains boundaries.

**Tech Stack:** MATLAB R2023a (`module load matlab`), EEGLAB 2026.0.0 at `/home/nga4004/tools/eeglab2026.0.0`, MATLAB Signal Processing Toolbox (`butter`, `filtfilt`), MATLAB Wavelet Toolbox (`swt`, `iswt`, `wthresh`, `ddencmp`).

**Spec:** `docs/superpowers/specs/2026-09-08-boundary-aware-preprocessing-design.md`

## Global Constraints

- **Repo root:** `/home/nga4004/gitRepos/RELAX-GEDAI`. Branch: `tmseeg`. All paths below are relative to the repo root unless absolute.
- **No test framework exists in this repo.** Tests are plain MATLAB scripts under `tests/`, named `test_*.m`, that use `assert` and are discovered by `tests/run_all_tests.m`.
- **Unit test command** (Tasks 1-6; fast, safe to run interactively on a login node):
  ```bash
  module load matlab
  cd /home/nga4004/gitRepos/RELAX-GEDAI
  matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
  ```
  `-batch` exits non-zero if any assertion fails, so the exit code is the pass/fail signal.
- **Never run computationally intensive work on a login node.** Login nodes are named `scu-login*`. Task 7's integration run must be submitted via `sbatch`/`srun` to the `sackler` or `scu-cpu` partition. Tasks 1-6 unit tests take seconds and are fine interactively.
- **Do not add EEGLAB as a dependency of the Task 1-6 unit tests.** They build EEG-like structs by hand (`struct('data',...,'srate',...,'pnts',...,'event',...,'trials',1)`). Only Task 7 needs EEGLAB loaded.
- **Backward compatibility is a hard requirement.** Every modified function must produce bitwise-identical output to its pre-change version when `RELAX_cfg.RespectBoundaries` is `0` or the flag is absent. `RELAX_Wrapper.m` and `RELAX_Wrapper_bart.m` are never modified.
- **New config field names, exact spelling:** `RELAX_cfg.RespectBoundaries`, `RELAX_cfg.Perform_PREP_bad_channel_rejection`, `RELAX_cfg.ComputeSER_ARR`.
- **Reflection convention:** odd (antisymmetric) reflection, `2*x(edge) - flip(interior)`, matching `filtfilt`'s internal convention. Never even/mirror reflection.
- **Commit after every task.** Do not squash tasks together.

---

### Task 1: Segment helpers and test infrastructure

**Files:**
- Create: `RELAX_get_segments.m`
- Create: `RELAX_reflect_pad.m`
- Create: `RELAX_apply_by_segments.m`
- Create: `tests/run_all_tests.m`
- Test: `tests/test_RELAX_get_segments.m`
- Test: `tests/test_RELAX_reflect_pad.m`
- Test: `tests/test_RELAX_apply_by_segments.m`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `segments = RELAX_get_segments(EEG)` — `EEG` needs `.pnts` and optionally `.event`. Returns `[nSegments x 2]` double, columns `[startSample endSample]`, 1-based inclusive.
  - `padded = RELAX_reflect_pad(x, p)` — `x` is `[nChannels x nSamples]`, `p` scalar. Returns `[nChannels x (nSamples + 2*pEffective)]` where `pEffective = min(p, nSamples-1)`.
  - `out = RELAX_apply_by_segments(data, segments, fn, p)` — `data` is `[nChannels x nSamples]`, `fn` a function handle `[nChan x N] -> [nChan x N]`. Returns `[nChannels x nSamples]`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_RELAX_get_segments.m`:

```matlab
function test_RELAX_get_segments()
% Segment derivation from boundary events.

    % No event field at all -> one whole-record segment
    EEG = struct('pnts', 200);
    assert(isequal(RELAX_get_segments(EEG), [1 200]), 'no event field should give one segment');

    % Empty event struct -> one whole-record segment
    EEG = struct('pnts', 200, 'event', struct('type', {}, 'latency', {}));
    assert(isequal(RELAX_get_segments(EEG), [1 200]), 'empty events should give one segment');

    % Events present but none are boundaries -> one whole-record segment
    EEG = struct('pnts', 200, 'event', struct('type', {'X','Y'}, 'latency', {50, 120}));
    assert(isequal(RELAX_get_segments(EEG), [1 200]), 'non-boundary events should give one segment');

    % Half-integer boundary (EEGLAB convention) splits without duplicating a sample
    EEG = struct('pnts', 200, 'event', struct('type', {'boundary'}, 'latency', {100.5}));
    assert(isequal(RELAX_get_segments(EEG), [1 100; 101 200]), 'half-integer boundary split');

    % Integer boundary splits without duplicating a sample
    EEG = struct('pnts', 200, 'event', struct('type', {'boundary'}, 'latency', {100}));
    assert(isequal(RELAX_get_segments(EEG), [1 100; 101 200]), 'integer boundary split');

    % Boundaries out of order in the event array are sorted
    EEG = struct('pnts', 300, 'event', struct('type', {'boundary','boundary'}, 'latency', {200.5, 100.5}));
    assert(isequal(RELAX_get_segments(EEG), [1 100; 101 200; 201 300]), 'boundaries should be sorted');

    % Numeric event types must not crash the type test
    EEG = struct('pnts', 200, 'event', struct('type', {1, 'boundary'}, 'latency', {10, 100.5}));
    assert(isequal(RELAX_get_segments(EEG), [1 100; 101 200]), 'numeric event types tolerated');

    % Segments always tile the record with no gaps and no overlaps
    EEG = struct('pnts', 1000, 'event', struct('type', {'boundary','boundary','boundary'}, ...
                                               'latency', {250.5, 500.5, 750.5}));
    segs = RELAX_get_segments(EEG);
    assert(segs(1,1) == 1, 'first segment starts at 1');
    assert(segs(end,2) == 1000, 'last segment ends at pnts');
    assert(all(segs(2:end,1) == segs(1:end-1,2) + 1), 'segments must tile without gap or overlap');

    fprintf('test_RELAX_get_segments: PASS\n');
end
```

Create `tests/test_RELAX_reflect_pad.m`:

```matlab
function test_RELAX_reflect_pad()
% Odd reflection padding.

    % A linear ramp must be continued exactly by odd reflection
    x = [1 2 3 4 5];
    padded = RELAX_reflect_pad(x, 2);
    assert(isequal(padded, [-1 0 1 2 3 4 5 6 7]), 'odd reflection must continue a linear ramp');

    % The original data must be preserved in the middle
    assert(isequal(padded(3:7), x), 'interior must be untouched');

    % p is clamped to nSamples-1
    x = [1 2 3];
    padded = RELAX_reflect_pad(x, 10);
    assert(size(padded,2) == 3 + 2*2, 'p must clamp to nSamples-1');

    % p of zero returns the input unchanged
    x = [1 2 3];
    assert(isequal(RELAX_reflect_pad(x, 0), x), 'p=0 returns input unchanged');

    % Multi-channel shape is preserved
    x = [1 2 3 4 5; 10 20 30 40 50];
    padded = RELAX_reflect_pad(x, 2);
    assert(isequal(size(padded), [2 9]), 'multi-channel output shape');
    assert(isequal(padded(2,:), [8 0 10 20 30 40 50 60 70]), 'second channel reflected independently');

    fprintf('test_RELAX_reflect_pad: PASS\n');
end
```

Create `tests/test_RELAX_apply_by_segments.m`:

```matlab
function test_RELAX_apply_by_segments()
% Per-segment application with reflection padding.

    data = reshape(1:20, 1, 20);
    segments = [1 10; 11 20];

    % Identity function must round-trip the data exactly
    out = RELAX_apply_by_segments(data, segments, @(d) d, 3);
    assert(isequal(out, data), 'identity fn must round-trip exactly');

    % A pointwise function must apply within every segment
    out = RELAX_apply_by_segments(data, segments, @(d) d*2, 3);
    assert(isequal(out, data*2), 'pointwise fn must apply across all segments');

    % Output size always matches input size
    assert(isequal(size(out), size(data)), 'output size must match input size');

    % ISOLATION: scaling segment 2 by 1e6 must not change segment 1's output at all.
    % A function that depends on the whole block it is given (mean subtraction) would
    % leak across segments if the implementation concatenated them.
    dataA = [ones(1,10), ones(1,10)];
    dataB = [ones(1,10), ones(1,10)*1e6];
    fn = @(d) d - mean(d,2);
    outA = RELAX_apply_by_segments(dataA, segments, fn, 3);
    outB = RELAX_apply_by_segments(dataB, segments, fn, 3);
    assert(isequal(outA(1:10), outB(1:10)), ...
        'segment 1 output must be bitwise unchanged when segment 2 changes');

    % Multi-channel data
    data2 = [reshape(1:20,1,20); reshape(21:40,1,20)];
    out2 = RELAX_apply_by_segments(data2, segments, @(d) d*3, 2);
    assert(isequal(out2, data2*3), 'multi-channel pointwise fn');

    % A single whole-record segment behaves like a plain call
    out3 = RELAX_apply_by_segments(data, [1 20], @(d) d*2, 3);
    assert(isequal(out3, data*2), 'single segment case');

    fprintf('test_RELAX_apply_by_segments: PASS\n');
end
```

Create `tests/run_all_tests.m`:

```matlab
function run_all_tests()
% Discover and run every tests/test_*.m, reporting PASS/FAIL per test.
% Errors at the end if any test failed, so `matlab -batch` exits non-zero.

    thisDir = fileparts(mfilename('fullpath'));
    testFiles = dir(fullfile(thisDir, 'test_*.m'));

    if isempty(testFiles)
        error('run_all_tests: no test_*.m files found in %s', thisDir);
    end

    failures = {};
    for k = 1:numel(testFiles)
        [~, testName] = fileparts(testFiles(k).name);
        try
            feval(testName);
        catch err
            failures{end+1} = sprintf('%s: %s', testName, err.message); %#ok<AGROW>
            fprintf('%s: FAIL - %s\n', testName, err.message);
        end
    end

    fprintf('\n%d/%d test files passed.\n', numel(testFiles) - numel(failures), numel(testFiles));

    if ~isempty(failures)
        error('run_all_tests: %d test file(s) failed:\n%s', numel(failures), strjoin(failures, '\n'));
    end
end
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: all three test files FAIL with "Unrecognized function or variable 'RELAX_get_segments'" (and the equivalent for the other two), and the command exits non-zero.

- [ ] **Step 3: Write the three helper implementations**

Create `RELAX_get_segments.m`:

```matlab
% RELAX_get_segments() - Return the sample ranges of contiguous data between
%                        boundary events.
%
% Usage:
%   >> segments = RELAX_get_segments(EEG);
%
% Inputs:
%   EEG - EEGLAB EEG structure. Requires .pnts; uses .event if present.
%
% Outputs:
%   segments - [nSegments x 2] matrix of [startSample endSample], 1-based and
%              inclusive, one row per stretch of contiguous data. Data with no
%              boundary events returns a single row [1 EEG.pnts], so callers
%              degrade to whole-record behaviour on ordinary continuous data.
%
% Boundary latencies are converted with the integer-safe form
%   starts = [1, floor(b)+1],  ends = [floor(b), pnts]
% which handles both EEGLAB's convention of placing boundaries at x.5 (between
% samples) and tools that write integer boundary latencies, without ever
% assigning the same sample to two segments.

function segments = RELAX_get_segments(EEG)

    nPoints = EEG.pnts;

    boundaries = [];
    if isfield(EEG,'event') && ~isempty(EEG.event)
        isBoundary = false(1,numel(EEG.event));
        for e = 1:numel(EEG.event)
            thisType = EEG.event(e).type;
            if (ischar(thisType) || isstring(thisType)) && strcmp(thisType,'boundary')
                isBoundary(e) = true;
            end
        end
        boundaries = sort([EEG.event(isBoundary).latency]);
    end

    % Ignore boundaries outside the data, and boundaries at the very first or last
    % sample, which do not split the record into two usable pieces:
    boundaries = boundaries(boundaries > 1 & boundaries < nPoints);

    starts = [1, floor(boundaries) + 1];
    ends   = [floor(boundaries), nPoints];

    valid = starts <= ends;
    segments = [starts(valid)', ends(valid)'];
end
```

Create `RELAX_reflect_pad.m`:

```matlab
% RELAX_reflect_pad() - Odd (antisymmetric) reflection padding about the edge samples.
%
% Usage:
%   >> padded = RELAX_reflect_pad(x, p);
%
% Inputs:
%   x - [nChannels x nSamples] data
%   p - number of samples to add at each end. Clamped to nSamples-1.
%
% Outputs:
%   padded - [nChannels x (nSamples + 2*p)] data
%
% Uses odd reflection, 2*x(edge) - flip(interior), which is the convention
% filtfilt applies internally. Even (mirror) reflection would place a cusp at the
% segment edge that a high pass filter then rings on; odd reflection continues the
% local trend across the edge instead, which matters here because the segments
% between boundary events are only ~8 s long.

function padded = RELAX_reflect_pad(x, p)

    nSamples = size(x,2);
    p = min(p, nSamples-1);

    if p < 1
        padded = x;
        return
    end

    leftPad  = 2*x(:,1)   - x(:, p+1:-1:2);
    rightPad = 2*x(:,end) - x(:, end-1:-1:end-p);
    padded = [leftPad, x, rightPad];
end
```

Create `RELAX_apply_by_segments.m`:

```matlab
% RELAX_apply_by_segments() - Apply a function to each segment of data independently,
%                             with reflection padding at each segment edge.
%
% Usage:
%   >> out = RELAX_apply_by_segments(data, segments, fn, p);
%
% Inputs:
%   data     - [nChannels x nSamples] data
%   segments - [nSegments x 2] matrix of [startSample endSample], as returned by
%              RELAX_get_segments
%   fn       - function handle taking [nChannels x N] and returning [nChannels x N]
%   p        - reflection pad width in samples, clamped per segment to segLen-1
%
% Outputs:
%   out - [nChannels x nSamples], the same size as data. Any samples not covered by
%         a segment are passed through unchanged.
%
% Each segment is padded, transformed and trimmed on its own, so no segment can
% influence any other. Callers must derive `segments` from the CURRENT EEG.event:
% eeg_eegrej adds boundaries part-way through the pipeline, so a segment list must
% never be cached across stages.

function out = RELAX_apply_by_segments(data, segments, fn, p)

    out = data;

    for s = 1:size(segments,1)
        segStart = segments(s,1);
        segEnd   = segments(s,2);
        segData  = data(:, segStart:segEnd);

        thisPad  = min(p, size(segData,2)-1);
        if thisPad < 0
            thisPad = 0;
        end

        padded   = RELAX_reflect_pad(segData, thisPad);
        filtered = fn(padded);

        out(:, segStart:segEnd) = filtered(:, thisPad+1 : thisPad+size(segData,2));
    end
end
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: `3/3 test files passed.` and exit code 0.

- [ ] **Step 5: Commit**

```bash
cd /home/nga4004/gitRepos/RELAX-GEDAI
git add RELAX_get_segments.m RELAX_reflect_pad.m RELAX_apply_by_segments.m tests/
git commit -m "Add boundary segment helpers and MATLAB test harness"
```

---

### Task 2: Config flags, startup guard, and disabling the unneeded stages

Turns off the two boundary-blind stages that are not needed — the dead SER/ARR marking block and PREP's `findNoisyChannels` — and adds the guard that errors if a known-unsafe stage is enabled on data containing boundaries.

**Files:**
- Create: `RELAX_check_boundary_support.m`
- Modify: `RELAX_Wrapper_slurm.m:324` (PREP call), `RELAX_Wrapper_slurm.m:367` (SER/ARR block guard)
- Modify: `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` (add the three new config fields)
- Test: `tests/test_RELAX_check_boundary_support.m`

**Interfaces:**
- Consumes: `RELAX_get_segments(EEG)` from Task 1.
- Produces: `RELAX_check_boundary_support(EEG, RELAX_cfg)` — returns nothing, throws `error` with identifier `RELAX:boundaryUnsafeStage` when an unsafe stage is enabled on data with boundaries.

- [ ] **Step 1: Write the failing test**

Create `tests/test_RELAX_check_boundary_support.m`:

```matlab
function test_RELAX_check_boundary_support()
% Startup guard for stages that cannot handle boundary events.

    EEG_with    = struct('pnts', 200, 'event', struct('type', {'boundary'}, 'latency', {100.5}));
    EEG_without = struct('pnts', 200, 'event', struct('type', {'X'}, 'latency', {100}));

    safe_cfg = struct( ...
        'Perform_PREP_bad_channel_rejection', 0, ...
        'Do_MWF_Once', 0, 'Do_MWF_Twice', 0, 'Do_MWF_Thrice', 0, ...
        'Run_GEDAI', 0, ...
        'NotchFilterType', 'ZaplinePlus', ...
        'Perform_wICA_on_ICLabel', 0, ...
        'ComputeSER_ARR', 0, ...
        'Perform_ICA_subtract', 0);

    % A fully safe config must not throw, with or without boundaries
    RELAX_check_boundary_support(EEG_with, safe_cfg);
    RELAX_check_boundary_support(EEG_without, safe_cfg);

    % Each unsafe stage must throw when boundaries are present
    unsafe_fields = {'Perform_PREP_bad_channel_rejection', 'Do_MWF_Once', 'Do_MWF_Twice', ...
                     'Do_MWF_Thrice', 'Run_GEDAI', 'Perform_wICA_on_ICLabel', 'ComputeSER_ARR'};
    for k = 1:numel(unsafe_fields)
        cfg = safe_cfg;
        cfg.(unsafe_fields{k}) = 1;
        threw = false;
        try
            RELAX_check_boundary_support(EEG_with, cfg);
        catch err
            threw = true;
            assert(strcmp(err.identifier, 'RELAX:boundaryUnsafeStage'), ...
                'wrong error identifier for %s', unsafe_fields{k});
        end
        assert(threw, 'guard must throw when %s is enabled with boundaries present', unsafe_fields{k});

        % ...and must NOT throw on data with no boundaries
        RELAX_check_boundary_support(EEG_without, cfg);
    end

    % The Butterworth notch is unsafe; ZaplinePlus and PMnotch are not checked here
    cfg = safe_cfg;
    cfg.NotchFilterType = 'Butterworth';
    threw = false;
    try
        RELAX_check_boundary_support(EEG_with, cfg);
    catch err
        threw = true;
        assert(strcmp(err.identifier, 'RELAX:boundaryUnsafeStage'), 'wrong identifier for notch');
    end
    assert(threw, 'guard must throw for Butterworth notch with boundaries present');

    % Perform_ICA_subtract is spatial-only and must NEVER be guarded
    cfg = safe_cfg;
    cfg.Perform_ICA_subtract = 1;
    RELAX_check_boundary_support(EEG_with, cfg);

    % Missing fields must be tolerated (treated as disabled), not error
    RELAX_check_boundary_support(EEG_with, struct('Run_GEDAI', 0));

    fprintf('test_RELAX_check_boundary_support: PASS\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: `test_RELAX_check_boundary_support` FAILs with "Unrecognized function or variable 'RELAX_check_boundary_support'". The three Task 1 tests still pass.

- [ ] **Step 3: Write the guard implementation**

Create `RELAX_check_boundary_support.m`:

```matlab
% RELAX_check_boundary_support() - Error if a boundary-blind pipeline stage is enabled
%                                  on data that contains boundary events.
%
% Usage:
%   >> RELAX_check_boundary_support(EEG, RELAX_cfg);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%   RELAX_cfg - RELAX configuration structure
%
% Outputs:
%   none. Throws error 'RELAX:boundaryUnsafeStage' if an unsafe stage is enabled.
%
% The stages listed here all process the record as one unbroken time series and
% would smear across discontinuities. See
% docs/superpowers/specs/2026-09-08-boundary-aware-preprocessing-design.md for the
% analysis of each. RELAX_cfg.Perform_ICA_subtract is deliberately NOT listed: it
% is ICA plus pop_subcomp, a purely spatial projection with no time-domain
% filtering and no wavelets, so it is unaffected by discontinuities.

function RELAX_check_boundary_support(EEG, RELAX_cfg)

    segments = RELAX_get_segments(EEG);
    if size(segments,1) <= 1
        return % no boundaries, nothing to guard against
    end

    unsafe = {};

    flagChecks = { ...
        'Perform_PREP_bad_channel_rejection', 'PREP findNoisyChannels tiles its correlation windows flatly from sample 1, so windows straddle boundaries'; ...
        'Do_MWF_Once',                        'MWF delay-stacking builds one lagged covariance across the whole record, so lags span boundaries'; ...
        'Do_MWF_Twice',                       'MWF delay-stacking builds one lagged covariance across the whole record, so lags span boundaries'; ...
        'Do_MWF_Thrice',                      'MWF delay-stacking builds one lagged covariance across the whole record, so lags span boundaries'; ...
        'Run_GEDAI',                          'GEDAI has not been checked for boundary handling'; ...
        'Perform_wICA_on_ICLabel',            'RELAX_wICA_on_ICLabel_artifacts applies wavelet transforms across the whole record'; ...
        'ComputeSER_ARR',                     'the SER/ARR marking block calls RELAX_muscle, RELAX_horizontaleye, RELAX_drift and RELAX_pad_brief_mask_periods, all of which are boundary-blind'};

    for k = 1:size(flagChecks,1)
        thisField = flagChecks{k,1};
        if isfield(RELAX_cfg, thisField) && RELAX_cfg.(thisField) == 1
            unsafe{end+1} = sprintf('  RELAX_cfg.%s: %s', thisField, flagChecks{k,2}); %#ok<AGROW>
        end
    end

    if isfield(RELAX_cfg,'NotchFilterType') && strcmp(RELAX_cfg.NotchFilterType,'Butterworth')
        unsafe{end+1} = '  RELAX_cfg.NotchFilterType=''Butterworth'': the Butterworth notch filters across the whole record';
    end

    if ~isempty(unsafe)
        error('RELAX:boundaryUnsafeStage', ...
            ['This data contains %d boundary events (%d segments), but the following enabled ' ...
             'stages do not handle discontinuities:\n%s\n' ...
             'Either disable them, or implement boundary handling for them first. See ' ...
             'docs/superpowers/specs/2026-09-08-boundary-aware-preprocessing-design.md.'], ...
            size(segments,1)-1, size(segments,1), strjoin(unsafe, '\n'));
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: `4/4 test files passed.`

- [ ] **Step 5: Gate the PREP call in the wrapper**

In `RELAX_Wrapper_slurm.m`, replace this block (currently at lines 321-331):

```matlab
    %% Clean flat channels and bad channels showing improbable data:
    % PREP pipeline: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4471356/
    addpath(genpath(PrepFolderLocation{1,1})); %  1.1.4: fix error where PREP seems to be removed from the path after an EEGLAB update
    noisyOut = findNoisyChannels(EEG);  
    EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject={};
    for x=1:size(noisyOut.noisyChannels.all,2) % loop through output of PREP's findNoisyChannels and take a record of noisy electrodes for deletion:
        PREPBasedChannelToReject{x}=EEG.chanlocs(noisyOut.noisyChannels.all(x)).labels;
        EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject = PREPBasedChannelToReject';
    end
```

with:

```matlab
    %% Clean flat channels and bad channels showing improbable data:
    % PREP pipeline: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4471356/
    % NOTE: findNoisyChannels states its own precondition "No segments of the EEG data
    % have been removed" and tiles its correlation windows flatly from sample 1, so on
    % data with boundary events a large fraction of windows straddle a discontinuity and
    % channel rejection becomes close to arbitrary. It is therefore skipped by default
    % here; RELAX_excluding_channels_and_epoching performs its own boundary-safe channel
    % rejection immediately afterwards. Set Perform_PREP_bad_channel_rejection=1 to
    % restore the original behaviour on data with no discontinuities.
    if ~isfield(RELAX_cfg,'Perform_PREP_bad_channel_rejection')
        RELAX_cfg.Perform_PREP_bad_channel_rejection = 1; % absent = original behaviour
    end
    EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject={};
    noisyOut = struct('noisyChannels', struct('all', []));
    if RELAX_cfg.Perform_PREP_bad_channel_rejection == 1
        addpath(genpath(PrepFolderLocation{1,1})); %  1.1.4: fix error where PREP seems to be removed from the path after an EEGLAB update
        noisyOut = findNoisyChannels(EEG);
        for x=1:size(noisyOut.noisyChannels.all,2) % loop through output of PREP's findNoisyChannels and take a record of noisy electrodes for deletion:
            PREPBasedChannelToReject{x}=EEG.chanlocs(noisyOut.noisyChannels.all(x)).labels;
            EEG.RELAXProcessingExtremeRejections.PREPBasedChannelToReject = PREPBasedChannelToReject';
        end
    else
        fprintf('Skipping PREP findNoisyChannels (RELAX_cfg.Perform_PREP_bad_channel_rejection=0).\n');
    end
```

Leave the two `fprintf`/`pop_select` lines that follow (reporting and deleting `noisyOut.noisyChannels.all`) exactly as they are — with `noisyOut.noisyChannels.all` empty they correctly report "No noisy channels detected by PREP" and delete nothing.

- [ ] **Step 6: Gate the dead SER/ARR block in the wrapper**

In `RELAX_Wrapper_slurm.m`, change the `if` on the line that currently reads:

```matlab
    if RELAX_cfg.computecleanedmetrics==1 && (RELAX_cfg.Do_MWF_Once==0 || RELAX_cfg.Do_MWF_Twice==0 || RELAX_cfg.Do_MWF_Thrice==0)
```

to:

```matlab
    % NOTE: this block exists only to build NoiseMaskFullLengthR1, whose sole consumer
    % is RELAX_metrics_final_SER_and_ARR - and that call is commented out further down
    % this file. It is off by default because RELAX_muscle, RELAX_horizontaleye,
    % RELAX_drift and RELAX_pad_brief_mask_periods are all boundary-blind and would need
    % boundary handling before this block can be re-enabled on discontinuous data.
    if ~isfield(RELAX_cfg,'ComputeSER_ARR'); RELAX_cfg.ComputeSER_ARR = 0; end
    if RELAX_cfg.ComputeSER_ARR==1 && RELAX_cfg.computecleanedmetrics==1 && (RELAX_cfg.Do_MWF_Once==0 || RELAX_cfg.Do_MWF_Twice==0 || RELAX_cfg.Do_MWF_Thrice==0)
```

- [ ] **Step 7: Add the config fields to the parameters script**

In `RELAX_SET_PARAMETERS_AND_RUN_slurm.m`, immediately after the line:

```matlab
RELAX_cfg.Run_GEDAI=0;  % IMK added
```

insert:

```matlab
% Boundary/discontinuity handling. The input files for this pipeline contain boundary
% events marking periods excised by TMS artifact removal, with only ~8 s of continuous
% data between them. See
% docs/superpowers/specs/2026-09-08-boundary-aware-preprocessing-design.md
RELAX_cfg.RespectBoundaries=1; % 1 = process each stretch of data between boundary events separately, with reflection-padded edges
RELAX_cfg.Perform_PREP_bad_channel_rejection=0; % 0 = skip PREP findNoisyChannels (boundary-blind); RELAX_excluding_channels_and_epoching still rejects bad channels boundary-safely
RELAX_cfg.ComputeSER_ARR=0; % 0 = skip the SER/ARR artifact marking block (its only consumer, RELAX_metrics_final_SER_and_ARR, is commented out, and the marking functions are boundary-blind)
```

- [ ] **Step 8: Call the guard from the wrapper**

In `RELAX_Wrapper_slurm.m`, immediately after the line:

```matlab
    RELAX_cfg.ms_per_sample=(1000/EEG.srate);
```

insert:

```matlab
    % Fail fast if this file contains discontinuities and a boundary-blind stage is on:
    RELAX_check_boundary_support(EEG, RELAX_cfg);
```

- [ ] **Step 9: Verify the wrapper still parses and the guard is wired**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); \
  checkcode('RELAX_Wrapper_slurm.m','-string'), \
  checkcode('RELAX_SET_PARAMETERS_AND_RUN_slurm.m','-string'), \
  run_all_tests"
```

Expected: `checkcode` reports no syntax errors (style warnings such as `AGROW` are acceptable and pre-existing), and `4/4 test files passed.`

Behavioural verification of the wrapper edits happens in Task 7's end-to-end run — there is no way to exercise them without the full EEGLAB stack.

- [ ] **Step 10: Commit**

```bash
cd /home/nga4004/gitRepos/RELAX-GEDAI
git add RELAX_check_boundary_support.m tests/test_RELAX_check_boundary_support.m \
        RELAX_Wrapper_slurm.m RELAX_SET_PARAMETERS_AND_RUN_slurm.m
git commit -m "Skip PREP and the dead SER/ARR block; add boundary-support startup guard"
```

---

### Task 3: Boundary-aware `RELAX_filtbutter`

**Files:**
- Modify: `RELAX_filtbutter.m` (signature, and the filtering block at lines 105-118)
- Modify: `RELAX_Wrapper_slurm.m` (call sites at lines 231, 233, 238, 240, 257, 263, 280, 285, 317, 648)
- Test: `tests/test_RELAX_filtbutter_boundaries.m`

**Interfaces:**
- Consumes: `RELAX_get_segments`, `RELAX_apply_by_segments` from Task 1.
- Produces: `EEG = RELAX_filtbutter(EEG, high, low, ord, type, causal_or_acausal_filter, respect_boundaries)` — the 7th argument is new, logical, and defaults to `false`. All six-argument call sites keep working unchanged.

- [ ] **Step 1: Write the failing test**

Create `tests/test_RELAX_filtbutter_boundaries.m`:

```matlab
function test_RELAX_filtbutter_boundaries()
% Segment-wise, reflection-padded Butterworth filtering.

    srate = 500;
    nSamples = 8000;
    t = (0:nSamples-1)/srate;

    % Two segments of clean 10 Hz alpha, with a large DC step across the boundary.
    % A boundary-blind highpass smears that step deep into both segments.
    signal = sin(2*pi*10*t);
    signal(4001:end) = signal(4001:end) + 500; % DC step at the boundary

    EEG = struct('data', signal, 'srate', srate, 'pnts', nSamples, 'trials', 1, ...
                 'event', struct('type', {'boundary'}, 'latency', {4000.5}));

    % --- Backward compatibility: 6-arg call is unchanged -------------------------
    EEG_legacy_6arg = RELAX_filtbutter(EEG, 1, 200, 4, 'bandpass', 'acausal');
    EEG_legacy_7arg = RELAX_filtbutter(EEG, 1, 200, 4, 'bandpass', 'acausal', false);
    assert(isequal(EEG_legacy_6arg.data, EEG_legacy_7arg.data), ...
        'respect_boundaries=false must be bitwise identical to the 6-argument call');

    % --- Segment path runs and preserves shape -----------------------------------
    EEG_seg = RELAX_filtbutter(EEG, 1, 200, 4, 'bandpass', 'acausal', true);
    assert(isequal(size(EEG_seg.data), size(signal)), 'output size must be preserved');
    assert(~any(isnan(EEG_seg.data)), 'output must contain no NaNs');

    % --- ISOLATION: segment 1 must not change when segment 2 changes -------------
    signalB = signal;
    signalB(4001:end) = signalB(4001:end) * 100;
    EEG_B = EEG; EEG_B.data = signalB;
    out_A = RELAX_filtbutter(EEG,   1, 200, 4, 'bandpass', 'acausal', true);
    out_B = RELAX_filtbutter(EEG_B, 1, 200, 4, 'bandpass', 'acausal', true);
    assert(isequal(out_A.data(1:4000), out_B.data(1:4000)), ...
        'segment 1 must be bitwise unchanged when segment 2 is scaled');

    % The boundary-blind path does NOT have that property - this is the bug being fixed
    blind_A = RELAX_filtbutter(EEG,   1, 200, 4, 'bandpass', 'acausal', false);
    blind_B = RELAX_filtbutter(EEG_B, 1, 200, 4, 'bandpass', 'acausal', false);
    assert(~isequal(blind_A.data(1:4000), blind_B.data(1:4000)), ...
        'sanity check: the legacy path is expected to leak across the boundary');

    % --- Segment-wise result matches filtering that segment on its own ------------
    seg1 = struct('data', signal(1:4000), 'srate', srate, 'pnts', 4000, 'trials', 1, 'event', []);
    seg1_alone = RELAX_filtbutter(seg1, 1, 200, 4, 'bandpass', 'acausal', true);
    assert(max(abs(out_A.data(1:4000) - seg1_alone.data)) < 1e-9, ...
        'a segment filtered in place must match that segment filtered standalone');

    % --- The 10 Hz signal survives in the interior of each segment ---------------
    % Compare against the known input away from the edges (1 s in from each end).
    interior = 501:3500;
    recovered = out_A.data(interior);
    expected  = signal(interior);
    assert(corr(recovered(:), expected(:)) > 0.99, ...
        'the 10 Hz component must survive the bandpass in the segment interior');

    % --- No boundaries present: segment path stays close to the legacy path -------
    EEG_nb = struct('data', sin(2*pi*10*t), 'srate', srate, 'pnts', nSamples, 'trials', 1, 'event', []);
    nb_legacy = RELAX_filtbutter(EEG_nb, 1, 200, 4, 'bandpass', 'acausal', false);
    nb_seg    = RELAX_filtbutter(EEG_nb, 1, 200, 4, 'bandpass', 'acausal', true);
    mid = 1001:7000; % 2 s in from each end, beyond any edge-treatment difference
    assert(max(abs(nb_legacy.data(mid) - nb_seg.data(mid))) < 1e-6, ...
        'with no boundaries the two paths must agree in the interior');

    % --- Bookkeeping: filtering must not disturb events, pnts or srate -------------
    assert(EEG_seg.pnts == EEG.pnts, 'pnts must be preserved');
    assert(numel(EEG_seg.event) == numel(EEG.event), 'event count must be preserved');
    assert(EEG_seg.event(1).latency == EEG.event(1).latency, 'boundary latency must be preserved');
    assert(strcmp(EEG_seg.event(1).type, 'boundary'), 'boundary event type must be preserved');

    % --- Highpass and lowpass types also work through the segment path ------------
    hp = RELAX_filtbutter(EEG, 1, [], 4, 'highpass', 'acausal', true);
    assert(isequal(size(hp.data), size(signal)) && ~any(isnan(hp.data)), 'highpass segment path');
    lp = RELAX_filtbutter(EEG, [], 45, 4, 'lowpass', 'acausal', true);
    assert(isequal(size(lp.data), size(signal)) && ~any(isnan(lp.data)), 'lowpass segment path');
    bs = RELAX_filtbutter(EEG, 56, 64, 4, 'bandstop', 'acausal', true);
    assert(isequal(size(bs.data), size(signal)) && ~any(isnan(bs.data)), 'bandstop segment path');

    fprintf('test_RELAX_filtbutter_boundaries: PASS\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: `test_RELAX_filtbutter_boundaries` FAILs — `RELAX_filtbutter` currently takes at most 6 arguments, so the 7-argument call errors with "Too many input arguments".

- [ ] **Step 3: Add the `respect_boundaries` argument to `RELAX_filtbutter`**

In `RELAX_filtbutter.m`, change the function signature line:

```matlab
function EEG = RELAX_filtbutter(EEG, high, low, ord, type, causal_or_acausal_filter)
```

to:

```matlab
function EEG = RELAX_filtbutter(EEG, high, low, ord, type, causal_or_acausal_filter, respect_boundaries)
```

Immediately after the existing block:

```matlab
if exist('causal_or_acausal_filter', 'var')==0
        causal_or_acausal_filter='acausal';
end
```

insert:

```matlab
if exist('respect_boundaries', 'var')==0
        respect_boundaries=false;
end
```

- [ ] **Step 4: Add the segment-wise filtering path**

In `RELAX_filtbutter.m`, replace this block (currently lines 105-118):

```matlab
data = double(EEG.data);
temp = NaN(size(data,1),EEG.pnts,size(data,3));
for x = 1:size(data,1) 
    dataIn = data(x,:);
    if strcmpi(causal_or_acausal_filter,'causal')
        dataFilt1 = filter(z1,p1,dataIn);
    end
    if strcmpi(causal_or_acausal_filter,'acausal')
        dataFilt1 = filtfilt(z1,p1,dataIn);
    end
    temp(x,:,:) = reshape(dataFilt1,1,EEG.pnts,[]);
end 

EEG.data = temp;
```

with:

```matlab
data = double(EEG.data);

if respect_boundaries && size(data,3)==1
    % Filter each stretch of data between boundary events separately, with odd
    % reflection padding at each segment edge. The padding width is set from the
    % frequency that governs the filter's transient length: the low edge for a high
    % pass or band pass, the cutoff for a low pass, and the stop-band width for a
    % band stop (a notch rings at its bandwidth, not at its centre frequency).
    segments = RELAX_get_segments(EEG);

    if strcmp(type,'lowpass')
        f_governing = low;
    elseif strcmp(type,'stop')
        f_governing = abs(low - high);
    else % 'bandpass' or 'highpass'
        f_governing = high;
    end
    if isempty(f_governing) || f_governing <= 0
        f_governing = Fs/2;
    end
    padSamples = round(3*Fs/f_governing);

    if strcmpi(causal_or_acausal_filter,'causal')
        filterFn = @(d) filter(z1,p1,d,[],2);
    else
        filterFn = @(d) filtfilt(z1,p1,d')';
    end

    temp = RELAX_apply_by_segments(data, segments, filterFn, padSamples);
    fprintf('Filtering applied separately to %d segment(s) between boundary events (pad = %d samples).\n', ...
        size(segments,1), padSamples);
else
    temp = NaN(size(data,1),EEG.pnts,size(data,3));
    for x = 1:size(data,1) 
        dataIn = data(x,:);
        if strcmpi(causal_or_acausal_filter,'causal')
            dataFilt1 = filter(z1,p1,dataIn);
        end
        if strcmpi(causal_or_acausal_filter,'acausal')
            dataFilt1 = filtfilt(z1,p1,dataIn);
        end
        temp(x,:,:) = reshape(dataFilt1,1,EEG.pnts,[]);
    end 
end

EEG.data = temp;
```

Note the legacy loop is left byte-for-byte intact inside the `else`. Do not vectorise it: the backward-compatibility assertion depends on it producing exactly what it produced before. The vectorised `filtfilt(z1,p1,d')'` is used only on the new segment path, where there is no bitwise baseline to preserve.

- [ ] **Step 5: Run test to verify it passes**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: `5/5 test files passed.`

- [ ] **Step 6: Opt the wrapper's call sites in**

In `RELAX_Wrapper_slurm.m`, add `RELAX_cfg.RespectBoundaries` as a 7th argument to every `RELAX_filtbutter` call. There are ten, at lines 231, 233, 238, 240, 257, 263, 280, 285, 317 and 648. For example, line 257:

```matlab
            EEG = RELAX_filtbutter( EEG, RELAX_cfg.HighPassFilter, [], 4, 'highpass', RELAX_cfg.causal_or_acausal_filter);
```

becomes:

```matlab
            EEG = RELAX_filtbutter( EEG, RELAX_cfg.HighPassFilter, [], 4, 'highpass', RELAX_cfg.causal_or_acausal_filter, RELAX_cfg.RespectBoundaries);
```

and line 231:

```matlab
                EEG = RELAX_filtbutter( EEG, freq-4, freq+4, 4, 'bandstop','acausal');
```

becomes:

```matlab
                EEG = RELAX_filtbutter( EEG, freq-4, freq+4, 4, 'bandstop','acausal', RELAX_cfg.RespectBoundaries);
```

Apply the same transformation to all ten. Also add a default near the top of the per-file loop, immediately BEFORE the `RELAX_check_boundary_support` call added in Task 2, so the field exists by the time anything reads it:

```matlab
    if ~isfield(RELAX_cfg,'RespectBoundaries'); RELAX_cfg.RespectBoundaries = 0; end
```

Verify all ten were changed:

```bash
cd /home/nga4004/gitRepos/RELAX-GEDAI
grep -c "RELAX_filtbutter(.*RespectBoundaries" RELAX_Wrapper_slurm.m
```

Expected output: `10`

- [ ] **Step 7: Verify the wrapper still parses**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "checkcode('RELAX_Wrapper_slurm.m','-string')"
```

Expected: no syntax errors reported.

- [ ] **Step 8: Commit**

```bash
cd /home/nga4004/gitRepos/RELAX-GEDAI
git add RELAX_filtbutter.m RELAX_Wrapper_slurm.m tests/test_RELAX_filtbutter_boundaries.m
git commit -m "Make RELAX_filtbutter filter each inter-boundary segment separately"
```

---

### Task 4: Boundary-aware `RELAX_blinks_IQR_method`

**Files:**
- Create: `RELAX_blink_mask_bounds.m`
- Modify: `RELAX_blinks_IQR_method.m` (the two `RELAX_filtbutter` calls at lines 90 and 92; the blink-marking loop at lines 146-163)
- Test: `tests/test_RELAX_blink_mask_bounds.m`

**Interfaces:**
- Consumes: `RELAX_get_segments` from Task 1; the 7-argument `RELAX_filtbutter` from Task 3.
- Produces: `[maskStart, maskEnd] = RELAX_blink_mask_bounds(latency, halfWidth, segments, nPoints)` — returns the blink mask window centred on `latency`, extending `halfWidth` samples each way, truncated to the segment containing `latency`. Both outputs are integers, and `maskStart <= maskEnd` always holds.

**Behaviour on the boundary path:** every detected blink is marked, exactly as on the legacy path. A blink whose ±400 ms window would run past the edge of its segment has that window truncated to the segment, rather than being dropped. The mask and the `EyeBlinkLeftBase` / `EyeBlinkRightBase` events use the same truncated bounds so they stay consistent with each other; `EyeBlinkMax` is never moved. Truncating the base events is safe because nothing reads them — every consumer of blink events (`RELAX_excluding_channels_and_epoching.m:146`, `RELAX_excluding_extreme_values.m:74`, `RELAX_metrics_blinks.m:25`) epochs around `EyeBlinkMax` only.

Because no blink is ever skipped, `IQRmethodDetectedBlinks` keeps the value the detection stage gave it; the legacy path's rule of zeroing it when a blink sits near the record edge does not apply on the boundary path, and the legacy branch is left untouched.

- [ ] **Step 1: Write the failing test**

Create `tests/test_RELAX_blink_mask_bounds.m`:

```matlab
function test_RELAX_blink_mask_bounds()
% Blink mask windows, truncated to the segment containing the blink.

    segments = [1 4000; 4001 8000];
    nPoints  = 8000;
    hw       = 400; % half-width in samples

    % A blink well inside a segment is not truncated at all
    [s, e] = RELAX_blink_mask_bounds(2000, hw, segments, nPoints);
    assert(s == 1600 && e == 2400, 'interior blink must keep its full window');

    % A blink near the end of segment 1 is truncated on the right
    [s, e] = RELAX_blink_mask_bounds(3900, hw, segments, nPoints);
    assert(s == 3500 && e == 4000, 'window must be truncated at the segment end');

    % A blink near the start of segment 2 is truncated on the left
    [s, e] = RELAX_blink_mask_bounds(4100, hw, segments, nPoints);
    assert(s == 4001 && e == 4500, 'window must be truncated at the segment start');

    % Exactly on a segment start: belongs to that segment, truncated to it
    [s, e] = RELAX_blink_mask_bounds(4001, hw, segments, nPoints);
    assert(s == 4001 && e == 4401, 'blink on a segment start');

    % Exactly on a segment end: belongs to the earlier segment
    [s, e] = RELAX_blink_mask_bounds(4000, hw, segments, nPoints);
    assert(s == 3600 && e == 4000, 'blink on a segment end');

    % Truncation at the very start of the record, with one whole-record segment
    [s, e] = RELAX_blink_mask_bounds(50, hw, [1 nPoints], nPoints);
    assert(s == 1 && e == 450, 'window must be truncated at the record start');

    % Truncation at the very end of the record
    [s, e] = RELAX_blink_mask_bounds(nPoints-50, hw, [1 nPoints], nPoints);
    assert(s == nPoints-450 && e == nPoints, 'window must be truncated at the record end');

    % A latency outside every segment falls back to the whole record
    [s, e] = RELAX_blink_mask_bounds(99999, hw, segments, nPoints);
    assert(s == nPoints-hw && e == nPoints, 'out-of-range latency clamps to the record');

    % Degenerate case: a one-sample segment still yields a valid single-sample window
    [s, e] = RELAX_blink_mask_bounds(10, hw, [10 10], nPoints);
    assert(s == 10 && e == 10, 'one-sample segment yields a one-sample window');

    % Invariants that must hold for every case above
    testLatencies = [1 50 2000 3900 4000 4001 4100 7999 8000];
    for lat = testLatencies
        [s, e] = RELAX_blink_mask_bounds(lat, hw, segments, nPoints);
        assert(s <= e, 'maskStart must never exceed maskEnd (latency %d)', lat);
        assert(s >= 1 && e <= nPoints, 'window must stay inside the record (latency %d)', lat);
        assert(s == round(s) && e == round(e), 'bounds must be integers (latency %d)', lat);
        % The window must lie wholly within one segment
        inOne = any(segments(:,1) <= s & segments(:,2) >= e);
        assert(inOne, 'window must lie inside a single segment (latency %d)', lat);
    end

    % Non-integer latencies (the caller passes o+I from a max search) are rounded
    [s, e] = RELAX_blink_mask_bounds(2000.4, hw, segments, nPoints);
    assert(s == round(2000.4-hw) && e == round(2000.4+hw), 'non-integer latency handled');

    fprintf('test_RELAX_blink_mask_bounds: PASS\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: `test_RELAX_blink_mask_bounds` FAILs with "Unrecognized function or variable 'RELAX_blink_mask_bounds'".

- [ ] **Step 3: Write the implementation**

Create `RELAX_blink_mask_bounds.m`:

```matlab
% RELAX_blink_mask_bounds() - Blink mask window, truncated to the segment containing
%                             the blink.
%
% Usage:
%   >> [maskStart, maskEnd] = RELAX_blink_mask_bounds(latency, halfWidth, segments, nPoints);
%
% Inputs:
%   latency   - sample index of the blink peak (may be non-integer; it is rounded)
%   halfWidth - how many samples to mark either side of the peak
%   segments  - [nSegments x 2] matrix from RELAX_get_segments
%   nPoints   - total samples in the record, used as the fallback range
%
% Outputs:
%   maskStart, maskEnd - inclusive integer bounds of the blink window, truncated to
%                        the segment containing `latency`. maskStart <= maskEnd always.
%
% Every blink gets marked; a blink close to a discontinuity simply gets a shorter
% window rather than being discarded. Truncating rather than dropping matters here
% because the segments between boundary events are only ~8 s long, so a large share of
% blinks sit within 400 ms of a segment edge.

function [maskStart, maskEnd] = RELAX_blink_mask_bounds(latency, halfWidth, segments, nPoints)

    segStart = 1;
    segEnd   = nPoints;

    match = find(segments(:,1) <= latency & segments(:,2) >= latency, 1, 'first');
    if ~isempty(match)
        segStart = segments(match,1);
        segEnd   = segments(match,2);
    end

    maskStart = max(round(latency - halfWidth), segStart);
    maskEnd   = min(round(latency + halfWidth), segEnd);

    % A blink peak always lies inside its own segment, so this cannot invert; the
    % guard only covers the fallback case of a latency outside every segment.
    if maskStart > maskEnd
        maskStart = min(maskStart, segEnd);
        maskEnd   = maskStart;
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: `6/6 test files passed.`

- [ ] **Step 5: Pass the flag through the blink filter calls**

In `RELAX_blinks_IQR_method.m`, insert the following immediately above the line:

```matlab
    % Check that blink channels are still available after Noisy channels
```

(that line is the first statement after the `if exist('RELAX_cfg','var')==1 ... end` defaults block):

```matlab
    if isfield(RELAX_cfg,'RespectBoundaries')
        RespectBoundaries = RELAX_cfg.RespectBoundaries;
    else
        RespectBoundaries = 0;
    end
```

Then change lines 90 and 92 from:

```matlab
        EEGEyeOnly = RELAX_filtbutter( EEGEyeOnly, 1, 25, 4, 'bandpass','acausal');
    elseif strcmp(RELAX_cfg.LowPassFilterAt_6Hz_BeforeDetectingBlinks,'yes')
        EEGEyeOnly = RELAX_filtbutter( EEGEyeOnly, 1, 6, 4, 'bandpass','acausal');
```

to:

```matlab
        EEGEyeOnly = RELAX_filtbutter( EEGEyeOnly, 1, 25, 4, 'bandpass','acausal', RespectBoundaries);
    elseif strcmp(RELAX_cfg.LowPassFilterAt_6Hz_BeforeDetectingBlinks,'yes')
        EEGEyeOnly = RELAX_filtbutter( EEGEyeOnly, 1, 6, 4, 'bandpass','acausal', RespectBoundaries);
```

`EEGEyeOnly` is produced by `pop_select(EEG,'nochannel',...)`, which preserves `EEG.event`, so the boundary events are still present for `RELAX_get_segments` to find.

- [ ] **Step 6: Truncate the blink marking loop to segments**

In `RELAX_blinks_IQR_method.m`, replace this block (currently lines 146-163):

```matlab
    % Mark the blink peak into the event list:
    if continuousEEG.RELAX.IQRmethodDetectedBlinks==1
        for x=1:size(BlinkRunIndex,2) % v1.1.3 update, removed an erroneous "- 1" after size(BlinkRunIndex,2), thankyou Mana Biabani for the update
            if (round(BlinkMaxLatency(1,x)-400/RELAX_cfg.ms_per_sample)>0) && (round(BlinkMaxLatency(1,x)+400/RELAX_cfg.ms_per_sample))<size(EEG.data,2) %v1.1.6 update to address issue for sample rates above 1000Hz
                continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkLeftBase';
                continuousEEG.event(size(continuousEEG.event,2)).latency=BlinkMaxLatency(1,x)-(400/RELAX_cfg.ms_per_sample);
                continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkMax';
                continuousEEG.event(size(continuousEEG.event,2)).latency=BlinkMaxLatency(1,x);
                continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkRightBase';
                continuousEEG.event(size(continuousEEG.event,2)).latency=BlinkMaxLatency(1,x)+(400/RELAX_cfg.ms_per_sample);
                continuousEEG.RELAXProcessing.Details.eyeblinkmask(round(BlinkMaxLatency(1,x)-(400/RELAX_cfg.ms_per_sample)):round(BlinkMaxLatency(1,x)+(400/RELAX_cfg.ms_per_sample)))=1;
            else
                continuousEEG.RELAX.IQRmethodDetectedBlinks=0; % v1.1.3 update, so later functions do not search for blinks that aren't marked, if the blink was only present right on the edge of the file. Thankyou Mana Biabani for the update
            end
        end
    end    
```

with:

```matlab
    % Mark the blink peak into the event list:
    if continuousEEG.RELAX.IQRmethodDetectedBlinks==1
        if RespectBoundaries
            % Every blink is marked. A blink whose 400 ms window would run past the edge
            % of the stretch of continuous data it sits in has that window truncated to
            % the segment, so no mask ever spans a discontinuity and no blink is lost.
            % The base events use the same truncated bounds as the mask; nothing reads
            % them (every consumer epochs around EyeBlinkMax), but keeping them
            % consistent avoids events that point outside their own segment.
            blinkSegments = RELAX_get_segments(continuousEEG);
            blinkHalfWidth = 400/RELAX_cfg.ms_per_sample;
            for x=1:size(BlinkRunIndex,2)
                [maskStart, maskEnd] = RELAX_blink_mask_bounds(BlinkMaxLatency(1,x), ...
                    blinkHalfWidth, blinkSegments, size(EEG.data,2));
                continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkLeftBase';
                continuousEEG.event(size(continuousEEG.event,2)).latency=maskStart;
                continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkMax';
                continuousEEG.event(size(continuousEEG.event,2)).latency=BlinkMaxLatency(1,x);
                continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkRightBase';
                continuousEEG.event(size(continuousEEG.event,2)).latency=maskEnd;
                continuousEEG.RELAXProcessing.Details.eyeblinkmask(maskStart:maskEnd)=1;
            end
            fprintf('Marked %d blink(s), truncated to their segments where needed.\n', size(BlinkRunIndex,2));
        else
            for x=1:size(BlinkRunIndex,2) % v1.1.3 update, removed an erroneous "- 1" after size(BlinkRunIndex,2), thankyou Mana Biabani for the update
                if (round(BlinkMaxLatency(1,x)-400/RELAX_cfg.ms_per_sample)>0) && (round(BlinkMaxLatency(1,x)+400/RELAX_cfg.ms_per_sample))<size(EEG.data,2) %v1.1.6 update to address issue for sample rates above 1000Hz
                    continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkLeftBase';
                    continuousEEG.event(size(continuousEEG.event,2)).latency=BlinkMaxLatency(1,x)-(400/RELAX_cfg.ms_per_sample);
                    continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkMax';
                    continuousEEG.event(size(continuousEEG.event,2)).latency=BlinkMaxLatency(1,x);
                    continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkRightBase';
                    continuousEEG.event(size(continuousEEG.event,2)).latency=BlinkMaxLatency(1,x)+(400/RELAX_cfg.ms_per_sample);
                    continuousEEG.RELAXProcessing.Details.eyeblinkmask(round(BlinkMaxLatency(1,x)-(400/RELAX_cfg.ms_per_sample)):round(BlinkMaxLatency(1,x)+(400/RELAX_cfg.ms_per_sample)))=1;
                else
                    continuousEEG.RELAX.IQRmethodDetectedBlinks=0; % v1.1.3 update, so later functions do not search for blinks that aren't marked, if the blink was only present right on the edge of the file. Thankyou Mana Biabani for the update
                end
            end
        end
    end    
```

- [ ] **Step 7: Verify the file still parses and tests still pass**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); \
  checkcode('RELAX_blinks_IQR_method.m','-string'), run_all_tests"
```

Expected: no syntax errors, `6/6 test files passed.`

- [ ] **Step 8: Commit**

```bash
cd /home/nga4004/gitRepos/RELAX-GEDAI
git add RELAX_blinks_IQR_method.m RELAX_blink_mask_bounds.m tests/test_RELAX_blink_mask_bounds.m
git commit -m "Truncate blink masks to their segment instead of spanning boundaries"
```

---

### Task 5: Per-segment wavelet denoising helper

Extracts the wavelet thresholding currently inlined in `RELAX_targeted_wICA.m` into a standalone, testable function that operates on one block of data at a time. Task 6 wires it in.

**Files:**
- Create: `RELAX_wavelet_denoise_segment.m`
- Test: `tests/test_RELAX_wavelet_denoise_segment.m`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `comp_out = RELAX_wavelet_denoise_segment(comp_in, wavelet_threshold, threshold_type)` — `comp_in` is a `[1 x N]` row vector, `wavelet_threshold` a scalar, `threshold_type` the string returned by `ddencmp` (`'s'` or `'h'`). Returns `[1 x N]`, the same length as `comp_in`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_RELAX_wavelet_denoise_segment.m`:

```matlab
function test_RELAX_wavelet_denoise_segment()
% Level-5 coif5 stationary wavelet thresholding of one block of data.

    rng(42);
    n = 4000; % not a multiple of 32, so internal padding is exercised
    t = (0:n-1)/500;
    x = sin(2*pi*10*t) + 0.1*randn(1,n);

    % Length is always preserved, including for non-multiples of 32
    y = RELAX_wavelet_denoise_segment(x, 0.5, 's');
    assert(isequal(size(y), size(x)), 'output length must equal input length');
    assert(~any(isnan(y)), 'output must contain no NaNs');

    % A length that IS a multiple of 32 also works
    x32 = x(1:3968); % 3968 = 124*32
    y32 = RELAX_wavelet_denoise_segment(x32, 0.5, 's');
    assert(isequal(size(y32), size(x32)), 'output length preserved for multiple-of-32 input');

    % A zero threshold is a near-perfect reconstruction (swt/iswt round-trip)
    y0 = RELAX_wavelet_denoise_segment(x, 0, 's');
    assert(max(abs(y0 - x)) < 1e-8, 'zero threshold must reconstruct the input');

    % A larger threshold removes more energy than a smaller one
    y_small = RELAX_wavelet_denoise_segment(x, 0.1, 's');
    y_large = RELAX_wavelet_denoise_segment(x, 1.0, 's');
    assert(sum(y_large.^2) < sum(y_small.^2), 'a larger threshold must remove more energy');

    % Hard thresholding is also accepted
    yh = RELAX_wavelet_denoise_segment(x, 0.5, 'h');
    assert(isequal(size(yh), size(x)), 'hard threshold type accepted');

    % Determinism
    assert(isequal(RELAX_wavelet_denoise_segment(x, 0.5, 's'), y), 'must be deterministic');

    fprintf('test_RELAX_wavelet_denoise_segment: PASS\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: `test_RELAX_wavelet_denoise_segment` FAILs with "Unrecognized function or variable 'RELAX_wavelet_denoise_segment'".

- [ ] **Step 3: Write the implementation**

Create `RELAX_wavelet_denoise_segment.m`:

```matlab
% RELAX_wavelet_denoise_segment() - Level-5 coif5 stationary wavelet thresholding of
%                                   one contiguous block of data.
%
% Usage:
%   >> comp_out = RELAX_wavelet_denoise_segment(comp_in, wavelet_threshold, threshold_type);
%
% Inputs:
%   comp_in           - [1 x nSamples] row vector (one independent component)
%   wavelet_threshold - scalar threshold, as returned by ddencmp
%   threshold_type    - 's' (soft) or 'h' (hard), as returned by ddencmp
%
% Outputs:
%   comp_out - [1 x nSamples], the same length as comp_in
%
% swt requires the signal length to be divisible by 2^5, so the block is zero-padded
% up to the next multiple of 32 and trimmed back afterwards. The threshold is passed
% in rather than computed here so that a caller processing several segments of the
% same component can apply one threshold derived from the whole component, instead of
% a different threshold per segment.

function comp_out = RELAX_wavelet_denoise_segment(comp_in, wavelet_threshold, threshold_type)

    check_padding_required = mod(numel(comp_in), 2^5);
    if check_padding_required ~= 0
        padding = zeros(1,(2^5)-check_padding_required);
    else
        padding = [];
    end

    padded_comp = [comp_in, padding];

    wavelet_transform = swt(padded_comp,5,'coif5');
    thresholded_wavelet_transform = wthresh(wavelet_transform,threshold_type,wavelet_threshold);
    comp_out = iswt(thresholded_wavelet_transform,'coif5');

    if ~isempty(padding)
        comp_out = comp_out(1:end-numel(padding));
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: `7/7 test files passed.`

- [ ] **Step 5: Commit**

```bash
cd /home/nga4004/gitRepos/RELAX-GEDAI
git add RELAX_wavelet_denoise_segment.m tests/test_RELAX_wavelet_denoise_segment.m
git commit -m "Add per-segment wavelet denoising helper for wICA"
```

---

### Task 6: Boundary-aware `RELAX_targeted_wICA`

**Files:**
- Modify: `RELAX_targeted_wICA.m` (wavelet block at lines 238-279; blink `filtfilt` at lines 285-292; blink mask padding loops at lines 316-330; muscle `filtfilt` at lines 332-342)
- Test: `tests/test_RELAX_wica_mask_padding.m`

**Interfaces:**
- Consumes: `RELAX_get_segments`, `RELAX_apply_by_segments` from Task 1; `RELAX_wavelet_denoise_segment` from Task 5.
- Produces: `RELAX_pad_blink_mask_block(mask, moving_mean_length)` in its own file, so the mask-padding logic is testable in isolation.

- [ ] **Step 1: Write the failing test**

Create `tests/test_RELAX_wica_mask_padding.m`:

```matlab
function test_RELAX_wica_mask_padding()
% Blink mask dilation and smoothing, applied to one block at a time.

    % A single flagged sample must dilate symmetrically by moving_mean_length
    mask = zeros(1,101);
    mask(51) = 1;
    out = RELAX_pad_blink_mask_block(mask, 10);
    assert(isequal(size(out), size(mask)), 'output size must match input size');
    assert(all(out(41:61) > 0), 'the dilated region must be non-zero');
    assert(out(51) > 0, 'the flagged sample itself must stay non-zero');
    assert(out(1) == 0 && out(101) == 0, 'far-away samples must stay zero');

    % An all-zero mask stays all zero
    assert(all(RELAX_pad_blink_mask_block(zeros(1,101), 10) == 0), 'zero mask stays zero');

    % An all-ones mask stays all ones
    assert(all(abs(RELAX_pad_blink_mask_block(ones(1,101), 10) - 1) < 1e-12), 'ones mask stays ones');

    % ISOLATION through RELAX_apply_by_segments: a blink at the end of segment 1
    % must not bleed into segment 2.
    segments = [1 100; 101 200];
    mask = zeros(1,200);
    mask(98) = 1; % 2 samples from the end of segment 1
    out = RELAX_apply_by_segments(mask, segments, @(d) RELAX_pad_blink_mask_block(d, 10), 0);
    assert(all(out(101:200) == 0), 'dilation must not cross a segment boundary');
    assert(any(out(88:100) > 0), 'dilation must still occur inside segment 1');

    fprintf('test_RELAX_wica_mask_padding: PASS\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: `test_RELAX_wica_mask_padding` FAILs with "Unrecognized function or variable 'RELAX_pad_blink_mask_block'".

- [ ] **Step 3: Extract the mask padding into its own function**

Create `RELAX_pad_blink_mask_block.m`:

```matlab
% RELAX_pad_blink_mask_block() - Dilate and smooth a binary blink mask over one block.
%
% Usage:
%   >> out = RELAX_pad_blink_mask_block(mask, moving_mean_length);
%
% Inputs:
%   mask               - [1 x nSamples] binary mask of blink periods
%   moving_mean_length - dilation radius in samples
%
% Outputs:
%   out - [1 x nSamples] mask dilated by moving_mean_length in each direction and then
%         smoothed with a centred moving average of the same half-width, giving the
%         tapered weighting that wICA multiplies blink components by.
%
% Applied one block at a time so that, when called through RELAX_apply_by_segments, a
% blink near the end of one segment cannot dilate into the next.

function out = RELAX_pad_blink_mask_block(mask, moving_mean_length)

    padded_blink_periods = double(mask);
    n = size(padded_blink_periods,2);

    for c = flip(1:max(0, n-(moving_mean_length+1)))
        if padded_blink_periods(1,c)==1
            padded_blink_periods(1,c:c+moving_mean_length)=1;
        end
    end
    for c = (moving_mean_length+1):n
        if padded_blink_periods(1,c)==1
            padded_blink_periods(1,c-moving_mean_length:c)=1;
        end
    end

    out = movmean(padded_blink_periods,[moving_mean_length moving_mean_length]);
end
```

- [ ] **Step 4: Run test to verify it passes**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: `8/8 test files passed.`

- [ ] **Step 5: Add the flag and segment list to `RELAX_targeted_wICA`**

In `RELAX_targeted_wICA.m`, immediately before the line:

```matlab
    check_padding_required = mod(size(Component,2),2^5);
```

insert:

```matlab
    % Boundary handling. Derive the segment list from the CURRENT EEG: eeg_eegrej runs
    % earlier in the wrapper and adds boundaries, so wICA sees more discontinuities than
    % the filtering stage did.
    if isfield(RELAX_cfg,'RespectBoundaries')
        RespectBoundaries = RELAX_cfg.RespectBoundaries;
    else
        RespectBoundaries = 0;
    end
    wICAsegments = RELAX_get_segments(EEG);
    wavelet_pad_samples = 512; % coif5 level-5 support is about 30*2^4 = 480 samples
```

- [ ] **Step 6: Route the wavelet thresholding through the segment helper**

In `RELAX_targeted_wICA.m`, replace this block (currently lines 247-265):

```matlab
    for comp_number = 1:size(Component,1)
        if ICsMostLikelyNotBrain(comp_number)==1 % wavelet enhance only on artifacts identified by ICLabel
            if ~isempty(padding)
                padded_comp = [Component(comp_number,:),padding]; % pad the component with zeros if required
            else
                padded_comp = Component(comp_number,:);
            end
            [wavelet_threshold,threshold_type,~] = ddencmp('den','wv',padded_comp); % automatically obtain wavelet enhancement threshold
            if ICsMostLikelyEye(comp_number)==1
                wavelet_threshold = wavelet_threshold*2; % increase threshold for blink components based on optimal results in our informal testing
            else
                wavelet_threshold = wavelet_threshold*1;
            end
            wavelet_transform = swt(padded_comp,5,'coif5'); % apply stationary wavelet transform to each component to reduce neural contribution to component
            thresholded_wavelet_transform = wthresh(wavelet_transform,threshold_type,wavelet_threshold); % remove negligible values by applying thresholding
            artifact_comp(comp_number,:) = iswt(thresholded_wavelet_transform,'coif5'); % use inverse wavelet transform to obtained the wavelet transformed component
            clear thresholded_wavelet_transform padded_comp wavelet_threshold threshold_type wavelet_transform 
        end
    end
```

with:

```matlab
    for comp_number = 1:size(Component,1)
        if ICsMostLikelyNotBrain(comp_number)==1 % wavelet enhance only on artifacts identified by ICLabel
            if ~isempty(padding)
                padded_comp = [Component(comp_number,:),padding]; % pad the component with zeros if required
            else
                padded_comp = Component(comp_number,:);
            end
            [wavelet_threshold,threshold_type,~] = ddencmp('den','wv',padded_comp); % automatically obtain wavelet enhancement threshold
            if ICsMostLikelyEye(comp_number)==1
                wavelet_threshold = wavelet_threshold*2; % increase threshold for blink components based on optimal results in our informal testing
            else
                wavelet_threshold = wavelet_threshold*1;
            end
            if RespectBoundaries
                % Transform each segment separately. The threshold stays the one derived
                % from the whole component above, so segmentation changes only where the
                % wavelet basis is applied, not how aggressively it thresholds.
                artifact_comp(comp_number,:) = RELAX_apply_by_segments(Component(comp_number,:), wICAsegments, ...
                    @(d) RELAX_wavelet_denoise_segment(d, wavelet_threshold, threshold_type), wavelet_pad_samples);
            else
                wavelet_transform = swt(padded_comp,5,'coif5'); % apply stationary wavelet transform to each component to reduce neural contribution to component
                thresholded_wavelet_transform = wthresh(wavelet_transform,threshold_type,wavelet_threshold); % remove negligible values by applying thresholding
                artifact_comp(comp_number,:) = iswt(thresholded_wavelet_transform,'coif5'); % use inverse wavelet transform to obtained the wavelet transformed component
                clear thresholded_wavelet_transform wavelet_transform
            end
            clear padded_comp wavelet_threshold threshold_type
        end
    end
```

Note: on the boundary path `artifact_comp(comp_number,:)` is assigned an unpadded row, whereas the legacy path assigns a padded one. The existing "remove padding" step further down (`artifact_comp = artifact_comp(:,1:end-numel(padding));`) must therefore only run on the legacy path. Change:

```matlab
    % remove padding
    if ~isempty(padding)
        artifact_comp = artifact_comp(:,1:end-numel(padding));
    end
```

to:

```matlab
    % remove padding (only the legacy whole-record path adds any)
    if ~isempty(padding) && ~RespectBoundaries
        artifact_comp = artifact_comp(:,1:end-numel(padding));
    end
```

and change the zero-fill for the "no artifact components" case:

```matlab
    if sum(ICsMostLikelyNotBrain)==0
        artifact_comp(1,:)=zeros(1,size(EEG.data,2));
        artifact_comp = [artifact_comp(:,:),padding]; % pad with zeros
    end
```

to:

```matlab
    if sum(ICsMostLikelyNotBrain)==0
        artifact_comp(1,:)=zeros(1,size(EEG.data,2));
        if ~RespectBoundaries
            artifact_comp = [artifact_comp(:,:),padding]; % pad with zeros
        end
    end
```

- [ ] **Step 7: Route the blink localisation filter through the segment helper**

In `RELAX_targeted_wICA.m`, replace:

```matlab
            [z1, p1] = butter(2, [0.5 25]./(EEG.srate/2), 'bandpass');
            dataIn=Component(comp_number,:)';
            dataIn=double(dataIn);
            dataFilt1 = filtfilt(z1,p1,dataIn);
            IC_filtered = dataFilt1';
```

with:

```matlab
            [z1, p1] = butter(2, [0.5 25]./(EEG.srate/2), 'bandpass');
            if RespectBoundaries
                % Pad width from the 0.5 Hz low edge, which governs the transient length
                IC_filtered = RELAX_apply_by_segments(double(Component(comp_number,:)), wICAsegments, ...
                    @(d) filtfilt(z1,p1,d')', round(3*EEG.srate/0.5));
            else
                dataIn=Component(comp_number,:)';
                dataIn=double(dataIn);
                dataFilt1 = filtfilt(z1,p1,dataIn);
                IC_filtered = dataFilt1';
            end
```

- [ ] **Step 8: Route the blink mask dilation through the segment helper**

In `RELAX_targeted_wICA.m`, replace:

```matlab
            padded_blink_periods=double(blink_periods);
            for c=flip(1:size(padded_blink_periods,2)-(moving_mean_length+1))
                if padded_blink_periods(1,c)==1 
                    padded_blink_periods(1,c:c+moving_mean_length)=1;
                end
            end
            for c=(moving_mean_length+1):size(padded_blink_periods,2)
                if padded_blink_periods(1,c)==1 
                    padded_blink_periods(1,c-moving_mean_length:c)=1;
                end
            end
            M(comp_number,:) = movmean(padded_blink_periods,[moving_mean_length moving_mean_length]);
```

with:

```matlab
            if RespectBoundaries
                % Pad width 0: this is a binary mask, so reflection padding is meaningless;
                % the point is only that dilation cannot cross a segment edge.
                M(comp_number,:) = RELAX_apply_by_segments(double(blink_periods), wICAsegments, ...
                    @(d) RELAX_pad_blink_mask_block(d, moving_mean_length), 0);
            else
                M(comp_number,:) = RELAX_pad_blink_mask_block(double(blink_periods), moving_mean_length);
            end
```

`RELAX_pad_blink_mask_block` contains exactly the loops and `movmean` being replaced, so the `else` branch is behaviour-preserving.

- [ ] **Step 9: Route the muscle highpass through the segment helper**

In `RELAX_targeted_wICA.m`, replace:

```matlab
        if ICsMostLikelyMuscle(comp_number)==1
            [z1, p1] = butter(2, 15./(EEG.srate/2), 'high');
            dataIn=Component(comp_number,:)';
            dataIn=double(dataIn);
            dataFilt1 = filtfilt(z1,p1,dataIn);
            artifact_comp(comp_number,:) = dataFilt1';
        end
```

with:

```matlab
        if ICsMostLikelyMuscle(comp_number)==1
            [z1, p1] = butter(2, 15./(EEG.srate/2), 'high');
            if RespectBoundaries
                % Pad width from the 15 Hz cutoff, which governs the transient length
                artifact_comp(comp_number,:) = RELAX_apply_by_segments(double(Component(comp_number,:)), wICAsegments, ...
                    @(d) filtfilt(z1,p1,d')', round(3*EEG.srate/15));
            else
                dataIn=Component(comp_number,:)';
                dataIn=double(dataIn);
                dataFilt1 = filtfilt(z1,p1,dataIn);
                artifact_comp(comp_number,:) = dataFilt1';
            end
        end
```

- [ ] **Step 10: Verify the file parses and all tests pass**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); \
  checkcode('RELAX_targeted_wICA.m','-string'), run_all_tests"
```

Expected: no syntax errors, `8/8 test files passed.`

- [ ] **Step 11: Commit**

```bash
cd /home/nga4004/gitRepos/RELAX-GEDAI
git add RELAX_targeted_wICA.m RELAX_pad_blink_mask_block.m tests/test_RELAX_wica_mask_padding.m
git commit -m "Apply wICA wavelets, filters and blink masks per inter-boundary segment"
```

---

### Task 7: Real-data fixture, epoching verification, end-to-end run, and documentation

**Files:**
- Create: `tests/make_boundary_fixture.m`
- Create: `tests/verify_epoching_drops_boundary_epochs.m`
- Create: `tests/run_integration_smoke.sh`
- Modify: `CLAUDE.md` (add a "Boundary/discontinuity handling" section)

**Interfaces:**
- Consumes: everything from Tasks 1-6.
- Produces: `tests/fixtures/boundary_fixture.set`, a cropped real recording used by the epoching verification and the smoke run.

- [ ] **Step 1: Write the fixture builder**

Create `tests/make_boundary_fixture.m`:

```matlab
function fixturePath = make_boundary_fixture()
% Build a small cropped real recording that contains boundary events, for the
% integration checks. Crops around the first boundary of a real input file so the
% fixture spans a dozen or so of the ~8 s inter-train segments.
%
% Returns the path to the written .set file.

    sourceDir  = '/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Aug26/tms_artifact_removed/';
    sourceFile = 'm121_dlpfc_day1_treatment_20221226_034250_interTrainIntervals.set';

    fixtureDir = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
    if ~exist(fixtureDir, 'dir'); mkdir(fixtureDir); end
    fixturePath = fullfile(fixtureDir, 'boundary_fixture.set');

    if exist(fixturePath, 'file')
        fprintf('Fixture already exists at %s\n', fixturePath);
        return
    end

    EEG = pop_loadset('filename', sourceFile, 'filepath', sourceDir);

    segments = RELAX_get_segments(EEG);
    assert(size(segments,1) > 2, 'source file must contain boundary events');

    % Start 5 s before the first boundary and keep 120 s, which spans the first
    % boundary and roughly a dozen of the short inter-train segments after it.
    firstBoundary = segments(1,2);
    startSample = max(1, firstBoundary - 5*EEG.srate);
    endSample   = min(EEG.pnts, startSample + 120*EEG.srate - 1);

    EEG = pop_select(EEG, 'point', [startSample endSample]);
    EEG = eeg_checkset(EEG);

    croppedSegments = RELAX_get_segments(EEG);
    fprintf('Fixture: %d samples, %d segments, %d boundaries.\n', ...
        EEG.pnts, size(croppedSegments,1), size(croppedSegments,1)-1);
    assert(size(croppedSegments,1) > 5, 'cropped fixture must retain several segments');

    EEG = pop_saveset(EEG, 'filename', 'boundary_fixture.set', 'filepath', fixtureDir);
    fprintf('Wrote fixture to %s\n', fixturePath);
end
```

- [ ] **Step 2: Build the fixture and confirm it has boundaries**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
srun -p scu-cpu -c 2 --mem 32G -t 00:30:00 \
  matlab -nodisplay -nosplash -batch \
  "addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; \
   addpath(pwd); addpath(fullfile(pwd,'tests')); make_boundary_fixture"
```

Expected: prints a segment count greater than 5 and writes `tests/fixtures/boundary_fixture.set`. If `eeglab nogui` errors, use `eeglab('nogui')`.

- [ ] **Step 3: Write the epoching verification**

This is spec Step 5 — it checks the claim that `pop_epoch` already protects the epoch-based stages, rather than changing any code.

Create `tests/verify_epoching_drops_boundary_epochs.m`:

```matlab
function verify_epoching_drops_boundary_epochs()
% Verify the assumption that pop_epoch discards every epoch containing a boundary
% event, which is what protects RELAX_epoching, RELAX_excluding_channels_and_epoching,
% RELAX_excluding_extreme_values, RELAX_muscle and RELAX_drift from discontinuities.
%
% Requires EEGLAB on the path and tests/fixtures/boundary_fixture.set to exist.

    fixturePath = fullfile(fileparts(mfilename('fullpath')), 'fixtures', 'boundary_fixture.set');
    assert(exist(fixturePath,'file') == 2, 'run make_boundary_fixture first');

    EEG = pop_loadset(fixturePath);
    segments = RELAX_get_segments(EEG);
    fprintf('Fixture has %d segments.\n', size(segments,1));
    assert(size(segments,1) > 1, 'fixture must contain boundaries for this check to mean anything');

    % Reproduce the epoching that RELAX_epoching performs
    EEG = eeg_regepochs(EEG,'recurrence',0.500,'eventtype','X','extractepochs','off');
    for e = 1:numel(EEG.event)
        EEG.event(e).originallatency = EEG.event(e).latency;
    end
    EEG = pop_selectevent(EEG, 'type', 'X', 'deleteevents', 'off');
    epoched = pop_epoch(EEG, {'X'}, [0 1.0], 'epochinfo', 'yes');

    % No surviving epoch may contain a boundary event
    boundaryInEpoch = 0;
    for ep = 1:numel(epoched.epoch)
        types = epoched.epoch(ep).eventtype;
        if ~iscell(types); types = {types}; end
        if any(strcmp(types, 'boundary'))
            boundaryInEpoch = boundaryInEpoch + 1;
        end
    end
    fprintf('Epochs retained: %d. Epochs containing a boundary: %d.\n', epoched.trials, boundaryInEpoch);
    assert(boundaryInEpoch == 0, ...
        'pop_epoch retained %d epoch(s) containing a boundary - the epoch-based stages are NOT protected and need their own fix', ...
        boundaryInEpoch);

    % Every retained epoch must lie wholly inside one segment
    epochLen = round(1.0 * epoched.srate);
    outside = 0;
    for ep = 1:numel(epoched.epoch)
        lat = epoched.epoch(ep).eventoriginallatency;
        if iscell(lat); lat = lat{1}; end
        inSegment = any(segments(:,1) <= lat & segments(:,2) >= lat + epochLen - 1);
        if ~inSegment
            outside = outside + 1;
        end
    end
    fprintf('Epochs not wholly inside a single segment: %d.\n', outside);
    assert(outside == 0, '%d epoch(s) straddle a segment boundary', outside);

    fprintf('verify_epoching_drops_boundary_epochs: PASS\n');
end
```

- [ ] **Step 4: Run the epoching verification**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
srun -p scu-cpu -c 2 --mem 32G -t 00:30:00 \
  matlab -nodisplay -nosplash -batch \
  "addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; \
   addpath(pwd); addpath(fullfile(pwd,'tests')); verify_epoching_drops_boundary_epochs"
```

Expected: `verify_epoching_drops_boundary_epochs: PASS`.

**If this FAILS**, stop and report it: the spec's assumption that the epoch-based stages are already protected is wrong, and `RELAX_epoching` / `RELAX_excluding_channels_and_epoching` need their own boundary handling as an additional task. Do not attempt that fix as part of this task.

- [ ] **Step 5: Write and run the end-to-end smoke test**

This exercises the wrapper edits from Task 2 (the PREP and SER/ARR gating) and the whole modified path together. It validates plumbing, not cleaning quality: 120 s of data is far too little for a well-conditioned 257-channel ICA, and that is fine — the check is that the pipeline runs to completion and produces the expected fields.

Create `tests/run_integration_smoke.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=relax-boundary-smoke
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=04:00:00
#SBATCH --partition=scu-cpu
#SBATCH --output=logs/boundary_smoke_%j.log
#SBATCH --error=logs/boundary_smoke_%j.error

# End-to-end smoke run of the boundary-aware pipeline on the cropped fixture.
# Submit with: sbatch tests/run_integration_smoke.sh

module load matlab

cd /home/nga4004/gitRepos/RELAX-GEDAI

matlab -nodisplay -nosplash -batch "\
addpath('/home/nga4004/tools/eeglab2026.0.0'); eeglab nogui; \
addpath('/home/nga4004/tools/eeglab2026.0.0/plugins/Fieldtrip-lite250523/'); \
addpath(pwd); addpath(fullfile(pwd,'tests')); \
run_integration_smoke"
```

Create `tests/run_integration_smoke.m`:

```matlab
function run_integration_smoke()
% Run RELAX_Wrapper_slurm end to end on the cropped fixture and assert the
% boundary-related outcomes. Validates plumbing, not cleaning quality.

    fixtureDir = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
    assert(exist(fullfile(fixtureDir,'boundary_fixture.set'),'file') == 2, ...
        'run make_boundary_fixture first');

    RELAX_cfg = struct();
    RELAX_cfg.caploc = [];
    RELAX_cfg.myPath = fixtureDir;
    RELAX_cfg.all_data_in_1_folder_or_BIDS_format = 'folder';
    RELAX_cfg.filename = [];
    RELAX_cfg.files = {'boundary_fixture.set'};
    RELAX_cfg.dirList = dir(fullfile(fixtureDir,'boundary_fixture.set'));
    RELAX_cfg.FilesToProcess = 1;
    RELAX_cfg.ParticipantTag = 'smoke';

    % The three new flags under test
    RELAX_cfg.RespectBoundaries = 1;
    RELAX_cfg.Perform_PREP_bad_channel_rejection = 0;
    RELAX_cfg.ComputeSER_ARR = 0;

    % Everything else copied from RELAX_SET_PARAMETERS_AND_RUN_slurm.m
    RELAX_cfg.Perform_targeted_wICA = 1;
    RELAX_cfg.Run_GEDAI = 0;
    RELAX_cfg.Do_MWF_Once = 0; RELAX_cfg.Do_MWF_Twice = 0; RELAX_cfg.Do_MWF_Thrice = 0;
    RELAX_cfg.Perform_wICA_on_ICLabel = 0;
    RELAX_cfg.Perform_ICA_subtract = 0;
    RELAX_cfg.ICA_method = 'picard';
    RELAX_cfg.Report_all_ICA_info = 'yes';
    RELAX_cfg.Clean_other_comps = 'no';
    RELAX_cfg.ICLabel_thresholds = [0.5 0.5 0.5 0.5 0.5 0.5 0.5];
    RELAX_cfg.computerawmetrics = 1;
    RELAX_cfg.computecleanedmetrics = 1;
    RELAX_cfg.MWFRoundToCleanBlinks = 2;
    RELAX_cfg.LowPassFilterAt_6Hz_BeforeDetectingBlinks = 'no';
    RELAX_cfg.ProbabilityDataHasNoBlinks = 0;
    RELAX_cfg.DriftSeverityThreshold = 10;
    RELAX_cfg.ProportionWorstEpochsForDrift = 0.30;
    RELAX_cfg.ExtremeVoltageShiftThreshold = 25;
    RELAX_cfg.ExtremeAbsoluteVoltageThreshold = 1000;
    RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold = 10;
    RELAX_cfg.ExtremeSingleChannelKurtosisThreshold = 10;
    RELAX_cfg.ExtremeAllChannelKurtosisThreshold = 10;
    RELAX_cfg.ExtremeDriftSlopeThreshold = -4;
    RELAX_cfg.ExtremeBlinkShiftThreshold = 10;
    RELAX_cfg.MinimumArtifactDuration = 1200;
    RELAX_cfg.MinimumBlinkArtifactDuration = 800;
    RELAX_cfg.BlinkElectrodes = {'E18'; 'E37'; 'E19'; 'E33'; 'E32'; 'E25'; 'E31'; 'E26'};
    RELAX_cfg.HEOGLeftpattern  = ["E46", "E54", "E252", "E47", "E248", "E241", "E244"];
    RELAX_cfg.HEOGRightpattern = ["E10", "E1", "E226", "E2", "E230", "E238", "E234"];
    RELAX_cfg.BlinkMaskFocus = 150;
    RELAX_cfg.HorizontalEyeMovementType = 2;
    RELAX_cfg.HorizontalEyeMovementThreshold = 2;
    RELAX_cfg.HorizontalEyeMovementThresholdIQR = 1.5;
    RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold = 25;
    RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow = 49;
    RELAX_cfg.HorizontalEyeMovementFocus = 200;
    RELAX_cfg.LowPassFilterBeforeMWF = 'yes';
    RELAX_cfg.DownSample = 'yes';
    RELAX_cfg.DownSample_to_X_Hz = 500;
    RELAX_cfg.FilterType = 'Butterworth';
    RELAX_cfg.causal_or_acausal_filter = 'acausal';
    RELAX_cfg.HighPassFilter = 1;
    RELAX_cfg.LowPassFilter = 200;
    RELAX_cfg.NotchFilterType = 'ZaplinePlus';
    RELAX_cfg.LineNoiseFrequency = 60;
    RELAX_cfg.ElectrodesToDelete = {'E257'};
    RELAX_cfg.electrodes_2_keep_but_not_clean = {''};
    RELAX_cfg.HighPassFilter_aux_elecs = 1;
    RELAX_cfg.LowPassFilter_aux_elecs = 200;
    RELAX_cfg.KeepAllInfo = 0;
    RELAX_cfg.saveextremesrejected = 1;
    RELAX_cfg.saveround1 = 0; RELAX_cfg.saveround2 = 0; RELAX_cfg.saveround3 = 0;
    RELAX_cfg.OnlyIncludeTaskRelatedEpochs = 0;
    RELAX_cfg.MuscleSlopeThreshold = -0.31;
    RELAX_cfg.MaxProportionOfDataCanBeMarkedAsMuscle = 0.50;
    RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel = 0.50;
    RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel = 0.25;
    RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0.25;
    RELAX_cfg.InterpolateRejectedElectrodesAfterCleaning = 'yes';
    RELAX_cfg.MWFDelayPeriod_for_eye_movements = 10;
    RELAX_cfg.MWFDelayPeriod_for_muscle_artifacts = 10;
    RELAX_cfg.MWF_delay_spacing_for_eye_movements = 16;
    RELAX_cfg.MWF_delay_spacing_for_muscle_artifacts = 2;

    RELAX_Wrapper_slurm(RELAX_cfg);

    % --- Assertions on the output --------------------------------------------------
    cleanedDir = fullfile(fixtureDir, 'RELAXProcessed', 'Cleaned_Data');
    cleaned = dir(fullfile(cleanedDir, '*RELAX*.set'));
    assert(~isempty(cleaned), 'no cleaned file was written to %s', cleanedDir);

    OUT = pop_loadset('filename', cleaned(1).name, 'filepath', cleanedDir);

    % PREP was skipped, so its channel list must be empty
    assert(isfield(OUT,'RELAXProcessingExtremeRejections'), 'missing RELAXProcessingExtremeRejections');
    assert(isempty(OUT.RELAXProcessingExtremeRejections.PREPBasedChannelToReject), ...
        'PREPBasedChannelToReject must be empty when PREP is skipped');

    % The SER/ARR block was skipped, so its only product must be absent
    assert(~isfield(OUT.RELAX,'NoiseMaskFullLengthR1'), ...
        'NoiseMaskFullLengthR1 must be absent when ComputeSER_ARR=0');

    % Boundaries must have survived cleaning
    outSegments = RELAX_get_segments(OUT);
    fprintf('Cleaned file has %d segments (%d boundaries).\n', ...
        size(outSegments,1), size(outSegments,1)-1);
    assert(size(outSegments,1) > 1, 'cleaned file lost all its boundary events');

    assert(OUT.RELAX.Data_has_been_cleaned == 1, 'Data_has_been_cleaned must be 1');
    assert(~any(isnan(OUT.data(:))), 'cleaned data must contain no NaNs');

    fprintf('Retained channels: %d\n', OUT.nbchan);
    fprintf('run_integration_smoke: PASS\n');
end
```

Submit it:

```bash
cd /home/nga4004/gitRepos/RELAX-GEDAI
mkdir -p logs
sbatch tests/run_integration_smoke.sh
```

Poll with `squeue -u $USER` until the job leaves the queue, then read `logs/boundary_smoke_<jobid>.log`.

Expected: `run_integration_smoke: PASS`, with the log showing "Skipping PREP findNoisyChannels", "Filtering applied separately to N segment(s)" for each filter call, and a non-zero retained-channel count. **Record the retained-channel count in your report** — it is the measure of how much less aggressive channel rejection is now that PREP is skipped.

- [ ] **Step 6: Document the behaviour in `CLAUDE.md`**

Append this section to `CLAUDE.md`, immediately after the "Skip-already-processed behavior" section:

```markdown
### Boundary/discontinuity handling

The TMS-EEG inputs to `RELAX_Wrapper_slurm.m` contain `boundary` events marking periods excised by
TMS artifact removal, leaving only ~8 s of continuous data between them. `RELAX_cfg.RespectBoundaries=1`
makes the pipeline process each inter-boundary segment independently, with odd reflection padding at
each segment edge. Design and analysis:
`docs/superpowers/specs/2026-09-08-boundary-aware-preprocessing-design.md`.

Helpers: `RELAX_get_segments.m` (segment list from `EEG.event`), `RELAX_reflect_pad.m` (odd reflection),
`RELAX_apply_by_segments.m` (apply a function per segment), `RELAX_check_boundary_support.m` (startup
guard). Callers must always re-derive the segment list from the *current* `EEG` — `eeg_eegrej` adds
boundaries part-way through the wrapper.

Boundary-aware without any change from us (verified against EEGLAB 2026.0.0 source): `pop_resample`
(splits at boundaries), `clean_data_with_zapline_plus_eeglab_wrapper` (passes `boundaries` to zapline),
`pop_eegfiltnew`/`firfilt` (splits at discontinuities), and `pop_epoch` (drops any epoch containing a
boundary — which is what protects `RELAX_epoching`, `RELAX_excluding_channels_and_epoching`,
`RELAX_excluding_extreme_values`, `RELAX_muscle` and `RELAX_drift`). `RELAX_ICA_subtract` is also safe:
it is a purely spatial projection.

Made boundary-aware here: `RELAX_filtbutter.m` (optional 7th argument), `RELAX_blinks_IQR_method.m`,
`RELAX_targeted_wICA.m`.

**Still boundary-blind — `RELAX_check_boundary_support` errors if any of these is enabled on data with
boundaries:** PREP `findNoisyChannels` (`Perform_PREP_bad_channel_rejection`, off by default — it tiles
correlation windows flatly from sample 1 against a 1% `badTimeThreshold`, so ~9-11% of windows straddle
a discontinuity on this data), MWF delay-stacking (`Do_MWF_*`), `Run_GEDAI`, the Butterworth notch,
`Perform_wICA_on_ICLabel`, and the SER/ARR marking block (`ComputeSER_ARR`, off by default — its only
consumer `RELAX_metrics_final_SER_and_ARR` is commented out, and `RELAX_muscle`, `RELAX_horizontaleye`,
`RELAX_drift` and `RELAX_pad_brief_mask_periods` would all need fixing first).

Tests live in `tests/`. Run them with:
`matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"`
```

- [ ] **Step 7: Run the full unit suite one last time**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: `8/8 test files passed.`

- [ ] **Step 8: Commit**

Do not commit the fixture `.set`/`.fdt` — they are large binaries derived from real recordings.

```bash
cd /home/nga4004/gitRepos/RELAX-GEDAI
echo "tests/fixtures/" >> .gitignore
git add .gitignore CLAUDE.md tests/make_boundary_fixture.m \
        tests/verify_epoching_drops_boundary_epochs.m \
        tests/run_integration_smoke.m tests/run_integration_smoke.sh
git commit -m "Add integration fixture, epoching verification and boundary handling docs"
```
