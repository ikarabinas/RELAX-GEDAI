# Electrode-Pop Detection and Repair Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a pre-filtering stage to `RELAX_Wrapper_slurm.m` that detects electrode pops (abrupt DC baseline steps) and repairs them by robust per-sub-segment detrending, before the 1 Hz Butterworth highpass can convert each step into multi-second bipolar ringing.

**Architecture:** Two new functions behind a stable `steps` struct - `RELAX_detect_steps` (median-filtered derivative + robust z + rise-time criterion) and `RELAX_repair_steps` (partition each channel at {boundary starts} U {its own step samples}, then fit and remove a low-order polynomial inside each sub-segment independently). Because the fitted trend carries the weighted mean, per-sub-segment detrending removes each sub-segment's offset as a side effect, so no separate step-subtraction pass exists. Chronically-popping channels escalate to rejection through a new shared channel-rejection budget.

**Tech Stack:** MATLAB (no package manager, no build). EEGLAB for the `EEG` struct. NoiseTools `nt_detrend` vendored into the repo. Tests are plain MATLAB functions using `assert`, discovered by `tests/run_all_tests.m`.

**Spec:** `docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md`

## Global Constraints

- **MATLAB only.** No build step, no package manager, no test framework beyond `tests/run_all_tests.m`.
- **Run the full suite with:** `matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"` from the repo root. Load MATLAB first with `module load matlab`.
- **Run one test with:** `matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_NAME"`
- **Never run MATLAB on a login node** (hostnames starting `scu-login`). Use `srun`/`sbatch` on the `sackler` or `scu-cpu` partition for anything touching real data. Unit tests on synthetic arrays are small enough to run in a `-batch` call directly.
- **Test file convention:** a file `tests/test_X.m` defines `function test_X()`, uses bare `assert(cond, 'message')`, and ends with `fprintf('test_X: PASS\n');`. It takes no arguments and returns nothing. Failure = a thrown error.
- **Only `RELAX_Wrapper_slurm.m` is modified.** Sibling wrappers (`RELAX_Wrapper.m`, `RELAX_Wrapper_bart.m`) must keep byte-identical behaviour. Any change to a shared function must be backward-compatible by default.
- **Absent config field semantics:** `RELAX_cfg.Perform_pop_artifact_repair` absent means **off** (new behaviour must not switch itself on for other datasets). Every other new `PopDetect_*` / `PopRepair_*` field absent means **use the stated default**.
- **Measured constants, do not recompute:** `c = [1.69 1.24 1.09]` for polynomial orders 1, 2, 3 in `f_3dB = c(n)*n/(2T)`.
- **`nt_detrend` takes `wsize` in samples**, not seconds, and forces it even (`nt_detrend.m:54`).
- **Segment lists are always re-derived** from the current `EEG` via `RELAX_get_segments`, never cached across stages.
- **Commit style:** end every commit message with the two trailer lines used in this repo (`Co-Authored-By:` and `Claude-Session:`), copied from `git log -1 --format=%B` of any recent commit.

---

### Task 1: Vendor NoiseTools

Vendors the 9-file runtime closure of `nt_detrend` so the pipeline does not depend on `~/TMS_response_variability/NoiseTools` being present, and puts it on the path from the SLURM entry point.

**Files:**
- Create: `vendor/noisetools/nt_detrend.m`, `nt_regw.m`, `nt_demean.m`, `nt_demean2.m`, `nt_fold.m`, `nt_unfold.m`, `nt_vecadd.m`, `nt_vecmult.m`, `nt_greetings.m` (copied verbatim)
- Create: `vendor/noisetools/VENDORING.md`
- Create: `RELAX_add_vendored_noisetools.m`
- Create: `tests/test_RELAX_add_vendored_noisetools.m`
- Modify: `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` (add the path call next to the existing `RELAX_add_vendored_zapline()` call)

**Interfaces:**
- Consumes: nothing.
- Produces: `RELAX_add_vendored_noisetools()` - no arguments, no return value, throws on failure. After it runs, `nt_detrend` resolves to `vendor/noisetools/nt_detrend.m`. Task 4 depends on this.

- [ ] **Step 1: Copy the nine files**

```bash
cd /home/nga4004/gitRepos/RELAX-GEDAI
mkdir -p vendor/noisetools
SRC=~/TMS_response_variability/NoiseTools
for f in nt_detrend nt_regw nt_demean nt_demean2 nt_fold nt_unfold nt_vecadd nt_vecmult nt_greetings; do
  cp "$SRC/$f.m" vendor/noisetools/
done
ls vendor/noisetools/
```

Expected: exactly nine `.m` files.

Do **not** copy `nt_detrend (original).m` or `nt_regw (original).m` - those are older upstream backups in the source tree, not local patches.

- [ ] **Step 2: Verify the closure is actually complete**

```bash
cd /home/nga4004/gitRepos/RELAX-GEDAI/vendor/noisetools
grep -ohE '\bnt_[a-z0-9_]+' *.m | sort -u > /tmp/called.txt
ls *.m | sed 's/\.m$//' | sort > /tmp/present.txt
comm -23 /tmp/called.txt /tmp/present.txt
```

Expected output: exactly these three names and nothing else -
```
nt_detrend2
nt_imagescc
nt_read_data
```
These three appear only inside the `if 0` demo block at the bottom of `nt_detrend.m` (from line 273 onward) and are never reached at runtime. `nt_detrend_helper` must NOT appear - it is a local subfunction at `nt_detrend.m:163`. If any other name appears, that file is a genuine missing dependency: copy it from `$SRC` and re-run this step.

- [ ] **Step 3: Write VENDORING.md**

```bash
cd /home/nga4004/gitRepos/RELAX-GEDAI
{
cat <<'HDR'
# Vendored NoiseTools (partial)

Source: `~/TMS_response_variability/NoiseTools` (the tree used by the TMS_response_variability
repo), copied 2026-09-14. Upstream project: http://audition.ens.fr/adc/NoiseTools
(public mirror: https://github.com/deron2/NoiseTools).

## This is a PARTIAL mirror, unlike vendor/zapline-plus

`vendor/zapline-plus` is a complete byte-for-byte mirror because RELAX shadows an
installed EEGLAB plugin with it, and a partial shadow would mix two versions. Nothing
shadows NoiseTools - it is not an EEGLAB plugin and is not otherwise installed - so
only the runtime closure of `nt_detrend` is vendored: 9 files instead of 138, and
6.2 MB of demo data and unrelated toolboxes left behind.

## Closure

`RELAX_repair_steps.m` calls exactly one NoiseTools entry point, `nt_detrend`. Its
runtime closure is:

    nt_detrend -> nt_regw, nt_demean
    nt_regw    -> nt_demean, nt_fold, nt_unfold
    nt_demean  -> nt_demean2, nt_fold, nt_unfold, nt_vecadd, nt_vecmult
    nt_fold    -> nt_greetings, nt_unfold
    nt_unfold  -> nt_fold, nt_greetings
    nt_vecmult -> nt_fold, nt_unfold

Three names referenced from `nt_detrend.m` are deliberately absent:

  - `nt_detrend_helper` is a LOCAL SUBFUNCTION at `nt_detrend.m:163`, not a file.
  - `nt_detrend2`, `nt_read_data` and `nt_imagescc` appear only inside the `if 0`
    demo block from `nt_detrend.m:273` onward and are unreachable at runtime.

To verify the closure after any change, run from this directory:

    grep -ohE '\bnt_[a-z0-9_]+' *.m | sort -u > /tmp/called.txt
    ls *.m | sed 's/\.m$//' | sort > /tmp/present.txt
    comm -23 /tmp/called.txt /tmp/present.txt

That must print exactly `nt_detrend2`, `nt_imagescc`, `nt_read_data` and nothing else.

## Why the local tree rather than the public mirror

The public mirror's `nt_regw.m` has its all-zero-weight guard commented out, so a
degenerate fit silently returns `b = 0` and propagates `NaN` through the detrended
output. The local tree raises `error('weights all zero')` instead. A loud failure is
required here: `RELAX_repair_steps` catches it per sub-segment and steps down its
treatment ladder, which is only possible if the failure is signalled at all.

## Customizations

NONE. All nine files are byte-identical copies of the source tree. `RELAX_repair_steps`
works around `nt_detrend`'s rough edges from the outside (try/catch around each fit,
never handing it an all-zero-weight window) rather than patching it, so a later re-sync
is a straight overwrite.

## Re-sync procedure

1. Copy the nine files again from the source tree (or from upstream).
2. Re-run the closure check above.
3. Run `tests/test_RELAX_repair_steps.m` - in particular the spectral-fidelity case,
   which fails if `nt_detrend`'s transfer function changed.
4. Update the checksum table below.

## Checksums at time of vendoring

HDR
printf '| file | lines | sha256 |\n| --- | --- | --- |\n'
for f in vendor/noisetools/*.m; do
  printf '| `%s` | %s | `%s` |\n' "$(basename $f)" "$(wc -l < $f)" "$(sha256sum $f | cut -c1-16)"
done
} > vendor/noisetools/VENDORING.md
head -5 vendor/noisetools/VENDORING.md; tail -12 vendor/noisetools/VENDORING.md
```

- [ ] **Step 4: Write the failing test**

Create `tests/test_RELAX_add_vendored_noisetools.m`:

```matlab
function test_RELAX_add_vendored_noisetools()
% The vendored NoiseTools closure is present, self-contained, and resolves on the path.

    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    vendorDir = fullfile(repoRoot, 'vendor', 'noisetools');

    assert(isfolder(vendorDir), 'vendor/noisetools must exist');

    expected = {'nt_detrend','nt_regw','nt_demean','nt_demean2','nt_fold', ...
                'nt_unfold','nt_vecadd','nt_vecmult','nt_greetings'};
    for k = 1:numel(expected)
        assert(isfile(fullfile(vendorDir, [expected{k} '.m'])), ...
            'vendored file %s.m is missing', expected{k});
    end

    d = dir(fullfile(vendorDir, '*.m'));
    assert(numel(d) == numel(expected), ...
        'vendor/noisetools should hold exactly %d .m files, found %d', numel(expected), numel(d));

    % The helper puts nt_detrend on the path and points at the vendored copy
    RELAX_add_vendored_noisetools();
    resolved = which('nt_detrend');
    assert(strncmp(resolved, vendorDir, numel(vendorDir)), ...
        'nt_detrend resolves to %s, not the vendored copy in %s', resolved, vendorDir);

    % And it actually runs: a pure linear ramp is removed by an order-1 fit
    x = (1:1000)' * 0.5;
    y = nt_detrend(x, 1);
    assert(max(abs(y)) < 1e-6, 'order-1 detrend should flatten a linear ramp, max residual %g', max(abs(y)));

    fprintf('test_RELAX_add_vendored_noisetools: PASS\n');
end
```

- [ ] **Step 5: Run the test to verify it fails**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_add_vendored_noisetools"
```

Expected: FAIL with `Unrecognized function or variable 'RELAX_add_vendored_noisetools'`.

- [ ] **Step 6: Write the path helper**

Create `RELAX_add_vendored_noisetools.m`:

```matlab
function RELAX_add_vendored_noisetools()
% RELAX_ADD_VENDORED_NOISETOOLS - put the vendored NoiseTools closure on the path.
%
% RELAX_repair_steps calls nt_detrend. Rather than depend on a NoiseTools checkout
% living next to this repo, vendor/noisetools holds the 9-file runtime closure of
% nt_detrend (see vendor/noisetools/VENDORING.md). Unlike RELAX_add_vendored_zapline
% this shadows nothing - NoiseTools is not an EEGLAB plugin - so the call order
% relative to eeglab() does not matter.

    vendorDir = fullfile(fileparts(mfilename('fullpath')), 'vendor', 'noisetools');
    if ~isfolder(vendorDir)
        error('RELAX_add_vendored_noisetools:missing', ...
            'Vendored NoiseTools not found at %s', vendorDir);
    end

    addpath(vendorDir);

    resolved = which('nt_detrend');
    if ~strncmp(resolved, vendorDir, numel(vendorDir))
        error('RELAX_add_vendored_noisetools:shadowed', ...
            ['nt_detrend resolves to\n  %s\ninstead of the vendored copy in\n  %s\n' ...
             'Another NoiseTools copy is earlier on the path - remove it, or move this ' ...
             'call later so its addpath wins.'], resolved, vendorDir);
    end
end
```

- [ ] **Step 7: Run the test to verify it passes**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_add_vendored_noisetools"
```

Expected: `test_RELAX_add_vendored_noisetools: PASS`

- [ ] **Step 8: Wire it into the SLURM entry point**

Find the existing `RELAX_add_vendored_zapline()` call in `RELAX_SET_PARAMETERS_AND_RUN_slurm.m`:

```bash
grep -n 'RELAX_add_vendored_zapline' RELAX_SET_PARAMETERS_AND_RUN_slurm.m
```

Add the new call on the line immediately after it:

```matlab
RELAX_add_vendored_noisetools();  % nt_detrend for RELAX_repair_steps, see vendor/noisetools/VENDORING.md
```

- [ ] **Step 9: Run the full suite**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: all test files pass, count increased by one from the pre-task baseline.

- [ ] **Step 10: Commit**

```bash
git add vendor/noisetools RELAX_add_vendored_noisetools.m tests/test_RELAX_add_vendored_noisetools.m RELAX_SET_PARAMETERS_AND_RUN_slurm.m
git commit -m "Vendor the nt_detrend runtime closure from NoiseTools

Nine files, the full runtime closure of nt_detrend, copied from the local
NoiseTools tree rather than the public mirror (whose nt_regw has its
all-zero-weight guard commented out and returns NaN instead of raising).

Partial rather than complete mirror because nothing shadows NoiseTools -
see vendor/noisetools/VENDORING.md for the closure analysis and the
verification command.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01TRpqXm5P93dw7q3f6wpWgQ"
```

---

### Task 2: Shared channel-rejection budget

Creates a running budget so the pop, impedance and PREP screens together cannot delete more than `MaxProportionOfElectrodesThatCanBeDeleted` of the original montage, and teaches the impedance screen to spend from it.

**Files:**
- Create: `RELAX_channel_budget.m`
- Create: `tests/test_RELAX_channel_budget.m`
- Modify: `RELAX_impedance_bad_channels.m` (optional 4th argument, around line 117)

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `EEG = RELAX_channel_budget(EEG, 'init', maxProportion)` - creates `EEG.RELAX.ChannelBudget`.
  - `[fits, EEG] = RELAX_channel_budget(EEG, 'spend', stageName, worstFirstLabels)` - `worstFirstLabels` is a cell array of channel labels ordered worst-first; `fits` is the prefix of it that the budget allows (a cell array, possibly empty).
  - `RELAX_impedance_bad_channels(chanlocs, maxKOhm, maxProportion, maxRejectableCount)` - the 4-argument form overrides the internal cap; the 3-argument form is unchanged.
  - Tasks 3 and 5 depend on both.

- [ ] **Step 1: Write the failing test**

Create `tests/test_RELAX_channel_budget.m`:

```matlab
function test_RELAX_channel_budget()
% Shared running channel-rejection budget across screens.

    mkEEG = @(n) struct('nbchan', n, 'chanlocs', struct('labels', ...
        arrayfun(@(k) sprintf('E%d', k), 1:n, 'UniformOutput', false)));

    % init records the montage size and the derived cap
    EEG = RELAX_channel_budget(mkEEG(100), 'init', 0.10);
    b = EEG.RELAX.ChannelBudget;
    assert(b.nOriginalChannels == 100, 'nOriginalChannels');
    assert(b.maxRejectable == 10, 'maxRejectable should be floor(0.10*100)');
    assert(b.nSpent == 0, 'nothing spent yet');
    assert(~b.exhausted, 'not exhausted yet');

    % A request inside the budget passes through untouched
    [fits, EEG] = RELAX_channel_budget(EEG, 'spend', 'pop', {'E1','E2','E3'});
    assert(isequal(fits, {'E1','E2','E3'}), 'in-budget request should pass through whole');
    assert(EEG.RELAX.ChannelBudget.nSpent == 3, 'nSpent should be 3');
    assert(isequal(EEG.RELAX.ChannelBudget.spentBy.pop, {'E1','E2','E3'}), 'attribution by stage');

    % A second screen sees the REDUCED remainder, and truncation keeps the worst first
    [fits, EEG] = RELAX_channel_budget(EEG, 'spend', 'impedance', ...
        {'E10','E11','E12','E13','E14','E15','E16','E17','E18','E19'});
    assert(numel(fits) == 7, 'only 7 of the 10 should fit into the remaining budget');
    assert(isequal(fits, {'E10','E11','E12','E13','E14','E15','E16'}), ...
        'truncation must keep the worst-first prefix');
    assert(EEG.RELAX.ChannelBudget.nSpent == 10, 'budget fully spent');
    assert(EEG.RELAX.ChannelBudget.exhausted, 'exhausted flag must be set');

    % A third screen gets nothing, and gets an empty cell rather than []
    [fits, EEG] = RELAX_channel_budget(EEG, 'spend', 'PREP', {'E20'});
    assert(iscell(fits) && isempty(fits), 'exhausted budget must return an empty cell');
    assert(EEG.RELAX.ChannelBudget.nSpent == 10, 'nSpent must not grow past the cap');

    % Spending nothing is legal and does not set the exhausted flag
    EEG2 = RELAX_channel_budget(mkEEG(100), 'init', 0.10);
    [fits, EEG2] = RELAX_channel_budget(EEG2, 'spend', 'pop', {});
    assert(iscell(fits) && isempty(fits), 'empty request returns empty cell');
    assert(~EEG2.RELAX.ChannelBudget.exhausted, 'empty request must not exhaust the budget');

    % A budget that rounds down to zero rejects everything requested
    EEG3 = RELAX_channel_budget(mkEEG(5), 'init', 0.10);
    assert(EEG3.RELAX.ChannelBudget.maxRejectable == 0, 'floor(0.10*5) is 0');
    [fits, ~] = RELAX_channel_budget(EEG3, 'spend', 'pop', {'E1'});
    assert(isempty(fits), 'a zero budget grants nothing');

    % Init must not disturb an existing EEG.RELAX
    EEG4 = mkEEG(50);
    EEG4.RELAX.SomethingElse = 42;
    EEG4 = RELAX_channel_budget(EEG4, 'init', 0.10);
    assert(EEG4.RELAX.SomethingElse == 42, 'init must not clobber other EEG.RELAX fields');

    % --- the impedance screen's backward compatibility ---
    % RELAX_impedance_bad_channels reads kOhm from chanlocs(k).impedance, and returns
    % a COLUMN cell of labels (see tests/test_RELAX_impedance_bad_channels.m).
    chanlocs = struct('labels', {'E1','E2','E3','E4','E5'}, ...
                      'impedance', {900, 800, 700, 10, 10});

    % 3-argument form: cap is floor(0.10*5) = 0, so nothing is rejected
    rej3 = RELAX_impedance_bad_channels(chanlocs, 500, 0.10);
    assert(isempty(rej3), '3-arg form must keep its existing floor(maxProportion*n) cap');

    % 4-argument form: the explicit count overrides, worst impedance first
    rej4 = RELAX_impedance_bad_channels(chanlocs, 500, 0.10, 2);
    assert(numel(rej4) == 2, '4-arg form should honour the explicit count');
    assert(all(ismember({'E1','E2'}, rej4)), 'the two worst impedances must be the ones kept');

    % An explicit count large enough to not bind returns everything over threshold
    rej5 = RELAX_impedance_bad_channels(chanlocs, 500, 0.10, 99);
    assert(numel(rej5) == 3, 'a non-binding count should return all 3 over-threshold channels');

    fprintf('test_RELAX_channel_budget: PASS\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_channel_budget"
```

Expected: FAIL with `Unrecognized function or variable 'RELAX_channel_budget'`.

- [ ] **Step 3: Write the budget helper**

Create `RELAX_channel_budget.m`:

```matlab
% RELAX_channel_budget() - Running budget shared by the channel-rejection screens.
%
% Usage:
%   >> EEG          = RELAX_channel_budget(EEG, 'init', maxProportion);
%   >> [fits, EEG]  = RELAX_channel_budget(EEG, 'spend', stageName, worstFirstLabels);
%
% Inputs:
%   EEG              - EEGLAB EEG structure. 'init' reads .nbchan; 'spend' reads and
%                       writes EEG.RELAX.ChannelBudget.
%   maxProportion    - proportion of the ORIGINAL montage that all participating screens
%                       may delete between them. maxRejectable = floor(maxProportion*nbchan).
%   stageName        - char, a valid MATLAB field name ('pop', 'impedance', 'PREP'), used
%                       to attribute the spend in .spentBy.
%   worstFirstLabels - cell array of channel labels, ordered WORST FIRST. When the budget
%                       cannot cover the whole list it is truncated from the end, so the
%                       channels sacrificed are the least-bad ones.
%
% Outputs:
%   fits - the prefix of worstFirstLabels the budget allows. Always a cell array, {} when
%          nothing fits, so callers can pass it straight to pop_select without a type test.
%   EEG  - EEG with .RELAX.ChannelBudget updated.
%
% Call 'init' immediately after EEG.allchan = EEG.chanlocs in the wrapper: that line
% records the reference montage the proportion is taken against, so the budget has to be
% created against the same snapshot, before any screen has removed a channel.
%
% Screens that do NOT participate, deliberately:
%   - the bridged-channel screen. A bridge is a physical short; a retained cluster member
%     holds a smeared mixture while its chanlocs entry claims a precise location, so
%     capping it would preserve channels known to be invalid.
%   - RELAX_excluding_channels_and_epoching, which keeps its own independent later cap.
%     Bringing it in means touching three sibling variants and is a separate change.

function [varargout] = RELAX_channel_budget(EEG, op, varargin)

    switch lower(op)
        case 'init'
            maxProportion = varargin{1};
            budget = struct( ...
                'nOriginalChannels', EEG.nbchan, ...
                'maxProportion',     maxProportion, ...
                'maxRejectable',     floor(maxProportion * EEG.nbchan), ...
                'nSpent',            0, ...
                'spentBy',           struct(), ...
                'exhausted',         false);
            EEG.RELAX.ChannelBudget = budget;  % leaves any other EEG.RELAX field alone
            varargout = {EEG};

        case 'spend'
            stageName        = varargin{1};
            worstFirstLabels = varargin{2};

            if ~isfield(EEG, 'RELAX') || ~isfield(EEG.RELAX, 'ChannelBudget')
                error('RELAX_channel_budget:notInitialised', ...
                    ['EEG.RELAX.ChannelBudget does not exist. Call ' ...
                     'RELAX_channel_budget(EEG, ''init'', maxProportion) before any screen spends.']);
            end

            if isempty(worstFirstLabels)
                worstFirstLabels = {};
            end
            worstFirstLabels = worstFirstLabels(:)';  % row, so the prefix slice is a row cell

            b = EEG.RELAX.ChannelBudget;
            remaining = max(0, b.maxRejectable - b.nSpent);
            nFits = min(numel(worstFirstLabels), remaining);
            fits = worstFirstLabels(1:nFits);

            b.nSpent = b.nSpent + nFits;
            b.spentBy.(stageName) = fits;
            if nFits < numel(worstFirstLabels)
                b.exhausted = true;
                warning('RELAX_channel_budget:exhausted', ...
                    ['Stage ''%s'' asked to reject %d channel(s) but only %d fit in the ' ...
                     'shared budget (%d of %d channels, %g of the montage). ' ...
                     'Rejecting the worst %d; the rest are recorded but retained.'], ...
                    stageName, numel(worstFirstLabels), nFits, ...
                    b.maxRejectable, b.nOriginalChannels, b.maxProportion, nFits);
            end
            EEG.RELAX.ChannelBudget = b;

            varargout = {fits, EEG};

        otherwise
            error('RELAX_channel_budget:unknownOp', ...
                'Unknown operation ''%s''. Expected ''init'' or ''spend''.', op);
    end
end
```

- [ ] **Step 4: Add the optional 4th argument to the impedance screen**

In `RELAX_impedance_bad_channels.m`, change the signature at line 50:

```matlab
function [rejectLabels, info] = RELAX_impedance_bad_channels(chanlocs, maxKOhm, maxProportion, maxRejectableCount)
```

Immediately after the signature line, add:

```matlab
    % Optional 4th argument: an explicit cap in channels, from the shared budget in
    % RELAX_channel_budget. Overrides the internal floor(maxProportion*nChannels) so
    % this screen spends only what earlier screens left. The 3-argument form is
    % unchanged, which is what keeps the sibling wrappers byte-identical.
    if nargin < 4
        maxRejectableCount = [];
    end
```

Then at line 117, replace:

```matlab
    maxRejectable = floor(maxProportion * nChannels);
```

with:

```matlab
    if isempty(maxRejectableCount)
        maxRejectable = floor(maxProportion * nChannels);
    else
        maxRejectable = maxRejectableCount;
    end
```

Add `maxRejectableCount` to the header comment block's input list, above the `function` line:

```matlab
%   maxRejectableCount - (optional) explicit cap in channels, overriding
%                        floor(maxProportion*nChannels). Supplied by RELAX_channel_budget
%                        so this screen spends only the shared budget's remainder. Omit
%                        for the original single-screen behaviour.
```

- [ ] **Step 5: Run the test to verify it passes**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_channel_budget"
```

Expected: `test_RELAX_channel_budget: PASS`

- [ ] **Step 6: Run the existing impedance test to prove nothing regressed**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_impedance_bad_channels"
```

Expected: PASS, unchanged. This test only ever calls the 3-argument form, so it is the backward-compatibility check.

- [ ] **Step 7: Commit**

```bash
git add RELAX_channel_budget.m tests/test_RELAX_channel_budget.m RELAX_impedance_bad_channels.m
git commit -m "Add a shared channel-rejection budget across screens

RELAX_channel_budget tracks a running count so the pop, impedance and PREP
screens together cannot delete more than MaxProportionOfElectrodesThatCanBeDeleted
of the original montage. Previously each screen independently capped itself at
that proportion, so three screens could remove three times it.

Truncation keeps the worst-first prefix, so an exhausted budget sacrifices the
least-bad candidates. RELAX_impedance_bad_channels gains an optional 4th argument
for the remaining count; its 3-argument form is unchanged.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01TRpqXm5P93dw7q3f6wpWgQ"
```

---

### Task 3: `RELAX_detect_steps`

Detects electrode pops with a median-filtered derivative, a robust z threshold, an absolute amplitude floor, a sustain check and a rise-time criterion.

**Files:**
- Create: `RELAX_detect_steps.m`
- Create: `tests/test_RELAX_detect_steps.m`

**Interfaces:**
- Consumes: `RELAX_get_segments(EEG)` from the boundary work - returns `[nSeg x 2]` inclusive 1-based sample ranges.
- Produces: `[steps, info] = RELAX_detect_steps(EEG, 'Name', value, ...)`.
  - `steps` is a struct with parallel column fields `.channel` (double), `.label` (cell), `.sample` (double), `.amplitude` (double, signed uV), `.z` (double), each `[N x 1]`, sorted by channel then sample. On no detection every field is `0x1` / `{}`, never bare `[]`.
  - `info` is a struct with `.available` (logical), `.status` (char), `.params` (struct), `.nSteps`, `.nChannelsAffected`, `.stepsPerChannel` `[nbchan x 1]`, `.stepRatePerMinute` `[nbchan x 1]`, `.badChannels` (cell), `.badChannelReason` (cell), `.degenerateChannels` (cell), `.nCommonModeSteps`, `.commonModeSamples`.
  - Task 4 consumes `steps`; Task 5 consumes both.

**Note for the implementer:** the spec's config block omits two parameters its own detection step 5 requires - `SustainMs` and `MinRiseFraction`. They are defined here and added to the spec's config table in Task 5 Step 8.

- [ ] **Step 1: Write the failing test**

Create `tests/test_RELAX_detect_steps.m`:

```matlab
function test_RELAX_detect_steps()
% Step detection: recovers real steps, rejects ramps and noise, degrades safely.

    srate = 1000;
    n     = 10000;          % 10 s
    rng(7);

    % ---- helper: build a 2-channel EEG with no boundaries ----
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

    % ---- 1. a known step is found at the right sample with the right amplitude ----
    x = randn(1, n) * 5;                 % 5 uV rms background
    x(4001:end) = x(4001:end) + 500;     % +500 uV step, new level begins at 4001
    EEG = mkEEG([x; randn(1, n) * 5]);

    [steps, info] = RELAX_detect_steps(EEG);

    assert(info.available, 'detector should be available on ordinary data: %s', info.status);
    onCh1 = steps.sample(steps.channel == 1);
    assert(numel(onCh1) == 1, 'expected exactly 1 step on channel 1, got %d', numel(onCh1));
    assert(abs(onCh1 - 4001) <= 2, 'step sample should be ~4001, got %d', onCh1);
    amp = steps.amplitude(steps.channel == 1);
    assert(abs(amp - 500) < 25, 'step amplitude should be ~+500 uV, got %g', amp);
    assert(steps.z(steps.channel == 1) > 20, 'z should clear the threshold');
    assert(isempty(steps.sample(steps.channel == 2)), 'clean channel 2 must yield no steps');
    assert(strcmp(steps.label{steps.channel == 1}, 'E1'), 'label should accompany the index');

    % sign is preserved for a downward step
    x2 = randn(1, n) * 5;
    x2(4001:end) = x2(4001:end) - 500;
    [steps2, ~] = RELAX_detect_steps(mkEEG(x2));
    assert(steps2.amplitude(1) < 0, 'a downward step must have negative amplitude');

    % ---- 2. NO detection on a smooth ramp (the exploratory detector's failure mode) ----
    ramp = randn(1, n) * 5 + linspace(0, 500, n);   % same 500 uV excursion, spread out
    [stepsRamp, ~] = RELAX_detect_steps(mkEEG(ramp));
    assert(isempty(stepsRamp.sample), ...
        'a smooth 500 uV ramp must not be detected as a step, got %d detections', ...
        numel(stepsRamp.sample));

    % ---- 3. no detection on noise alone ----
    [stepsNoise, ~] = RELAX_detect_steps(mkEEG(randn(2, n) * 20));
    assert(isempty(stepsNoise.sample), ...
        'plain noise must not trigger detections, got %d', numel(stepsNoise.sample));

    % ---- 4. refractory: two steps 5 ms apart collapse to one ----
    x3 = randn(1, n) * 5;
    x3(4001:end) = x3(4001:end) + 400;
    x3(4006:end) = x3(4006:end) + 400;
    [steps3, ~] = RELAX_detect_steps(mkEEG(x3), 'RefractoryMs', 20);
    assert(numel(steps3.sample) == 1, ...
        'two steps inside the refractory window must give 1 detection, got %d', numel(steps3.sample));

    % ...but 500 ms apart they stay separate
    x4 = randn(1, n) * 5;
    x4(3001:end) = x4(3001:end) + 400;
    x4(3501:end) = x4(3501:end) - 400;
    [steps4, ~] = RELAX_detect_steps(mkEEG(x4), 'RefractoryMs', 20);
    assert(numel(steps4.sample) == 2, ...
        'two well-separated steps must give 2 detections, got %d', numel(steps4.sample));

    % ---- 5. flat channel: no steps, no NaN, listed as degenerate ----
    flat = zeros(1, n);
    [steps5, info5] = RELAX_detect_steps(mkEEG([flat; randn(1,n)*5]));
    assert(isempty(steps5.sample(steps5.channel == 1)), 'a flat channel must yield no steps');
    assert(ismember('E1', info5.degenerateChannels), 'flat channel must be listed as degenerate');
    assert(~any(isnan(steps5.z)), 'no NaN z-scores');

    % ---- 6. absolute amplitude floor suppresses a tiny step on a very quiet channel ----
    quiet = randn(1, n) * 0.01;
    quiet(4001:end) = quiet(4001:end) + 2;    % 2 uV: huge in z, below a 10 uV floor
    [steps6, ~] = RELAX_detect_steps(mkEEG(quiet), 'MinAmplitudeMicrovolts', 10);
    assert(isempty(steps6.sample), 'a 2 uV step must be suppressed by a 10 uV floor');
    [steps6b, ~] = RELAX_detect_steps(mkEEG(quiet), 'MinAmplitudeMicrovolts', 1);
    assert(numel(steps6b.sample) == 1, 'the same step must survive a 1 uV floor');

    % ---- 7. boundaries: no detection manufactured at the discontinuity itself ----
    xb = randn(1, n) * 5;
    xb(5001:end) = xb(5001:end) + 800;     % a big jump that IS the boundary, not a pop
    EEGb = mkEEG(xb);
    EEGb.event = struct('type', {'boundary'}, 'latency', {5000.5});
    [stepsB, ~] = RELAX_detect_steps(EEGb);
    assert(isempty(stepsB.sample), ...
        'a level change coinciding with a boundary is a splice, not a pop; got %d detections', ...
        numel(stepsB.sample));

    % a genuine pop well inside the second segment is still found
    xb2 = randn(1, n) * 5;
    xb2(5001:end) = xb2(5001:end) + 800;
    xb2(7001:end) = xb2(7001:end) + 400;
    EEGb2 = mkEEG(xb2);
    EEGb2.event = struct('type', {'boundary'}, 'latency', {5000.5});
    [stepsB2, ~] = RELAX_detect_steps(EEGb2);
    assert(numel(stepsB2.sample) == 1 && abs(stepsB2.sample - 7001) <= 2, ...
        'the pop inside the segment must still be found');

    % ---- 8. Population A: a chronically popping channel is flagged for rejection ----
    xa = randn(1, n) * 5;
    for s = 500:500:9500                    % 19 steps in 10 s = 114/min
        xa(s:end) = xa(s:end) + 300 * (-1)^(s/500);
    end
    [~, infoA] = RELAX_detect_steps(mkEEG([xa; randn(1,n)*5]), ...
        'StepsPerMinuteToRejectChannel', 30);
    assert(ismember('E1', infoA.badChannels), 'chronically popping channel must be flagged');
    assert(~ismember('E2', infoA.badChannels), 'the clean channel must not be flagged');
    assert(infoA.stepRatePerMinute(1) > 30, 'rate should exceed the threshold');
    assert(numel(infoA.badChannelReason) == numel(infoA.badChannels), 'one reason per bad channel');

    % ---- 9. empty result has the right shape, not a bare [] ----
    [stepsE, infoE] = RELAX_detect_steps(mkEEG(zeros(2, n) + randn(2, n)*1e-9));
    assert(isnumeric(stepsE.sample) && isempty(stepsE.sample) && size(stepsE.sample,2) <= 1, ...
        'empty .sample must be a 0x1 numeric');
    assert(iscell(stepsE.label), 'empty .label must still be a cell');
    assert(infoE.nSteps == 0 && infoE.nChannelsAffected == 0, 'zero counts on an empty result');

    % ---- 10. params are echoed so the saved struct is self-describing ----
    [~, infoP] = RELAX_detect_steps(mkEEG(randn(2,n)*5), 'ZThreshold', 33);
    assert(infoP.params.ZThreshold == 33, 'params must echo the effective settings');

    fprintf('test_RELAX_detect_steps: PASS\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_detect_steps"
```

Expected: FAIL with `Unrecognized function or variable 'RELAX_detect_steps'`.

- [ ] **Step 3: Write the detector**

Create `RELAX_detect_steps.m`:

```matlab
% RELAX_detect_steps() - Detect electrode-pop step artifacts in continuous EEG.
%
% Usage:
%   >> [steps, info] = RELAX_detect_steps(EEG);
%   >> [steps, info] = RELAX_detect_steps(EEG, 'ZThreshold', 20, 'RefractoryMs', 20);
%
% An electrode pop is a near-instantaneous change in a single channel's DC baseline,
% caused by a momentary change in the electrode-electrolyte contact. It is a STEP, not a
% transient: the baseline moves and stays moved. Detection follows FieldTrip's
% ft_artifact_jump recipe (median filter, absolute derivative, z threshold), with two
% additions this data needed: an absolute microvolt floor, and a rise-time criterion
% without which smooth ramps are flagged as steps.
%
% Run this BEFORE any filtering. A zero-phase highpass converts a step into a bipolar
% excursion spread over seconds in BOTH time directions, after which no step remains to
% find.
%
% Inputs:
%   EEG - EEGLAB EEG structure (uses .data, .nbchan, .pnts, .srate, .chanlocs, .event)
%   Name-value pairs (all optional):
%     MedianFilterOrder    - median filter length in samples, suppresses single-sample
%                             spikes without smearing a genuine step edge (default 9,
%                             FieldTrip's value)
%     ZThreshold           - robust z threshold on the absolute derivative (default 20,
%                             FieldTrip's value)
%     MinAmplitudeMicrovolts - absolute floor on the level change (default 10). Stops a
%                             very quiet channel generating detections from its own noise
%                             floor, where a tiny step is enormous in z.
%     RefractoryMs         - minimum spacing between steps on one channel (default 20)
%     SustainMs            - window either side used to measure the pre- and post-step
%                             level (default 100). The level change must persist across
%                             it - that is what makes this a step detector rather than a
%                             spike detector.
%     MinRiseFraction      - fraction of the total level change that must occur in the
%                             single largest derivative sample (default 0.5). This is the
%                             rise-time criterion: a genuine pop delivers most of its
%                             change in one sample, a ramp spreads it over thousands.
%     StepsPerMinuteToRejectChannel - per-channel step rate above which the channel is a
%                             chronically bad contact rather than a repairable one
%                             (default 30). Derived from survey data: the clearest bad
%                             contact ran at ~72 steps/min while sporadically popping
%                             channels averaged well under 1/min, so 30 sits in the gap.
%     CommonModeBinMs      - bin width for the reference-pop co-occurrence count (default 25)
%     CommonModeChannelFraction - fraction of channels in one bin above which the bin is
%                             counted as a common-mode (reference) step (default 0.25).
%                             REPORTING ONLY - such steps are still repaired per channel.
%
% Outputs:
%   steps - struct of parallel [N x 1] columns, one entry per step, sorted by channel then
%           sample: .channel (row index into EEG.data), .label (channel label, so the
%           struct survives a pop_select), .sample (FIRST sample of the new level),
%           .amplitude (signed uV, post-level minus pre-level), .z. An empty result has
%           0x1 / {} fields, never a bare [], so numel(steps.sample) is always valid.
%   info  - .available, .status, .params, .nSteps, .nChannelsAffected, .stepsPerChannel,
%           .stepRatePerMinute, .badChannels, .badChannelReason, .degenerateChannels,
%           .nCommonModeSteps, .commonModeSamples
%
% Detection runs inside each inter-boundary segment separately (via RELAX_get_segments),
% so the splice discontinuity at a boundary is never mistaken for a pop. The robust
% z-scale, however, is computed per channel across the WHOLE recording: an 8 s segment is
% too short for a stable scale estimate, and MAD is insensitive to the pops themselves.

function [steps, info] = RELAX_detect_steps(EEG, varargin)

    p = inputParser;
    p.addParameter('MedianFilterOrder', 9);
    p.addParameter('ZThreshold', 20);
    p.addParameter('MinAmplitudeMicrovolts', 10);
    p.addParameter('RefractoryMs', 20);
    p.addParameter('SustainMs', 100);
    p.addParameter('MinRiseFraction', 0.5);
    p.addParameter('StepsPerMinuteToRejectChannel', 30);
    p.addParameter('CommonModeBinMs', 25);
    p.addParameter('CommonModeChannelFraction', 0.25);
    p.parse(varargin{:});
    opt = p.Results;

    steps = emptySteps();
    info  = struct('available', false, 'status', '', 'params', opt, ...
                   'nSteps', 0, 'nChannelsAffected', 0, ...
                   'stepsPerChannel', zeros(EEG.nbchan,1), ...
                   'stepRatePerMinute', zeros(EEG.nbchan,1), ...
                   'badChannels', {{}}, 'badChannelReason', {{}}, ...
                   'degenerateChannels', {{}}, ...
                   'nCommonModeSteps', 0, 'commonModeSamples', []);

    refractory = max(1, round(opt.RefractoryMs * EEG.srate / 1000));
    sustain    = max(2, round(opt.SustainMs    * EEG.srate / 1000));

    segments = RELAX_get_segments(EEG);
    % A segment must hold a pre-window, a post-window and the step itself to be usable.
    minSegLen = 2*sustain + 2;
    usable = (segments(:,2) - segments(:,1) + 1) >= minSegLen;
    if ~any(usable)
        info.status = sprintf(['no inter-boundary segment reaches the %d samples needed for a ' ...
            '%g ms sustain window either side'], minSegLen, opt.SustainMs);
        return
    end
    segments = segments(usable, :);

    labels = {EEG.chanlocs.labels};
    durationMinutes = (EEG.pnts / EEG.srate) / 60;

    allCh = []; allSamp = []; allAmp = []; allZ = [];
    degenerate = {};

    for ch = 1:EEG.nbchan
        x = double(EEG.data(ch, :));

        % 1. median-filter and differentiate, inside each segment only
        d = nan(1, EEG.pnts - 1);   % d(i) is the transition between sample i and i+1
        for s = 1:size(segments,1)
            a = segments(s,1); b = segments(s,2);
            xm = movmedian(x(a:b), opt.MedianFilterOrder);
            d(a : b-1) = abs(diff(xm));
        end
        dValid = d(~isnan(d));

        % 2. robust scale across the whole recording
        med = median(dValid);
        mad = median(abs(dValid - med));
        if mad <= 0
            % A flat, dead or heavily quantised channel. No scale exists, so no z exists.
            degenerate{end+1} = labels{ch}; %#ok<AGROW>
            continue
        end
        z = (d - med) / (1.4826 * mad);

        % 3. candidates over threshold, strongest first, refractory enforced greedily
        cand = find(z > opt.ZThreshold);
        if isempty(cand); continue; end
        [~, order] = sort(z(cand), 'descend');
        cand = cand(order);

        taken = false(1, EEG.pnts);
        for k = 1:numel(cand)
            i = cand(k);                       % transition index
            lo = max(1, i - refractory);
            hi = min(EEG.pnts, i + refractory);
            if any(taken(lo:hi)); continue; end

            stepSample = i + 1;                % first sample of the NEW level

            % the sustain windows must fit inside THIS segment
            s = find(segments(:,1) <= i & segments(:,2) >= stepSample, 1);
            if isempty(s); continue; end
            a = segments(s,1); b = segments(s,2);
            if (i - a + 1) < sustain || (b - stepSample + 1) < sustain; continue; end

            preLevel  = median(x(i - sustain + 1 : i));
            postLevel = median(x(stepSample : stepSample + sustain - 1));
            amplitude = postLevel - preLevel;

            % 4. absolute floor
            if abs(amplitude) < opt.MinAmplitudeMicrovolts; continue; end

            % 5. rise-time criterion: most of the change in one sample.
            %    A ramp delivers its change over thousands of samples and fails here.
            if abs(x(stepSample) - x(i)) < opt.MinRiseFraction * abs(amplitude); continue; end

            taken(lo:hi) = true;
            allCh(end+1,1)   = ch;          %#ok<AGROW>
            allSamp(end+1,1) = stepSample;  %#ok<AGROW>
            allAmp(end+1,1)  = amplitude;   %#ok<AGROW>
            allZ(end+1,1)    = z(i);        %#ok<AGROW>
        end
    end

    info.degenerateChannels = degenerate;
    info.available = true;
    info.status    = 'ok';

    if isempty(allCh)
        return
    end

    % sort by channel then sample
    [~, ord] = sortrows([allCh allSamp]);
    steps.channel   = allCh(ord);
    steps.sample    = allSamp(ord);
    steps.amplitude = allAmp(ord);
    steps.z         = allZ(ord);
    steps.label     = reshape(labels(steps.channel), [], 1);

    info.nSteps            = numel(steps.sample);
    info.stepsPerChannel   = accumarray(steps.channel, 1, [EEG.nbchan 1]);
    info.stepRatePerMinute = info.stepsPerChannel / max(durationMinutes, eps);
    info.nChannelsAffected = sum(info.stepsPerChannel > 0);

    % Population A: chronically bad contacts, ordered worst-first for the shared budget
    badIdx = find(info.stepRatePerMinute > opt.StepsPerMinuteToRejectChannel);
    [~, worstFirst] = sort(info.stepRatePerMinute(badIdx), 'descend');
    badIdx = badIdx(worstFirst);
    info.badChannels = reshape(labels(badIdx), [], 1);
    info.badChannelReason = reshape(arrayfun(@(k) sprintf('%.1f steps/min exceeds %g', ...
        info.stepRatePerMinute(k), opt.StepsPerMinuteToRejectChannel), ...
        badIdx, 'UniformOutput', false), [], 1);

    % Common-mode (reference) steps: reporting only, per design decision 10
    binSamples = max(1, round(opt.CommonModeBinMs * EEG.srate / 1000));
    binIdx = floor((steps.sample - 1) / binSamples) + 1;
    uniqueBins = unique(binIdx);
    nPerBin = arrayfun(@(b) numel(unique(steps.channel(binIdx == b))), uniqueBins);
    commonBins = uniqueBins(nPerBin > opt.CommonModeChannelFraction * EEG.nbchan);
    info.nCommonModeSteps  = numel(commonBins);
    info.commonModeSamples = (commonBins - 1) * binSamples + 1;
end

function s = emptySteps()
    s = struct('channel', zeros(0,1), 'label', {cell(0,1)}, 'sample', zeros(0,1), ...
               'amplitude', zeros(0,1), 'z', zeros(0,1));
end
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_detect_steps"
```

Expected: `test_RELAX_detect_steps: PASS`

If case 2 (the ramp) fails, `MinRiseFraction` is doing nothing - check that `abs(x(stepSample) - x(i))` uses the **unfiltered** `x`, not the median-filtered trace. If case 7 (boundaries) fails, the derivative is being taken across the boundary - check that `d` is left `NaN` at index `b` of each segment.

- [ ] **Step 5: Commit**

```bash
git add RELAX_detect_steps.m tests/test_RELAX_detect_steps.m
git commit -m "Add RELAX_detect_steps: electrode-pop detection

Median-filtered derivative plus robust z (FieldTrip ft_artifact_jump recipe),
with two additions this dataset needed: an absolute microvolt floor so a very
quiet channel cannot generate detections from its own noise floor, and a
rise-time criterion without which smooth ramps are flagged as steps.

Detection runs inside each inter-boundary segment so a splice is never taken
for a pop, while the robust z-scale is computed per channel across the whole
recording, where an 8 s segment would be too short to estimate it stably.

Channels above StepsPerMinuteToRejectChannel are reported as chronically bad
contacts, ordered worst-first for the shared rejection budget.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01TRpqXm5P93dw7q3f6wpWgQ"
```

---

### Task 4: `RELAX_repair_steps`

Partitions each channel at {boundary starts} U {its own step samples} and detrends each sub-segment in isolation, which removes both the drift and the step offset in one operation.

**Files:**
- Create: `RELAX_repair_steps.m`
- Create: `tests/test_RELAX_repair_steps.m`
- Create: `tests/test_RELAX_repair_steps_boundaries.m`

**Interfaces:**
- Consumes: `steps` from `RELAX_detect_steps` (Task 3); `nt_detrend` from `vendor/noisetools` (Task 1); `RELAX_get_segments`.
- Produces: `[EEG, info] = RELAX_repair_steps(EEG, steps, 'Name', value, ...)`. `info` has `.available`, `.status`, `.params`, `.nChannelsRepaired`, `.nSubSegments`, `.subSegmentsPerChannel` `[nbchan x 1]`, `.treatmentCounts` (struct with scalar fields `order3oa`, `order3`, `order2`, `order1`, `meanOnly`, `fallback`), `.minSubSegmentDuration`, `.medianSubSegmentDuration`, `.fractionSamplesMeanOnly`, `.lowFreqVarianceRemoved` `[nbchan x 1]`. Task 5 consumes all of it.

- [ ] **Step 1: Write the failing behaviour test**

Create `tests/test_RELAX_repair_steps.m`:

```matlab
function test_RELAX_repair_steps()
% Step repair: removes steps and drift, preserves signal above the cutoff, degrades safely.

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

    fprintf('test_RELAX_repair_steps: PASS\n');
end
```

- [ ] **Step 2: Write the failing boundary test**

Create `tests/test_RELAX_repair_steps_boundaries.m`:

```matlab
function test_RELAX_repair_steps_boundaries()
% Repair never fits across a boundary: each segment is processed as if standalone.

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

    fprintf('test_RELAX_repair_steps_boundaries: PASS\n');
end
```

- [ ] **Step 3: Run both tests to verify they fail**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_repair_steps"
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_repair_steps_boundaries"
```

Expected: both FAIL with `Unrecognized function or variable 'RELAX_repair_steps'`.

- [ ] **Step 4: Write the repairer**

Create `RELAX_repair_steps.m`:

```matlab
% RELAX_repair_steps() - Repair electrode-pop steps by per-sub-segment robust detrending.
%
% Usage:
%   >> [EEG, info] = RELAX_repair_steps(EEG, steps);
%   >> [EEG, info] = RELAX_repair_steps(EEG, steps, 'MaxCutoffHz', 1.0);
%
% Each channel is partitioned at the union of {inter-boundary segment starts} and {that
% channel's own step samples}. Every sub-segment is then detrended in complete isolation -
% no padding, no overlap across a cut.
%
% That isolation is the whole mechanism. nt_regw's returned trend CARRIES THE WEIGHTED
% MEAN (nt_regw.m lines 31 and 43/58: `mn = y - nt_demean(y,w)` is added back into `z`),
% so `x - z` has zero weighted mean over the fit interval. Detrending each inter-step
% sub-segment therefore removes that sub-segment's offset as a side effect, and a separate
% step-subtraction pass (nt_destep) would find nothing left to do.
%
% The step itself is never fitted: it sits ON a partition boundary, not inside a fit
% window.
%
% Inputs:
%   EEG   - EEGLAB EEG structure (uses .data, .nbchan, .pnts, .srate, .event)
%   steps - struct from RELAX_detect_steps. Pass one with empty .sample to apply
%            drift-only detrending with no step cuts.
%   Name-value pairs (all optional):
%     MaxPolyOrder      - highest polynomial order used (default 3)
%     MaxCutoffHz       - the invariant: this stage never attenuates above this frequency
%                          (default 1.0, just under the pipeline's own 1 Hz highpass, so
%                          the stage can only remove signal the highpass was going to
%                          remove anyway)
%     MaxWindowSeconds  - nt_detrend's overlap-add window for long sub-segments (default
%                          2.0, measured -3 dB 0.841 Hz at order 3). Keeps the cutoff from
%                          sliding below 0.5 Hz on a long segment.
%     DetrendThresh     - nt_detrend outlier threshold in sd (default 3)
%     DetrendIterations - nt_detrend reweighting iterations (default 3)
%     MergeSamples      - a step within this many samples of a segment edge is merged into
%                          it rather than creating a sliver sub-segment (default 20)
%
% Outputs:
%   EEG  - same structure with .data repaired
%   info - .available, .status, .params, .nChannelsRepaired, .nSubSegments,
%          .subSegmentsPerChannel, .treatmentCounts (order3oa/order3/order2/order1/
%          meanOnly/fallback), .minSubSegmentDuration, .medianSubSegmentDuration,
%          .fractionSamplesMeanOnly, .lowFreqVarianceRemoved
%
% TREATMENT LADDER. For a sub-segment of duration T, the largest order n <= MaxPolyOrder
% satisfying c(n)*n/(2T) <= MaxCutoffHz is used, with c = [1.69 1.24 1.09] MEASURED for
% orders 1-3. (The familiar n/(2T) rule of thumb underestimates the true -3 dB point by up
% to 70% at order 1, which is why the constants are here rather than the heuristic.) At the
% defaults that is:
%
%     T > 2.000 s           order 3, wsize = 2.0 s overlap-add     0.84 Hz
%     1.635 <= T <= 2.000   order 3, single window                 0.82-1.00 Hz
%     1.240 <= T < 1.635    order 2, single window                 0.76-1.00 Hz
%     0.845 <= T < 1.240    order 1, single window                 0.68-1.00 Hz
%     T < 0.845 s           robust mean subtraction only           DC only
%
% HARD vs SOFT SUBDIVISION. Cuts at boundaries and steps are HARD: the signal genuinely
% discontinues there. Subdividing a long sub-segment to respect MaxWindowSeconds is SOFT,
% and goes through nt_detrend's own wsize overlap-add (50% overlap, triangular weighting,
% continuous reconstructed trend). Hard-cutting a smooth stretch every 2 s would MANUFACTURE
% step discontinuities at 2 s intervals - exactly the input this stage exists to remove.
%
% Channels with no detected steps are still partitioned by boundaries and still detrended,
% so every channel receives the same drift treatment and the stage introduces no systematic
% difference between repaired and unrepaired channels.

function [EEG, info] = RELAX_repair_steps(EEG, steps, varargin)

    p = inputParser;
    p.addParameter('MaxPolyOrder', 3);
    p.addParameter('MaxCutoffHz', 1.0);
    p.addParameter('MaxWindowSeconds', 2.0);
    p.addParameter('DetrendThresh', 3);
    p.addParameter('DetrendIterations', 3);
    p.addParameter('MergeSamples', 20);
    p.parse(varargin{:});
    opt = p.Results;

    C = [1.69 1.24 1.09];   % measured f_3dB = C(n)*n/(2T); see the spec's response tables

    counts = struct('order3oa', 0, 'order3', 0, 'order2', 0, 'order1', 0, ...
                    'meanOnly', 0, 'fallback', 0);
    info = struct('available', false, 'status', '', 'params', opt, ...
                  'nChannelsRepaired', 0, 'nSubSegments', 0, ...
                  'subSegmentsPerChannel', zeros(EEG.nbchan,1), ...
                  'treatmentCounts', counts, ...
                  'minSubSegmentDuration', NaN, 'medianSubSegmentDuration', NaN, ...
                  'fractionSamplesMeanOnly', 0, ...
                  'lowFreqVarianceRemoved', zeros(EEG.nbchan,1));

    if isempty(which('nt_detrend'))
        info.status = 'nt_detrend not on the path - call RELAX_add_vendored_noisetools()';
        return
    end

    segments   = RELAX_get_segments(EEG);
    wsize      = 2 * floor(round(opt.MaxWindowSeconds * EEG.srate) / 2);  % samples, even
    origClass  = class(EEG.data);
    data       = double(EEG.data);
    repaired   = data;

    allDurations  = [];
    meanOnlySamps = 0;

    for ch = 1:EEG.nbchan
        % --- build this channel's partition ---
        chSteps = steps.sample(steps.channel == ch);
        cuts    = segments(:,1);
        for k = 1:numel(chSteps)
            % merge a step that sits on top of a segment edge rather than slivering it
            if all(abs(chSteps(k) - segments(:,1)) > opt.MergeSamples) && ...
               all(abs(chSteps(k) - (segments(:,2)+1)) > opt.MergeSamples)
                cuts(end+1) = chSteps(k); %#ok<AGROW>
            end
        end
        cuts = unique(cuts);

        subStarts = cuts;
        subEnds   = [cuts(2:end) - 1; EEG.pnts];
        % a cut may fall inside a later boundary segment; clip each sub-segment to its own
        for k = 1:numel(subStarts)
            s = find(segments(:,1) <= subStarts(k) & segments(:,2) >= subStarts(k), 1);
            subEnds(k) = min(subEnds(k), segments(s,2));
        end
        keep = subEnds >= subStarts;
        subStarts = subStarts(keep);
        subEnds   = subEnds(keep);

        info.subSegmentsPerChannel(ch) = numel(subStarts);

        for k = 1:numel(subStarts)
            a = subStarts(k); b = subEnds(k);
            x = data(ch, a:b)';                 % nt_detrend wants time on dim 1
            T = numel(x) / EEG.srate;
            allDurations(end+1,1) = T;          %#ok<AGROW>

            order = chooseOrder(T, C, opt.MaxPolyOrder, opt.MaxCutoffHz);

            if any(isnan(x)) || (max(x) - min(x)) == 0
                % A NaN would propagate through the whole fit; a constant stretch has no
                % trend to find and a zero scale for the reweighting.
                order = 0;
            end

            useOverlapAdd = (order == opt.MaxPolyOrder) && (numel(x) > wsize) && (wsize > 0);

            done = false;
            while ~done
                if order == 0
                    y = robustDemean(x, opt.DetrendThresh, opt.DetrendIterations);
                    counts.meanOnly = counts.meanOnly + 1;
                    meanOnlySamps = meanOnlySamps + numel(x);
                    done = true;
                else
                    try
                        if useOverlapAdd
                            y = nt_detrend(x, order, [], [], opt.DetrendThresh, ...
                                           opt.DetrendIterations, wsize);
                            counts.order3oa = counts.order3oa + 1;
                        else
                            y = nt_detrend(x, order, [], [], opt.DetrendThresh, ...
                                           opt.DetrendIterations);
                            switch order
                                case 3, counts.order3 = counts.order3 + 1;
                                case 2, counts.order2 = counts.order2 + 1;
                                case 1, counts.order1 = counts.order1 + 1;
                            end
                        end
                        if any(isnan(y)); error('RELAX:nanFromDetrend', 'nt_detrend returned NaN'); end
                        done = true;
                    catch
                        % nt_detrend's degenerate-window failure is a bare error('!') with
                        % no identifier. Step DOWN the ladder rather than killing a
                        % 257-channel SLURM task hours into a run.
                        counts.fallback = counts.fallback + 1;
                        useOverlapAdd = false;
                        order = order - 1;
                        if order < 1; order = 0; end
                    end
                end
            end

            repaired(ch, a:b) = y';
        end
    end

    removed = data - repaired;
    for ch = 1:EEG.nbchan
        v = var(data(ch, ~isnan(data(ch,:))));
        if v > 0
            info.lowFreqVarianceRemoved(ch) = var(removed(ch, ~isnan(removed(ch,:)))) / v;
        end
    end

    EEG.data = cast(repaired, origClass);

    info.available               = true;
    info.status                  = 'ok';
    info.nChannelsRepaired       = EEG.nbchan;
    info.nSubSegments            = sum(info.subSegmentsPerChannel);
    info.treatmentCounts         = counts;
    info.minSubSegmentDuration   = min(allDurations);
    info.medianSubSegmentDuration = median(allDurations);
    info.fractionSamplesMeanOnly = meanOnlySamps / (EEG.nbchan * EEG.pnts);
end

function order = chooseOrder(T, C, maxOrder, maxCutoffHz)
% Largest order n <= maxOrder whose measured -3 dB point stays at or below maxCutoffHz.
% 0 means "no polynomial fits in this duration" - robust mean subtraction only.
    order = 0;
    for n = 1:min(maxOrder, numel(C))
        if C(n) * n / (2*T) <= maxCutoffHz
            order = n;
        end
    end
end

function y = robustDemean(x, thresh, niter)
% Iteratively reweighted mean subtraction: the order-0 rung of the ladder, and the
% fallback when nt_detrend cannot fit. Written out rather than calling nt_detrend(x,0,...)
% because order 0 builds an empty polynomial basis, which nt_regw cannot regress on.
    finite = isfinite(x);
    if ~any(finite)
        y = x;
        return
    end
    w = double(finite);
    for k = 1:niter
        m = sum(w(finite) .* x(finite)) / sum(w(finite));
        r = x - m;
        s = sqrt(sum(w(finite) .* r(finite).^2) / sum(w(finite)));
        if ~isfinite(s) || s == 0; break; end
        wNew = double(abs(r) < thresh * s) .* finite;
        if ~any(wNew); break; end   % never leave an all-zero weight vector behind
        w = wNew;
    end
    m = sum(w(finite) .* x(finite)) / sum(w(finite));
    y = x - m;
end
```

- [ ] **Step 5: Run both tests to verify they pass**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_repair_steps"
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_RELAX_repair_steps_boundaries"
```

Expected: both PASS.

If the 2 Hz case of the spectral-fidelity test fails, check `chooseOrder` - the comparison must be `<=` against `MaxCutoffHz` and `C` must be `[1.69 1.24 1.09]`, not the `n/(2T)` heuristic. If the boundary bit-identity assertion fails, a fit is running across a boundary: check the clipping loop that limits each `subEnds(k)` to its own segment's end.

- [ ] **Step 6: Run the full suite**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: all test files pass.

- [ ] **Step 7: Commit**

```bash
git add RELAX_repair_steps.m tests/test_RELAX_repair_steps.m tests/test_RELAX_repair_steps_boundaries.m
git commit -m "Add RELAX_repair_steps: per-sub-segment robust detrending

Partitions each channel at {boundary starts} U {its own step samples} and
detrends each sub-segment in isolation. Because nt_regw's returned trend carries
the weighted mean, this removes each sub-segment's offset as a side effect - so
there is no separate step-subtraction pass, and the step itself is never fitted
because it sits on a partition boundary rather than inside a fit window.

The treatment ladder picks the largest polynomial order whose MEASURED -3 dB
point stays at or below MaxCutoffHz, using c = [1.69 1.24 1.09] rather than the
n/(2T) rule of thumb, which underestimates by up to 70% at order 1. Every rung
lands in [0.68, 1.00] Hz, just under the pipeline's own highpass, so the stage
can only remove signal the highpass was going to remove anyway.

Long sub-segments subdivide SOFT via nt_detrend's overlap-add rather than by
further hard cuts, which would manufacture new step discontinuities every 2 s.
Each fit is wrapped so nt_detrend's bare error('!') on a degenerate window steps
down the ladder instead of killing a 257-channel SLURM task.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01TRpqXm5P93dw7q3f6wpWgQ"
```

---

### Task 5: Wire the stage into the SLURM wrapper

Inserts the stage, adds the config fields, spends the shared budget, records the covariates, and updates the docs.

**Files:**
- Modify: `RELAX_Wrapper_slurm.m` (insert after line 237; budget init; impedance call ~line 367; stale comment ~line 364; covariate fields before line 512)
- Modify: `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` (new `RELAX_cfg` fields)
- Modify: `CLAUDE.md`
- Modify: `docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md` (add the two config fields Task 3 introduced)

**Interfaces:**
- Consumes: `RELAX_detect_steps`, `RELAX_repair_steps`, `RELAX_channel_budget`, `RELAX_impedance_bad_channels`'s 4-argument form.
- Produces: `EEG.RELAX.PopScreen`, `EEG.RELAX.PopRepair`, `EEG.RELAX.ChannelBudget`, and the scalar `RELAXProcessingExtremeRejections.Pop*` covariates that reach the aggregate CSV.

**Note on where covariates go:** the aggregate row is built by `struct2table(epochedEEG.RELAXProcessingExtremeRejections, 'AsArray', true)` at `RELAX_Wrapper_slurm.m:513`. Only **scalars and cell arrays** may go on that struct - a `[nbchan x 1]` vector or a nested struct produces an unusable column. Full info structs live on `EEG.RELAX.*` instead, which is not table-ified.

- [ ] **Step 1: Add the config fields to the entry point**

In `RELAX_SET_PARAMETERS_AND_RUN_slurm.m`, next to the other `Perform_*` flags, add:

```matlab
%% Electrode-pop detection and repair (runs before any filtering, at native sample rate).
% A pop is an abrupt DC baseline step on one channel. The 1 Hz acausal Butterworth turns
% each one into a multi-second bipolar excursion in BOTH time directions, so repair has to
% happen before the filter sees it. See
% docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md
RELAX_cfg.Perform_pop_artifact_repair             = 1;
% Detector - the first four are FieldTrip ft_artifact_jump's own defaults
RELAX_cfg.PopDetect_MedianFilterOrder             = 9;
RELAX_cfg.PopDetect_ZThreshold                    = 20;
RELAX_cfg.PopDetect_MinAmplitudeMicrovolts        = 10;
RELAX_cfg.PopDetect_RefractoryMs                  = 20;
RELAX_cfg.PopDetect_SustainMs                     = 100;
RELAX_cfg.PopDetect_MinRiseFraction               = 0.5;
RELAX_cfg.PopDetect_StepsPerMinuteToRejectChannel = 30;
% Repair
RELAX_cfg.PopRepair_MaxPolyOrder                  = 3;
RELAX_cfg.PopRepair_MaxCutoffHz                   = 1.0;
RELAX_cfg.PopRepair_MaxWindowSeconds              = 2.0;
RELAX_cfg.PopRepair_DetrendThresh                 = 3;
RELAX_cfg.PopRepair_DetrendIterations             = 3;
```

- [ ] **Step 2: Initialise the budget and insert the stage in the wrapper**

In `RELAX_Wrapper_slurm.m`, immediately after line 237 (`EEG.allchan=EEG.chanlocs;`) and before the `%% Notch filter data:` block, insert:

```matlab
    % The shared channel-rejection budget is taken against the montage EEG.allchan just
    % recorded, so it must be created here, before any screen has removed a channel.
    if ~isfield(RELAX_cfg,'MaxProportionOfElectrodesThatCanBeDeleted')
        RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0.10;
    end
    EEG = RELAX_channel_budget(EEG, 'init', RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted);

    %% Detect and repair electrode pops:
    % A pop is an abrupt change in one channel's DC baseline. This runs BEFORE the notch
    % and band-pass filters, at the native sample rate, because a zero-phase highpass
    % converts a step into a multi-second bipolar excursion extending in BOTH time
    % directions - after which there is no step left to find and the damage has already
    % spread to data recorded before the pop occurred.
    %
    % Absent = off: this is new behaviour and must not switch itself on for a dataset that
    % has never been checked for pops.
    if ~isfield(RELAX_cfg,'Perform_pop_artifact_repair')
        RELAX_cfg.Perform_pop_artifact_repair = 0;
    end
    EEG.RELAXProcessingExtremeRejections.PopBasedChannelToReject={};
    if RELAX_cfg.Perform_pop_artifact_repair == 1
        [PopSteps, PopScreenInfo] = RELAX_detect_steps(EEG, ...
            'MedianFilterOrder',             getfielddef(RELAX_cfg,'PopDetect_MedianFilterOrder',9), ...
            'ZThreshold',                    getfielddef(RELAX_cfg,'PopDetect_ZThreshold',20), ...
            'MinAmplitudeMicrovolts',        getfielddef(RELAX_cfg,'PopDetect_MinAmplitudeMicrovolts',10), ...
            'RefractoryMs',                  getfielddef(RELAX_cfg,'PopDetect_RefractoryMs',20), ...
            'SustainMs',                     getfielddef(RELAX_cfg,'PopDetect_SustainMs',100), ...
            'MinRiseFraction',               getfielddef(RELAX_cfg,'PopDetect_MinRiseFraction',0.5), ...
            'StepsPerMinuteToRejectChannel', getfielddef(RELAX_cfg,'PopDetect_StepsPerMinuteToRejectChannel',30));
        EEG.RELAX.PopScreen = PopScreenInfo;

        if ~PopScreenInfo.available
            warning('Pop screen unavailable for %s (%s) - no pops detected or repaired.', ...
                FileName, PopScreenInfo.status);
        else
            fprintf('Pop screen: %d step(s) on %d channel(s); %d common-mode bin(s).\n', ...
                PopScreenInfo.nSteps, PopScreenInfo.nChannelsAffected, PopScreenInfo.nCommonModeSteps);

            % Population A: chronically bad contacts. Repair is meaningless on a channel
            % with no stable baseline, so it is rejected rather than repaired. badChannels
            % arrives ordered worst-first, which is what the budget's truncation needs.
            [PopChannelToReject, EEG] = RELAX_channel_budget(EEG, 'spend', 'pop', PopScreenInfo.badChannels);
            EEG.RELAXProcessingExtremeRejections.PopBasedChannelToReject = PopChannelToReject;
            if ~isempty(PopChannelToReject)
                fprintf('Pop screen: removing %d chronically popping channel(s): %s\n', ...
                    numel(PopChannelToReject), strjoin(PopChannelToReject, ', '));
                keep = ~ismember(PopSteps.label, PopChannelToReject);
                PopSteps = structfun(@(f) f(keep), PopSteps, 'UniformOutput', false);
                EEG = pop_select(EEG,'nochannel',PopChannelToReject);
                EEG = eeg_checkset( EEG );
                % RELAX_detect_steps returns channel INDICES; re-map them after pop_select
                [~, PopSteps.channel] = ismember(PopSteps.label, {EEG.chanlocs.labels});
            end

            [EEG, PopRepairInfo] = RELAX_repair_steps(EEG, PopSteps, ...
                'MaxPolyOrder',      getfielddef(RELAX_cfg,'PopRepair_MaxPolyOrder',3), ...
                'MaxCutoffHz',       getfielddef(RELAX_cfg,'PopRepair_MaxCutoffHz',1.0), ...
                'MaxWindowSeconds',  getfielddef(RELAX_cfg,'PopRepair_MaxWindowSeconds',2.0), ...
                'DetrendThresh',     getfielddef(RELAX_cfg,'PopRepair_DetrendThresh',3), ...
                'DetrendIterations', getfielddef(RELAX_cfg,'PopRepair_DetrendIterations',3));
            EEG.RELAX.PopRepair = PopRepairInfo;
            if ~PopRepairInfo.available
                warning('Pop repair unavailable for %s (%s) - steps detected but NOT repaired.', ...
                    FileName, PopRepairInfo.status);
            else
                fprintf(['Pop repair: %d sub-segment(s), median %.2f s, min %.2f s; ' ...
                         '%.1f%% of samples mean-only, %d fallback(s).\n'], ...
                    PopRepairInfo.nSubSegments, PopRepairInfo.medianSubSegmentDuration, ...
                    PopRepairInfo.minSubSegmentDuration, ...
                    100*PopRepairInfo.fractionSamplesMeanOnly, ...
                    PopRepairInfo.treatmentCounts.fallback);
            end

            % Per-file covariates. RELAX thresholds are computed from within-recording
            % statistics, so a poppy session gets cleaned harder than a quiet one in a way
            % that tracks its electrode quality. A downstream model predicting treatment
            % response from session-to-session EEG change could learn the cleaning instead
            % of the neural signal, so how hard this stage worked has to be recorded per
            % file - see CLAUDE.md, "Saving covariates for downstream ML classification
            % tasks". Scalars and cells only: this struct is struct2table'd below.
            EEG.RELAXProcessingExtremeRejections.PopNSteps                = PopScreenInfo.nSteps;
            EEG.RELAXProcessingExtremeRejections.PopNChannelsAffected     = PopScreenInfo.nChannelsAffected;
            EEG.RELAXProcessingExtremeRejections.PopMedianStepAmplitude   = median(abs(PopSteps.amplitude));
            EEG.RELAXProcessingExtremeRejections.PopNChannelsRejected     = numel(PopChannelToReject);
            EEG.RELAXProcessingExtremeRejections.PopNCommonModeSteps      = PopScreenInfo.nCommonModeSteps;
            if PopRepairInfo.available
                EEG.RELAXProcessingExtremeRejections.PopMedianSubSegmentDuration = PopRepairInfo.medianSubSegmentDuration;
                EEG.RELAXProcessingExtremeRejections.PopFractionSamplesMeanOnly  = PopRepairInfo.fractionSamplesMeanOnly;
                EEG.RELAXProcessingExtremeRejections.PopMedianLowFreqVarianceRemoved = median(PopRepairInfo.lowFreqVarianceRemoved);
                EEG.RELAXProcessingExtremeRejections.PopNOrder3OA  = PopRepairInfo.treatmentCounts.order3oa;
                EEG.RELAXProcessingExtremeRejections.PopNOrder3    = PopRepairInfo.treatmentCounts.order3;
                EEG.RELAXProcessingExtremeRejections.PopNOrder2    = PopRepairInfo.treatmentCounts.order2;
                EEG.RELAXProcessingExtremeRejections.PopNOrder1    = PopRepairInfo.treatmentCounts.order1;
                EEG.RELAXProcessingExtremeRejections.PopNMeanOnly  = PopRepairInfo.treatmentCounts.meanOnly;
                EEG.RELAXProcessingExtremeRejections.PopNFallback  = PopRepairInfo.treatmentCounts.fallback;
            end
        end
    else
        fprintf('Skipping electrode-pop repair (RELAX_cfg.Perform_pop_artifact_repair=0).\n');
    end
```

At the very bottom of `RELAX_Wrapper_slurm.m`, after the main function's `end`, add the small helper it uses:

```matlab
function v = getfielddef(s, name, default)
% Field value, or a default when the field is absent. Keeps the long name-value lists
% above readable without an isfield test per parameter.
    if isfield(s, name)
        v = s.(name);
    else
        v = default;
    end
end
```

- [ ] **Step 3: Make the impedance screen spend from the shared budget**

In `RELAX_Wrapper_slurm.m`, replace the stale comment at lines 364-366:

```matlab
        % Capped by the same budget RELAX_excluding_channels_and_epoching spends, so the
        % two channel-rejection stages cannot between them delete more of the montage
        % than MaxProportionOfElectrodesThatCanBeDeleted allows.
```

with:

```matlab
        % Spends from the shared budget created at the top of this loop, so the pop,
        % impedance and PREP screens cannot between them delete more of the montage than
        % MaxProportionOfElectrodesThatCanBeDeleted allows. (Before the budget existed
        % each screen independently capped itself at that proportion, so three screens
        % could remove three times it.) RELAX_excluding_channels_and_epoching keeps its
        % own later, independent cap - see RELAX_channel_budget.m for why.
```

Then change the call at line 367 to pass the remaining budget:

```matlab
        popBudget = EEG.RELAX.ChannelBudget;
        [ImpedanceBasedChannelToReject, ImpedanceScreenInfo] = RELAX_impedance_bad_channels( ...
            EEG.chanlocs, RELAX_cfg.MaxImpedanceKOhm, ...
            RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted, ...
            max(0, popBudget.maxRejectable - popBudget.nSpent));
```

And immediately after the existing `EEG.RELAX.ImpedanceScreen = ImpedanceScreenInfo;` line, record the spend so PREP sees the reduced remainder:

```matlab
        [~, EEG] = RELAX_channel_budget(EEG, 'spend', 'impedance', ImpedanceBasedChannelToReject);
```

- [ ] **Step 4: Run the full suite**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: all test files pass.

- [ ] **Step 5: Syntax-check the two modified scripts**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); checkcode('RELAX_Wrapper_slurm.m'); checkcode('RELAX_SET_PARAMETERS_AND_RUN_slurm.m'); disp('CHECKCODE DONE')"
```

Expected: `CHECKCODE DONE`, with no `Parse error`. Pre-existing style warnings are fine; a parse error is not.

- [ ] **Step 6: Run the integration smoke test**

```bash
bash tests/run_integration_smoke.sh
```

Expected: the same outcome as before this task. If this script does not exist or does not apply, skip and note it.

- [ ] **Step 7: Update CLAUDE.md**

Three edits.

(a) In the "Pipeline architecture" stage list, insert a new stage **between stage 1 (Load & prep) and stage 2 (Notch filter)** - the stage rejects channels but runs before any filtering, so it does not belong beside the other channel-rejection screens at stage 4. Renumber the following stages:

```markdown
2. **Electrode-pop detection & repair** (this fork's addition) - gated by
   `RELAX_cfg.Perform_pop_artifact_repair`. `RELAX_detect_steps.m` finds abrupt DC baseline
   steps with a median-filtered derivative plus a robust z threshold, an absolute microvolt
   floor and a rise-time criterion; `RELAX_repair_steps.m` partitions each channel at
   {boundary starts} U {its own step samples} and detrends each sub-segment in isolation,
   which removes the step offset and the drift in one operation. Runs before any filtering
   because a zero-phase highpass turns a step into a multi-second bipolar excursion in both
   time directions. Chronically popping channels are rejected instead of repaired, capped by
   the shared budget in `RELAX_channel_budget.m`. Design:
   `docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md`.
```

(b) In the "External dependencies" section, add alongside the zapline-plus paragraph:

```markdown
NoiseTools is **not** an external dependency: the 9-file runtime closure of `nt_detrend` is
vendored under `vendor/noisetools/` (see `vendor/noisetools/VENDORING.md`). Unlike
`vendor/zapline-plus/` this is a deliberately PARTIAL mirror, because nothing shadows
NoiseTools - it is not an EEGLAB plugin. Entry points must call
`RELAX_add_vendored_noisetools()`.
```

(c) In the "Naming convention for pipeline variants" section, add:

```markdown
The electrode-pop stage exists **only** in `RELAX_Wrapper_slurm.m`. `RELAX_Wrapper.m` and
`RELAX_Wrapper_bart.m` have no pop detection or repair, and no shared channel-rejection
budget - their screens still cap themselves independently at
`MaxProportionOfElectrodesThatCanBeDeleted` each. Porting the stage means adding the budget
init, the stage block, and the impedance screen's 4th argument to those wrappers too.
```

- [ ] **Step 8: Add the two missing config fields to the spec**

The spec's config block omits `SustainMs` and `MinRiseFraction`, which its own detection
step 5 ("require the level change to be sustained", "rise-time criterion") requires. In
`docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md`, in the Configuration
block, add after the `PopDetect_RefractoryMs` line:

```matlab
RELAX_cfg.PopDetect_SustainMs                     = 100;
RELAX_cfg.PopDetect_MinRiseFraction               = 0.5;
```

And in the detection description, after the numbered step 5, add:

```markdown
`SustainMs` (default 100) sets the window either side used to measure the pre- and
post-step level - the level change must persist across it, which is what makes this a step
detector rather than a spike detector. `MinRiseFraction` (default 0.5) is the rise-time
criterion itself: the fraction of the total level change that must occur in the single
largest derivative sample. A genuine pop delivers most of its change in one sample; a ramp
spreads it over thousands and fails.
```

- [ ] **Step 9: Commit**

```bash
git add RELAX_Wrapper_slurm.m RELAX_SET_PARAMETERS_AND_RUN_slurm.m CLAUDE.md docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md
git commit -m "Wire electrode-pop repair into the SLURM wrapper

Inserts the stage between the channel deletion at line 237 and the notch filter,
so it runs before any filtering at the native sample rate. Initialises the shared
channel-rejection budget against the same montage snapshot EEG.allchan records,
and makes the impedance screen spend from it rather than applying its own
independent cap.

Per-file covariates land on RELAXProcessingExtremeRejections, which reaches the
aggregate CSV: how hard this stage worked scales with how poppy the session was,
so it has to be recorded or a downstream model could learn the cleaning instead
of the neural signal.

Documents in CLAUDE.md that the stage exists only in the SLURM wrapper, and adds
the two detector parameters (SustainMs, MinRiseFraction) the spec's config block
had omitted.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01TRpqXm5P93dw7q3f6wpWgQ"
```

---

### Task 6: Real-file validation and the deferred ringing report

Runs the stage on the five surveyed files and produces the measurement promised by design decision 7, which determines whether guard-band interpolation is added.

**Files:**
- Create: `validate_pop_repair.m` (repo root, alongside `find_missing_files.m` and `recover_metrics.m`)
- Create: `docs/superpowers/reports/2026-09-14-pop-repair-validation.md`

**Interfaces:**
- Consumes: everything from Tasks 1-5.
- Produces: a written report. No code depends on it.

- [ ] **Step 1: Write the validation script**

Create `validate_pop_repair.m`:

```matlab
function validate_pop_repair(outFile)
% VALIDATE_POP_REPAIR - Measure what the electrode-pop stage does on real files.
%
% Answers the question deferred by decision 7 of
% docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md: guard-band samples are
% currently left as recorded, so how much ringing survives at a repaired step site after
% the 1 Hz highpass? If it is large, guard-band interpolation gets added.
%
% Run on a compute node, never a login node:
%   srun -p sackler -c 4 --mem=32G matlab -nodisplay -nosplash -batch \
%     "addpath(pwd); validate_pop_repair('pop_validation.mat')"

    if nargin < 1; outFile = 'pop_validation.mat'; end

    eeglab nogui;
    RELAX_add_vendored_noisetools();

    dataDir = '~/TMS_response_variability/store/mdd_dlpfc/preprocessed_Aug26/tms_artifact_removed';
    files = {'m159_dlpfc_day1_treatment_20230227_123711_interTrainIntervals.set', ...
             'm497_dlpfc_day1_treatment_20240728_102919_interTrainIntervals.set', ...
             'm208_dlpfc_day1_treatment_20230717_120446_interTrainIntervals.set', ...
             'm583_dlpfc_day1_treatment_20250804_124646_interTrainIntervals.set', ...
             'm155_dlpfc_day1_treatment_20230313_012312_interTrainIntervals.set'};

    results = struct([]);
    for k = 1:numel(files)
        fprintf('\n===== %s\n', files{k});
        EEG = pop_loadset('filename', files{k}, 'filepath', dataDir);

        [steps, dInfo] = RELAX_detect_steps(EEG);
        fprintf('detected %d steps on %d channels; %d flagged as bad contacts\n', ...
            dInfo.nSteps, dInfo.nChannelsAffected, numel(dInfo.badChannels));

        rawEEG = EEG;
        [EEG, rInfo] = RELAX_repair_steps(EEG, steps);

        % --- the ringing measurement ---
        % Highpass both the raw and the repaired data exactly as the pipeline will, then
        % compare the energy in a 2 s window centred on each repaired step site against a
        % matched step-free control window on the same channel.
        hpRaw  = pop_eegfiltnew(rawEEG, 1, []);
        hpRep  = pop_eegfiltnew(EEG,    1, []);

        halfWin = round(1.0 * EEG.srate);
        ringRaw = []; ringRep = []; ctrlRep = [];
        for s = 1:numel(steps.sample)
            c = steps.channel(s); i = steps.sample(s);
            if i-halfWin < 1 || i+halfWin > EEG.pnts; continue; end
            ringRaw(end+1) = std(hpRaw.data(c, i-halfWin:i+halfWin)); %#ok<AGROW>
            ringRep(end+1) = std(hpRep.data(c, i-halfWin:i+halfWin)); %#ok<AGROW>
            j = mod(i + round(EEG.pnts/2), EEG.pnts - 2*halfWin) + halfWin + 1;
            ctrlRep(end+1) = std(hpRep.data(c, j-halfWin:j+halfWin)); %#ok<AGROW>
        end

        r = struct();
        r.file                  = files{k};
        r.detect                = dInfo;
        r.repair                = rInfo;
        r.ringingRawMedian      = median(ringRaw);
        r.ringingRepairedMedian = median(ringRep);
        r.controlMedian         = median(ctrlRep);
        r.residualRatio         = median(ringRep) / median(ctrlRep);
        results = [results, r]; %#ok<AGROW>

        fprintf(['ringing at step sites after 1 Hz HP: raw %.1f uV -> repaired %.1f uV; ' ...
                 'step-free control %.1f uV; residual ratio %.2fx\n'], ...
            r.ringingRawMedian, r.ringingRepairedMedian, r.controlMedian, r.residualRatio);
        fprintf('ladder: oa=%d o3=%d o2=%d o1=%d mean=%d fallback=%d; %.1f%% samples mean-only\n', ...
            rInfo.treatmentCounts.order3oa, rInfo.treatmentCounts.order3, ...
            rInfo.treatmentCounts.order2, rInfo.treatmentCounts.order1, ...
            rInfo.treatmentCounts.meanOnly, rInfo.treatmentCounts.fallback, ...
            100*rInfo.fractionSamplesMeanOnly);
    end

    save(outFile, 'results');
    fprintf('\nSaved %s\n', outFile);
end
```

- [ ] **Step 2: Run it on a compute node**

```bash
module load matlab
cd /home/nga4004/gitRepos/RELAX-GEDAI
srun -p sackler -c 4 --mem=32G matlab -nodisplay -nosplash -batch \
  "addpath(pwd); validate_pop_repair('$CLAUDE_JOB_DIR/tmp/pop_validation.mat')" \
  2>&1 | tee "$CLAUDE_JOB_DIR/tmp/pop_validation.log"
```

Expected: five `=====` blocks, each printing detection counts, a residual ratio, and a ladder breakdown. A `residualRatio` near 1.0 means the repaired step site looks like ordinary data after the highpass; substantially above 1.0 means ringing survives.

- [ ] **Step 3: Write the report**

Create `docs/superpowers/reports/2026-09-14-pop-repair-validation.md` containing, from the log:

1. A table of the five files: steps detected, channels affected, channels rejected as bad contacts, `fractionSamplesMeanOnly`, ladder rung counts, fallback count.
2. The ringing table: raw / repaired / control median uV and `residualRatio` per file.
3. **The decision on guard-band interpolation**, stated explicitly. If every `residualRatio` is below about 1.5, record that the guard band is being left as recorded and why. If any file is materially above that, recommend adding interpolation and say which files drove the call.
4. A check on `StepsPerMinuteToRejectChannel = 30` - the one default not borrowed from FieldTrip. Compare the per-channel rate distribution against the threshold and say whether 30 still sits in an empty gap on real data.
5. Any detector false positives visible in the counts, compared against the exploratory survey numbers in the spec's "Measured properties" table (that survey had no rise-time criterion, so the production counts should be **lower**; if they are higher, something is wrong).

- [ ] **Step 4: Commit**

```bash
git add validate_pop_repair.m docs/superpowers/reports/2026-09-14-pop-repair-validation.md
git commit -m "Add pop-repair validation script and report

Measures, on the five surveyed files, what the stage actually does: detection
counts against the exploratory survey, ladder rung distribution, and the
residual ringing at repaired step sites after the 1 Hz highpass compared against
matched step-free control windows.

That last measurement is the one deferred by decision 7 of the design spec: the
guard band is currently left as recorded, and this is the evidence for whether
interpolation needs adding.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01TRpqXm5P93dw7q3f6wpWgQ"
```
