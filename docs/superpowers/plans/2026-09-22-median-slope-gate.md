# Median-Slope Gate for the Step Detector — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans to implement this plan task-by-task, with the
> agent roles defined under "Execution protocol" below. Steps use checkbox (`- [ ]`) syntax
> for tracking.

**Goal:** Stop `RELAX_detect_steps` from accepting slow drift and slow waves as electrode
pops, by gating candidates on how *fast* the median level moves across the transition
allowance rather than on how *far* it moved over a 250 ms lag — and decide the thresholds
from real data before they ship.

**Architecture:** The current statistic `D` is already a difference of a centred median
filter at a lag of `TransitionMs + WindowMs`. This plan adds a second reading of the *same*
filter at a lag of `TransitionMs` alone (`Dslope`), scores it as an excess over the local
baseline slope in units of local robust sd, and makes that the local gate; the channel-wide
gate moves to `Dslope`; the absolute floor stays on `D`. Bracketing, repair, ranking, eye
flagging and event grouping are untouched. The sample-level instability criterion, whose
statistic disappears with the old local gate, is removed first as a self-contained deletion.
A survey runs the pinned old detector and the new one in one process on the five survey
files, and the thresholds are chosen from that before any default is committed.

**Tech Stack:** MATLAB (`module load matlab`) on the SCU SLURM cluster, EEGLAB 2026.0.0 at
`~/tools/eeglab2026.0.0`. Tests via `sbatch --wait run_one_test.sh <test>` and
`sbatch --wait run_all_tests.sh`; stage validation via `sbatch run_validate_pop_repair.sh`.

**Spec:** `docs/superpowers/specs/2026-09-22-median-slope-gate-design.md`. It extends the
detector described in `docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md`
(addendum 2026-09-17, "Detection") and supersedes that addendum's gating statistic.

## Global Constraints

- **This plan makes no commits.** The working tree already carries the uncommitted hold-drift
  and reference-pop work, and that is fine — the user's plan is a **single commit** covering
  the whole electrode-pop detection/repair line of work, made once all of it is done. Build on
  the dirty tree: never `git stash`, never `git checkout --` a tracked file, never
  `git commit`. Each task ends with an evaluator checkpoint that reads the diff instead.
- **Because nothing is committed, the diff is the only record.** Every task's checkpoint writes
  its findings — suite counts, log paths, decisions, deferrals — into the task notes kept in
  this session; Tasks 3-5 write theirs into deliverable files
  (`docs/superpowers/reports/2026-09-23-slope-gate-survey.md`). Nothing that matters lives only
  in a scrollback.
- Never run MATLAB on a login node (`scu-login*`) and never through the MATLAB MCP tool.
  Every MATLAB invocation goes through `sbatch` on the `scu-cpu` partition
  (CLAUDE.local.md). `mkdir -p logs` once before the first `sbatch`.
- **No threshold ships before Task 4C.** Task 3A sets one detection-only operating point in
  all four places (function, wrapper, entry script, `validate_pop_repair`). Its `TransitionMs`
  4 and `ZThreshold`/`ChannelZThreshold` 4 are the provisional middle of the grid; Task 4C
  replaces them with the values the scored survey supports. A default the survey did not
  support is never committed.
- `SlopeNoiseFloorMicrovolts` stays `0.1` and is **not** swept, but at lags of 1-10 samples it
  may do the work the local spread was meant to do. Task 3A measures how often it binds
  (`info.floorBoundFraction`), and Task 3B Step 6 stops for a user decision if it binds on most
  samples.
- **Blind scoring is the user's** (2026-09-23). The evaluator never shows the user a
  candidate's config, lag, z, stratum or weight before `scores.csv` is complete, and looks at
  rendered candidates only to check rendering quality.
- **All other evaluation happens in this session** (standing ruling 2026-09-21): reading logs,
  viewing figures, proposing thresholds from the pre-registered rule, deciding pass/fail. The
  user decides the operating point and the screens. No evaluator subagent is ever dispatched.
- Survey and validation run on the four Sep26 files named in the Revision note before Task 3A,
  not the Aug26 `tms_artifact_removed` files. The one exception to the sonnet-only rule for review agents is
  Task 6, where the user has chosen `complex-code-reviewer` (opus) for the final review.
- Follow the repo's function style: leading `%` help block, `inputParser` options, `info`
  struct with `.available` / `.status`, 0x1 empty columns (never a bare `[]`), 1-based
  inclusive sample indices. Prefer code that explains itself over comments (user CLAUDE.md).
- Tests are plain functions `tests/test_<name>.m` that `assert` and print `PASS`; the runner
  discovers them by name. Every test file seeds `rng`. **A test's expectation is never
  loosened to make an implementation pass**; the evaluator rules on any conflict.
- `tests/test_RELAX_bridged_channels.m` fails on this branch for reasons predating this plan
  ("cutoff frequencies must be within (0,1)"); the suite target is therefore **all but that
  one**. Record the pass count at the start of Task 1 and hold it.
- Always `git --no-pager` for any git command an agent runs: a pager invoked from inside
  MATLAB's `system()` hangs with no output.
- Do not modify `RELAX_Wrapper.m` or `RELAX_Wrapper_bart.m` — the pop stage exists only in
  `RELAX_Wrapper_slurm.m`.
- **The next production run must write to a fresh `RELAXProcessed/` directory.** Task 1
  removes a field from `RELAXProcessingExtremeRejections`, and the wrapper appends rows to
  reloaded aggregate `.mat` checkpoints. This is stated in `CLAUDE.md` by Task 5; nobody
  launches a production array from this branch before then.

---

## Design

### The statistic

For odd `W` with `h = floor(W/2)`, `movmedian(x, [W-1 0])` at index `i` is exactly the
centred median-filter output at `i - h`. Writing `mc` for that centred output:

```
mc(t)     = mPre(t + h)                  % an index shift of a value already computed
Dslope(t) = mc(t + g) - mc(t - g)        % lag 2g ~ TransitionMs
D(t)      = mPost(t + g) - mPre(t - g)   % unchanged: mc(t+g+h) - mc(t-g-h), lag TransitionMs + WindowMs
```

A lag-`L` difference of a median filter reports `A * min(1, L/T)` for a step of amplitude `A`
and rise `T`, and `s * L` for drift of slope `s`. Shortening the lag from `TransitionMs +
WindowMs` to `TransitionMs` is what buys drift rejection.

`Dslope` reads `x` over `t-g-h .. t+g+h`, strictly narrower than `D`'s `t-g-2h .. t+g+2h`, so
the existing `validMask` margins and `minSegLen` stay correct. The two median windows overlap
at the defaults (100 ms lag, 150 ms windows); this is harmless because a median filter is
edge-preserving — `mc` has already switched to the post-level once its centre passes the
transition — and it reduces the noise of the difference, because overlapping supports share
samples.

### The gate

```
b(t) = median(|Dslope|)      over a centred span of SlopeBaselineMs   (omitnan)
s(t) = 1.4826 * MAD(Dslope)  over the same span                       (omitnan)
z(t) = (|Dslope(t)| - b(t)) / max(s(t), SlopeNoiseFloorMicrovolts)

candidate  <=>  z > ZThreshold
            &   |Dslope| / channelScaleSlope > ChannelZThreshold
            &   |D| >= MinAmplitudeMicrovolts
```

- **Excess over a local baseline, not a ratio to variability alone.** Steady drift produces an
  almost-constant non-zero `Dslope` whose MAD is near zero; normalising only by variability
  would make `z` explode on the case this plan exists to reject. Subtracting `b` handles
  steady drift, dividing by `s` handles noisy drift.
- **A wide centred span, not the flanking windows.** `Dslope` at `u = t - g` equals
  `mc(t) - mc(t - 2g)`, which already contains half the step: windows flanking `t` are
  contaminated by the event being scored. Over ~2 s the contaminated samples are a ~15-20 %
  minority and both the median and the MAD have 50 % breakdown.
- **The floor stays on `D`.** `D`'s longer lag measures the settled level change, which is
  what "50 microvolts" means; keeping it also stops the new gate promoting small abrupt
  wiggles into repairable intervals.
- `channelScaleSlope` is `1.4826 * MAD(Dslope)` over locally-quiet samples, by the same
  two-pass construction as today: `loc` (the provisional local mask) still exists to select
  `valid & ~loc`; what disappears is its *export* as `info.unstableFraction`.

### What this trades away

A step whose rise `T` exceeds `2g` reports only `A * 2g/T`. At `TransitionMs = 100` that
attenuates the slow tail of this dataset's pops (rises run out to ~100 ms), which is why the
Task 3 sweep includes `TransitionMs` itself (100-200 ms).

**The risk to watch in the sweep:** a genuine step sitting inside a drifting stretch has its
own `b` raised by the drift around it. `SlopeBaselineMs` is the knob — too narrow and the step
contaminates its own baseline, too wide and the gate stops being local — and `ZThreshold` and
`SlopeBaselineMs` trade against each other on exactly this case. Task 2's "step riding on
drift" test is the red line; Task 3 reports, for every retained/lost step, `b` and `s` at the
step so the trade is visible rather than inferred.

> **2026-09-23:** the lag range above (100-200 ms) and the lost/retained reporting are
> superseded. The lag is now swept over 1-10 samples, and the survey scores candidates directly.
> See "Revision 2026-09-23" before Task 3A. At those lags a candidate run spans roughly the rise
> plus the lag rather than `2g - T`, so the next section holds only for the old lags.

### Run shape is unchanged

`Dslope` sits at the full amplitude for every `t` whose two read-points straddle the
transition, so a candidate run spans about `2g - T` samples (~100 ms at defaults).
`signedRuns`, `RefractoryMs`, `mid0`, `refineCandidate`, level re-measurement,
`plateauDrift`, `holdDrift`, `varianceRatio`, event grouping and the common-mode flag are all
untouched. `steps.z` and `steps.channelZ` keep their names and acquire the new meaning;
`steps.z` becomes a signed excess, so the row takes `max(z(rs:re))`, not `max(abs(...))`.

---

## Execution protocol

| Role | Agent | Does | Never does |
|---|---|---|---|
| **Agent** | `relax-implementer` (sonnet, effort medium — **effort high for Task 2**) | the whole task: writes the test edits exactly as the task specifies, submits with `sbatch --wait`, reports the Red log verbatim, then implements, then re-runs the task's test and the full suite and reports both logs | changes a test's expectation to make an implementation pass; commits anything |
| **Evaluator** | **this session only** (Opus, no subagent) | reads every log and every diff, diffs each test file against the plan's text, decides pass/fail, sends the agent back with precise findings, views the survey figures, adjudicates steps, picks the thresholds | trusts a "PASS" claim without reading the log; delegates any of this to an agent |

One agent per task does both the test and the implementation. The tester/implementer split
exists to stop an implementer quietly loosening a test it cannot pass; three things replace it
here, and **all three are required**:

1. **The test code is given verbatim in this plan.** The agent transcribes it; it does not
   design expectations. Task 2 Step 1's characterisation is the one place a test's *content* is
   chosen from evidence, and the evaluator chooses it, not the agent.
2. **Red before green is a reported artifact.** Every task has a step that runs the new test
   before the implementation exists and pastes the failure log. A task that never shows a Red
   log has not been done.
3. **The evaluator diffs every test file against this plan before accepting the task.** Any
   deviation from the text above — a loosened tolerance, a deleted assertion, a changed
   expectation — is a finding, not a judgement call.

Task 2 runs at **effort high**: it is index arithmetic on shifted median filters
(`mc(t) = mPre(t+h)`, `Dslope(t) = mc(t+g) - mc(t-g)`), where an off-by-`h` error produces a
detector that still runs, still finds steps, and is silently wrong. Tasks 1, 3, 4 and 5 are
mechanical or evaluator-led and run at medium.

Every subagent prompt contains: the task text copied from this plan, "MATLAB only via `sbatch`
on `scu-cpu`; never on a login node; never via the MATLAB MCP", "this plan commits nothing —
never run `git commit`, `git stash` or `git checkout --`", the exact commands, and "report the
full log path and the last 40 lines of the log verbatim". Subagents start cold.

Useful commands (all from `/home/nga4004/gitRepos/RELAX-GEDAI`):

```bash
mkdir -p logs
JOB=$(sbatch --wait --parsable run_one_test.sh test_RELAX_detect_steps); tail -40 logs/one_test_${JOB}.log
JOB=$(sbatch --wait --parsable run_all_tests.sh); tail -60 logs/run_all_tests_${JOB}.log
```

A passing single test ends with `ONE_TEST_PASS`; a failing one shows the MATLAB error and
`sbatch --wait` returns non-zero.

---

### Task 1: Retire the sample-level instability criterion

Pure deletion, done first so the gate swap in Task 2 lands on a detector with one fewer
statistic to carry. The criterion's statistic (`info.unstableFraction`, the fraction of
samples passing the **old local gate**) has no counterpart under the new gate, and the
telegraph-noise channels it was a safety net for already reject by step rate at the defaults —
which the existing test asserts on the line directly above the instability assertion.

**Files:**
- Modify: `RELAX_detect_steps.m` (help block lines 72-76 and 107; `inputParser` line 147;
  `info` initialiser line 160; assignment line 196)
- Modify: `RELAX_rank_pop_channels.m` (help lines 15-16 and 19; `U` line 43; `scores`
  line 44; `case 3` line 63)
- Modify: `RELAX_Wrapper_slurm.m` (config default line 349; pass-through line 453;
  `PopNRejectedByUnstable` line 541; covariate assignment line 621)
- Modify: `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` (line 366)
- Modify: `validate_pop_repair.m` (lines 173 and 194)
- Test: `tests/test_RELAX_detect_steps.m` (field list line 46-47; telegraph section
  lines 343-352)
- Test: `tests/test_RELAX_rank_pop_channels.m` (lines 39 and 46)
- Test: `tests/test_RELAX_protect_eye_channels.m` (fixtures lines 13 and 42)
- Test: `tests/assert_smoke_outcome.m` (covariate list line 33)
- Create: `slope_gate_survey/legacy/RELAX_detect_steps_level.m` (this detector, pinned for
  Task 3's comparison arm)

**Interfaces:**
- Produces: `RELAX_detect_steps` no longer accepts `MaxUnstableFractionToRejectChannel` and
  its `info` has no `unstableFraction` field. `RELAX_rank_pop_channels` scores channels on
  two criteria (`stepRatePerMinute / R`, `transitionFraction / F`) and emits only the
  `'steps/min exceeds'`, `'of samples in transition exceeds'` and `'flat: ...'` reasons.
  `EEG.RELAXProcessingExtremeRejections` no longer has `PopNChannelsRejectedByUnstable`.

- [ ] **Step 1: Record the suite baseline (Agent)**

```bash
mkdir -p logs
JOB=$(sbatch --wait --parsable run_all_tests.sh); tail -20 logs/run_all_tests_${JOB}.log
```

Report the `N/M test files passed` line verbatim. Expected: every file but
`test_RELAX_bridged_channels` passes. That count is the target for every later run.

- [ ] **Step 2: Write the failing tests (Agent)**

In `tests/test_RELAX_detect_steps.m`, remove `'unstableFraction'` from the field list so the
section reads:

```matlab
    for f = {'stepsPerChannel','stepRatePerMinute','transitionFraction','channelScaleD'}
        assert(isequal(size(info.(f{1})), [EEG.nbchan 1]), 'info.%s must be nbchan x 1', f{1});
    end
    assert(~isfield(info, 'unstableFraction'), ...
        'the sample-level instability criterion was retired on 2026-09-22');
```

Replace the telegraph section (the block beginning `% telegraph every 250 ms`) with:

```matlab
    % telegraph every 250 ms: every transition is abrupt and large, so the channel rejects on
    % step rate; the sample-level instability criterion it used to need was retired 2026-09-22
    x = randn(1, n60) * 5; sq = 300 * (mod(floor(tt / 250), 2) == 0); x = x + sq;
    [~, info] = RELAX_detect_steps(mkEEG(x));
    assert(isequal(info.badChannels, {'E1'}), 'telegraph channel must still be rejected');
    assert(contains(info.badChannelReason{1}, 'steps/min'), ...
        'telegraph channel must now reject by rate, got: %s', info.badChannelReason{1});
    assert(isequal(size(info.badChannels), size(info.badChannelReason)));
```

In `tests/test_RELAX_rank_pop_channels.m`, delete the assertion at line 39 that compares
`info3.unstableFraction`, keeping the `info3.nSteps == info.nSteps` half:

```matlab
    assert(info3.nSteps == info.nSteps);
```

and replace the line-46 assertion (which asserts an `'unstable'` reason) with a rate-based
one on the same fixture:

```matlab
    assert(isequal(info6.badChannels, {'E1'}) && contains(info6.badChannelReason{1}, 'steps/min'));
```

In `tests/test_RELAX_protect_eye_channels.m`, change both fixtures so no reason string
mentions instability — line 13's third entry and line 42's only entry become transition-
fraction reasons, which the function treats identically:

```matlab
                  'badChannelReason', {{'30.0 steps/min exceeds 6'; '12.0 steps/min exceeds 6'; '0.30 of samples in transition exceeds 0.02'}});
```

```matlab
    info4 = struct('badChannels', {{'E10'}}, 'badChannelReason', {{'0.30 of samples in transition exceeds 0.02'}});
```

In `tests/assert_smoke_outcome.m` line 33, delete `'PopNChannelsRejectedByUnstable',` from
the covariate name list.

- [ ] **Step 3: Run the tests to verify they fail (Agent)**

```bash
JOB=$(sbatch --wait --parsable run_one_test.sh test_RELAX_detect_steps); tail -40 logs/one_test_${JOB}.log
```

Expected: FAIL on `the sample-level instability criterion was retired on 2026-09-22` (the
field is still present). Report the log verbatim.

- [ ] **Step 4: Delete the criterion (Agent)**

`RELAX_detect_steps.m`: delete the `MaxUnstableFractionToRejectChannel` `addParameter` line,
the `'unstableFraction', zeros(nb,1), ...` entry from the `info` initialiser, the
`info.unstableFraction(ch) = mean(loc(valid));` assignment, the parameter's help paragraph
(lines 72-76), and the `.unstableFraction` mention in the `info` help at line 107. `loc`
stays — `gateCandidates` still returns it and still uses it to select quiet samples.

`RELAX_rank_pop_channels.m`: delete the `U = ...` line, drop the third column so

```matlab
    scores = [info.stepRatePerMinute / R, info.transitionFraction / F];
```

delete `case 3` from the reason `switch`, and update the help paragraph at lines 15-16 (which
describes the criterion) and the `info` argument list at line 19.

`RELAX_Wrapper_slurm.m`: delete the `'PopDetect_MaxUnstableFractionToRejectChannel', 0.25, ...`
default, the `'MaxUnstableFractionToRejectChannel', ...` pass-through argument, the
`PopNRejectedByUnstable = sum(contains(PopRejectReason, 'unstable'));` line, and the
`EEG.RELAXProcessingExtremeRejections.PopNChannelsRejectedByUnstable = ...` assignment.

`RELAX_SET_PARAMETERS_AND_RUN_slurm.m`: delete line 366.

`validate_pop_repair.m`: delete the `cfg.PopDetect_MaxUnstableFractionToRejectChannel = 0.25;`
line and the matching pass-through argument.

- [ ] **Step 5: Run the task's test and the full suite (Agent)**

```bash
JOB=$(sbatch --wait --parsable run_one_test.sh test_RELAX_detect_steps); tail -40 logs/one_test_${JOB}.log
JOB=$(sbatch --wait --parsable run_all_tests.sh); tail -60 logs/run_all_tests_${JOB}.log
```

Expected: `ONE_TEST_PASS`, and the suite back at the Step 1 count. Report both logs verbatim.

- [ ] **Step 6: Grep for survivors (Agent)**

```bash
grep -rn "unstableFraction\|MaxUnstableFraction\|RejectedByUnstable\|locally unstable" \
  --include=*.m . | grep -v "^./.superpowers" | grep -v "^./variance_ratio_survey"
```

Expected: no hits outside `docs/`. Report the output verbatim.

- [ ] **Step 7: Pin this detector as Task 3's comparison arm (Agent)**

Task 3 measures the new gate against this one, and because this plan commits nothing, the
comparison arm must be snapshotted now: `git show HEAD:RELAX_detect_steps.m` would hand back
the detector as it stood *before* Task 1, still carrying the instability criterion, and the
survey would then be measuring two changes at once.

```bash
mkdir -p slope_gate_survey/legacy
cp RELAX_detect_steps.m slope_gate_survey/legacy/RELAX_detect_steps_level.m
```

Change the function name on its first code line, and only there:

```matlab
function [steps, info] = RELAX_detect_steps_level(EEG, varargin)
```

Nothing else in the file changes: it calls the shared `RELAX_get_segments`,
`RELAX_flat_channels` and `RELAX_rank_pop_channels`, all compatible with Task 1's
two-criterion ranking. The file lives only under `slope_gate_survey/legacy/`, put on the path
by the survey script and never by an entry point, so it cannot shadow the production function.

- [ ] **Step 8: Checkpoint (Evaluator)**

Read the whole diff this task produced:

```bash
git --no-pager diff -- RELAX_detect_steps.m slope_gate_survey/ RELAX_rank_pop_channels.m RELAX_Wrapper_slurm.m RELAX_SET_PARAMETERS_AND_RUN_slurm.m validate_pop_repair.m tests/test_RELAX_detect_steps.m tests/test_RELAX_rank_pop_channels.m tests/test_RELAX_protect_eye_channels.m tests/assert_smoke_outcome.m
git --no-pager status --short -- RELAX_detect_steps.m RELAX_rank_pop_channels.m RELAX_Wrapper_slurm.m RELAX_SET_PARAMETERS_AND_RUN_slurm.m validate_pop_repair.m tests/test_RELAX_detect_steps.m tests/test_RELAX_rank_pop_channels.m tests/test_RELAX_protect_eye_channels.m tests/assert_smoke_outcome.m
```

(`--no-pager` always: a pager called from inside MATLAB's `system()` hangs silently.)
Record in the task notes what changed, the suite count, and anything deferred. **Do not commit.** This plan produces no commits — the whole electrode-pop detection/repair line of work lands as one commit the user makes at the end.

---

### Task 2: The median-slope statistic and gate

**Files:**
- Modify: `RELAX_detect_steps.m` (help block lines 8-28 and 40-43 and 105-113;
  `inputParser` lines 130-151; `info` initialiser lines 157-168; per-channel loop lines
  194-198; `signedRuns` call line 202; `derivedSamples` lines 282-293; `levelStatistic`
  lines 305-318; `gateCandidates` lines 320-334; `refineCandidate` row assembly line 425)
- Test: `tests/test_RELAX_detect_steps.m` (new section after the existing degenerate-channel
  sections; field-name change in the section-1a list)
- Create: `scratch_characterize_gate.m` in the scratchpad directory only — **not** committed

**Interfaces:**
- Consumes: Task 1's detector (no `unstableFraction`, two-criterion ranking).
- Produces:
  - `RELAX_detect_steps(EEG, ..., 'SlopeBaselineMs', 2000, 'SlopeNoiseFloorMicrovolts', 0.1)`
  - `info.channelScaleSlope` (nbchan x 1 double) replacing `info.channelScaleD`
  - `info.params.derived.baseline` (odd sample count)
  - `steps.z` is now the signed local excess `(|Dslope| - b) / max(s, floor)`;
    `steps.channelZ` is `|Dslope| / channelScaleSlope`
  - `levelStatistic(x, W, g, valid, baselineN, noiseFloor) -> [D, Dslope, z]`
  - `gateCandidates(D, Dslope, z, valid, opt) -> [cand, chanZ, scaleSlope]`

- [ ] **Step 1: Characterise both statistics on synthetic shapes (Evaluator, via the agent)**

Before writing a red test that asserts "the old detector finds this and the new one must
not", establish empirically which shapes the **current** detector accepts. Write
`scratch_characterize_gate.m` in the session scratchpad (not the repo):

```matlab
function scratch_characterize_gate()
    srate = 1000; n = 20000; tt = 0:n-1; rng(7);
    mk = @(x) struct('data', x, 'nbchan', 1, 'pnts', numel(x), 'srate', srate, 'trials', 1, ...
        'event', struct('type', {}, 'latency', {}), 'chanlocs', struct('labels', {'E1'}));
    shapes = struct('name', {}, 'x', {});
    base = randn(1, n) * 5;
    shapes(end+1) = struct('name', 'ramp 200uV / 2s',   'x', base + 200 * min(1, max(0, (tt - 8000) / 2000)));
    shapes(end+1) = struct('name', 'ramp 400uV / 4s',   'x', base + 400 * min(1, max(0, (tt - 6000) / 4000)));
    shapes(end+1) = struct('name', 'sine 400uV 0.5Hz',  'x', base + 400 * sin(2*pi*0.5*tt/srate));
    shapes(end+1) = struct('name', 'sine 300uV 0.2Hz',  'x', base + 300 * sin(2*pi*0.2*tt/srate));
    shapes(end+1) = struct('name', 'step 200uV / 1ms',  'x', base + 200 * (tt >= 10000));
    shapes(end+1) = struct('name', 'sigmoid 100uV/50ms','x', base + 100 ./ (1 + exp(-(tt - 10000) / 10)));
    shapes(end+1) = struct('name', 'sigmoid 300uV/400ms','x', base + 300 ./ (1 + exp(-(tt - 10000) / 80)));
    shapes(end+1) = struct('name', 'step 200uV on 0.2Hz drift', ...
        'x', base + 300 * sin(2*pi*0.2*tt/srate) + 200 * (tt >= 10000));
    for k = 1:numel(shapes)
        [s, info] = RELAX_detect_steps(mk(shapes(k).x));
        fprintf('%-28s steps=%d', shapes(k).name, numel(s.sample));
        if ~isempty(s.sample)
            fprintf('  |A|=%.0f rise=%.0fms plateauDrift=%.2f z=%.1f chanZ=%.1f', ...
                max(abs(s.amplitude)), max(s.riseMs), max(s.plateauDrift), max(s.z), max(s.channelZ));
        end
        fprintf('  rejected: %s\n', strjoin(cellfun(@(f) sprintf('%s=%d', f, info.nRejected.(f)), ...
            fieldnames(info.nRejected)', 'UniformOutput', false), ' '));
    end
end
```

Run it with a one-off `sbatch` job. **The evaluator reads this output and records, in the
task notes, which drift shapes the current detector accepts.** Step 2's red test asserts
rejection for those shapes specifically; a shape the current detector already rejects proves
nothing and must not be used as the red test.

- [ ] **Step 2: Write the failing tests (Agent)**

In `tests/test_RELAX_detect_steps.m` section 1a, change the field-name list to the new name:

```matlab
    for f = {'stepsPerChannel','stepRatePerMinute','transitionFraction','channelScaleSlope'}
```

Add a new section after the degenerate-channel sections. **Use the drift shapes Step 1 showed
the current detector accepting**; the shapes below are the expected set and the evaluator
confirms or substitutes them from Step 1's log:

```matlab
    %% ---- 1c. the gate scores how fast the median level moves, not how far ----
    n20 = 20000; t20 = 0:n20-1; base20 = randn(1, n20) * 5;

    % drift the old level gate accepted: a 2 s ramp is not a step
    x = base20 + 200 * min(1, max(0, (t20 - 8000) / 2000));
    steps = RELAX_detect_steps(mkEEG(x));
    assert(isempty(steps.sample), 'a 200 uV ramp over 2 s is drift, %d step(s) found', numel(steps.sample));

    % a slow wave is not a step either
    x = base20 + 400 * sin(2*pi*0.5*t20/srate);
    steps = RELAX_detect_steps(mkEEG(x));
    assert(isempty(steps.sample), 'a 0.5 Hz slow wave is drift, %d step(s) found', numel(steps.sample));

    % a real pop is still found, with its amplitude
    x = base20 + 200 * (t20 >= 10000);
    steps = RELAX_detect_steps(mkEEG(x));
    assert(numel(steps.sample) == 1, 'one step expected, got %d', numel(steps.sample));
    assert(abs(steps.amplitude - 200) < 20, 'amplitude %g should be near 200', steps.amplitude);
    assert(abs(steps.sample - 10000) < 100, 'step located at %d, expected near 10000', steps.sample);

    % a 50 ms sigmoid is a pop
    x = base20 + 100 ./ (1 + exp(-(t20 - 10000) / 10));
    steps = RELAX_detect_steps(mkEEG(x));
    assert(numel(steps.sample) == 1, 'a 50 ms sigmoid is a pop, got %d step(s)', numel(steps.sample));

    % a genuine step riding on drift survives the local baseline subtraction
    x = base20 + 300 * sin(2*pi*0.2*t20/srate) + 200 * (t20 >= 10000);
    steps = RELAX_detect_steps(mkEEG(x));
    assert(numel(steps.sample) == 1, 'a step on drifting baseline must still be found, got %d', ...
        numel(steps.sample));
    assert(abs(steps.sample - 10000) < 150, 'step on drift located at %d', steps.sample);

    % two steps 300 ms apart do not mask each other through the shared baseline span
    x = base20 + 150 * (t20 >= 9000) + 150 * (t20 >= 9300);
    steps = RELAX_detect_steps(mkEEG(x));
    assert(numel(steps.sample) == 2, 'two steps 300 ms apart, got %d', numel(steps.sample));

    % the lag is TransitionMs: a rise longer than it is attenuated by 2g/T and falls under the
    % floor, until TransitionMs is widened to match
    x = base20 + 300 ./ (1 + exp(-(t20 - 10000) / 80));       % ~400 ms rise
    steps = RELAX_detect_steps(mkEEG(x));
    assert(isempty(steps.sample), 'a 400 ms rise is outside a 100 ms lag, got %d', numel(steps.sample));
    steps = RELAX_detect_steps(mkEEG(x), 'TransitionMs', 400, 'MaxTransitionMs', 600);
    assert(numel(steps.sample) == 1, 'widening TransitionMs to 400 must recover it, got %d', ...
        numel(steps.sample));

    % new options are echoed, and the channel scale is the slope statistic's
    [~, info] = RELAX_detect_steps(mkEEG(base20), 'SlopeBaselineMs', 1500);
    assert(info.params.SlopeBaselineMs == 1500 && info.params.SlopeNoiseFloorMicrovolts == 0.1);
    assert(mod(info.params.derived.baseline, 2) == 1, 'the baseline span must be an odd sample count');
    assert(~isfield(info, 'channelScaleD') && isfield(info, 'channelScaleSlope'));
```

- [ ] **Step 3: Run the tests to verify they fail (Agent)**

```bash
JOB=$(sbatch --wait --parsable run_one_test.sh test_RELAX_detect_steps); tail -40 logs/one_test_${JOB}.log
```

Expected: FAIL — the field-list assertion trips first (`info.channelScaleSlope` does not
exist). Report the log verbatim.

- [ ] **Step 4: Implement the statistic and the gate (Agent)**

In `derivedSamples`, add the baseline span (odd, like `W` and `smooth`):

```matlab
    d.baseline = ms2n(opt.SlopeBaselineMs);
    if mod(d.baseline, 2) == 0; d.baseline = d.baseline + 1; end
```

Replace `levelStatistic` with:

```matlab
function [D, Dslope, z] = levelStatistic(x, W, g, valid, baselineN, noiseFloor)
% D(t): the settled level change across a TransitionMs + WindowMs lag, which sets the sign and
% meets the absolute floor. Dslope(t): the same centred median filter differenced across the
% transition allowance alone, so a level change is scored by how fast it happened. z(t): that
% slope in excess of the local baseline slope, in local robust sd - drift raises its own
% baseline and scores nothing. NaN outside valid.
    n = numel(x);
    h = floor(W / 2);
    mPre  = movmedian(x, [W-1 0], 'omitnan');       % median of x(t-W+1 .. t)
    mPost = movmedian(x, [0 W-1], 'omitnan');       % median of x(t .. t+W-1)
    mc    = [mPre(h+1:end) nan(1, h)];              % centred median filter: mc(t) = mPre(t+h)
    D = nan(1, n); Dslope = nan(1, n); z = nan(1, n);
    t = find(valid);
    D(t)      = mPost(t + g) - mPre(t - g);
    Dslope(t) = mc(t + g) - mc(t - g);
    b = movmedian(abs(Dslope), baselineN, 'omitnan');
    s = 1.4826 * movmad(Dslope, baselineN, 'omitnan');
    z(t) = (abs(Dslope(t)) - b(t)) ./ max(s(t), noiseFloor);
end
```

Replace `gateCandidates` with:

```matlab
function [cand, chanZ, scaleSlope] = gateCandidates(D, Dslope, z, valid, opt)
    loc = valid & z > opt.ZThreshold & abs(D) >= opt.MinAmplitudeMicrovolts;
    Dq = Dslope(valid & ~loc);
    Dq = Dq(isfinite(Dq));
    if isempty(Dq)
        scaleSlope = NaN;
    else
        scaleSlope = 1.4826 * median(abs(Dq - median(Dq)));
    end
    if ~(scaleSlope > 0)
        % A channel whose quiet samples have no spread: the channel gate has nothing to say.
        scaleSlope = eps;
    end
    chanZ = abs(Dslope) / scaleSlope;
    cand = loc & chanZ > opt.ChannelZThreshold;
end
```

In the per-channel loop, replace the two calls and the two `info` assignments:

```matlab
        [D, Dslope, z] = levelStatistic(x, d.W, d.g, valid, d.baseline, opt.SlopeNoiseFloorMicrovolts);
        [cand, chanZ, scaleSlope] = gateCandidates(D, Dslope, z, valid, opt);
        info.channelScaleSlope(ch) = scaleSlope;
        if ~any(cand); continue; end
```

(the `info.unstableFraction` line is already gone from Task 1), and pass `Dslope` to the run
finder so runs are signed by the slope:

```matlab
            [rs, re] = signedRuns(cand, Dslope, sgn, d.refr);
```

In the `inputParser`, add the two options next to the existing thresholds:

```matlab
    p.addParameter('SlopeBaselineMs', 2000);
    p.addParameter('SlopeNoiseFloorMicrovolts', 0.1);
```

In the `info` initialiser, rename the field:

```matlab
                   'channelScaleSlope', nan(nb,1), ...
```

In `refineCandidate`'s row assembly, `z` is now a signed excess, so take its maximum
directly:

```matlab
                 'amplitude', A, 'z', max(z(rs:re)), 'channelZ', max(chanZ(rs:re)), ...
```

Rewrite the help block's detector paragraph (lines 8-28) to describe the median-slope gate as
stated under "Design" above, document `SlopeBaselineMs` and `SlopeNoiseFloorMicrovolts` in the
name-value list, change `.channelScaleD` to `.channelScaleSlope` in the `info` description,
restate `.z` and `.channelZ` in the `steps` description, and add to "Known limits" that a rise
longer than `TransitionMs` is attenuated by `2g/T`.

- [ ] **Step 5: Run the task's test and the full suite (Agent)**

```bash
JOB=$(sbatch --wait --parsable run_one_test.sh test_RELAX_detect_steps); tail -40 logs/one_test_${JOB}.log
JOB=$(sbatch --wait --parsable run_all_tests.sh); tail -60 logs/run_all_tests_${JOB}.log
```

Expected: `ONE_TEST_PASS` and the suite at the Task 1 Step 1 count. Report both verbatim.
`tests/test_RELAX_refpop_wiring.m` and `tests/test_RELAX_detect_reference_pops.m` exercise the
same detector through the reference path; if either fails, that is a real finding, not a
fixture to loosen.

- [ ] **Step 6: Re-run the characterisation and compare (Evaluator, via the agent)**

Re-run `scratch_characterize_gate.m`. Expected: every drift shape now reports `steps=0`, every
pop shape still reports one step with an amplitude near its true value. The evaluator reads
the before/after pair and records it in the task notes — this is the first evidence the change
does what it claims.

- [ ] **Step 7: Checkpoint (Evaluator)**

Read the whole diff this task produced:

```bash
git --no-pager diff -- RELAX_detect_steps.m tests/test_RELAX_detect_steps.m
git --no-pager status --short -- RELAX_detect_steps.m tests/test_RELAX_detect_steps.m
```

(`--no-pager` always: a pager called from inside MATLAB's `system()` hangs silently.)
Record in the task notes what changed, the suite count, and anything deferred. **Do not commit.** This plan produces no commits — the whole electrode-pop detection/repair line of work lands as one commit the user makes at the end.

---

### Revision 2026-09-23 — why Tasks 3 and 4 were rewritten

The first Task 3 survey (and the exploratory rounds under `slope_gate_survey/` and
`docs/superpowers/reports/slope_gate_sample*/`) could not say whether most candidates were pops or
drift/other artifacts that later stages handle, for four reasons, each fixed below:

1. **Downstream screens shadowed detection.** Candidates were scored after bracketing,
   `MaxTransitionMs`, plateau drift and the 50 uV floor had already removed some of them, so a
   gate parameter's effect could not be separated from the screens'. Now: detection runs
   **alone** (refinement off, drift screens `Inf`), candidates are scored, and every screen is
   then judged against those scores (Task 4B).
2. **The survey measured a configuration production does not run.** Function defaults, wrapper
   defaults, `validate_pop_repair` and the survey scripts each carried their own values. Now:
   one operating point in all four places (Task 3A).
3. **The lag was far too long.** `TransitionMs` 100 lets 1/f drift build up across the lag.
   Now: the lag is swept over 1/2/4/6/10 samples, including an **uncentred one-sample
   difference** that the symmetric `g` could not express (Task 3A).
4. **Old-vs-new counts are neither precision nor recall.** The legacy level-gate arm and the
   lost/gained test are dropped. Now: candidates from every grid point are pooled into events,
   the **user scores a stratified sample blind**, and each config's precision and relative recall
   are estimated from those shared labels (Tasks 3B–4A).

Decisions (user, 2026-09-23):
- uncentred one-sample difference; `TransitionMs` grid {1, 2, 4, 6, 10};
- ~~the full 4×4 grid of `ZThreshold` × `ChannelZThreshold` over {3, 4, 5, 7}~~ — replaced
  2026-09-24 by `ZThreshold` {7, 12, 20, 35, 50} × `ChannelZThreshold` {10, 35} (see below);
- `SlopeBaselineMs` 1000, with 500 tested;
- refinement off by default;
- survey on the Sep26 files;
- the raw-noise gate is tested in 4B, not added to detection;
- all blind scoring by the user.

**Addendum 2026-09-24 — what the first real-data runs of 3A showed, and what changed.**

1. *The spread collapsed at short lags.* Two neighbouring centred medians are exactly equal
   whenever the entering and leaving samples fall on the same side of the median. So at lags of
   1-4 samples most of `Dslope` is exactly 0, its MAD is 0, and two things fail:
   - the local `z` came to be set by the 0.1 uV floor, turning the gate into an absolute slope
     threshold;
   - the channel scale dropped to about 1e-4, so the ChannelZ gate passed everything.

   On synthetic signals, a 2 Hz wave gave 34-141 false steps (jobs 19890846, 19890951). The fix
   is `SlopeSpreadQuantile` (default 0.9; 0.85 since 2026-09-24, see the note after Task 3A). The spread `s` and the channel scale are now the
   0.9-quantile of `|Dslope - b|`, divided by the Gaussian quantile so they stay in sd units;
   q = 0.5 is exactly the old `1.4826 * MAD`. After the fix, the channel scale stays at 0.04-1.2,
   the noise floor no longer changes any count, and waves pass only at lags 1-2.
2. *A window longer than half a wave's period turns the wave into edges.* A centred median does
   this to waves of 1/(2W)-1/W Hz: 1.7-3.3 Hz (delta) at W 301, and 60 Hz line noise at W 25
   (1.5 cycles). **`WindowMs` is swept over {7, 151, 301}** (user, 2026-09-23; W 7 replaced W 51
   on 2026-09-24 to test what a much shorter window does). W 25 was dropped
   after real data gave it the worst channel rejection and common-mode counts. 51, 151 and 301
   ms all hold whole or near-whole 60 Hz cycles.
3. *The z grid was far too low.* At Z = CZ = 4 every (W, lag) gives 1,800-37,000 steps/min per
   file, 60-95 % of them common-mode (array 19890952, 19892151). Part of the common-mode share
   is dense steps chaining through adjacent 50 ms bins.

   Survival curves (`docs/superpowers/reports/slope_gate_cm_look/survival_all_*.png`) show:
   - on z, a steep fall to about 10-15 and then a shallower tail on m159/m208, while m497/m583
     fall off faster from 20-40;
   - on channel z, a curve that steepens from about 20-40 and never levels off.

   Grid (user, 2026-09-24): **`ZThreshold` {7, 12, 20, 35, 50} × `ChannelZThreshold` {10, 35}**.
4. *The user's diagnosis of the shared events at W 151 / lag 4 / Z = CZ = 4* (renders
   `slope_gate_cm_look/<id>_{cm,shared}NN.png`):
   - Every common-mode event is eye movement, drift or oscillation. **Not one is a real pop.**
   - Of the 20 smaller shared bursts (2 or more channels, not common-mode), only m583 shared03
     and m208 shared05 carry a clear step near the event time, and m208 shared05's step is at
     +25 ms rather than 0.

   Two consequences:
   - Common-mode excision in the repair stage may be removing non-pops. Task 4B re-checks what
     the common-mode flag catches at the chosen operating point.
   - The step location can sit ~25 ms off the edge. This may be the median-filter shoulder
     (a clean pop at W ≥ 151 also yields a secondary candidate 20-45 samples from the edge).
     Shoulder suppression is deferred until realistic thresholds show whether it matters
     (user, 2026-09-24).
5. `info.rawScale` now subtracts each segment's own median. The stitched Sep26 files carry a
   different DC offset per segment, and the whole-channel MAD read about 2 mV.

Survey files (all 1 kHz, 17–21 min, 64–65 segments, median 7.9 s; checked by job 19890611),
in `/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Sep26/rest_treat_joined/`:
`m159_dlpfc_day1_stitched.set`, `m497_dlpfc_day1_stitched.set`, `m208_dlpfc_day1_stitched.set`,
`m583_dlpfc_day1_stitched.set`.

---

### Task 3A: One operating point, an exact lag, and detection-only defaults

Runs at **effort high**: it re-indexes the same shifted median filters Task 2 built, and an
off-by-one there still runs, still finds steps, and is silently wrong.

**Files:**
- Modify: `RELAX_detect_steps.m` (`derivedSamples`, `validMask`, `minSegLen`, `levelStatistic`,
  `gateCandidates` call, `refineCandidate` and `rawCandidate` window positions, `emptySteps`,
  `info` initialiser, `inputParser` defaults, help block)
- Modify: `RELAX_Wrapper_slurm.m` (`PopDetect_*` defaults lines ~339-349; the detector call
  ~442-452 and the reference-pop `'Detect'` list ~401-408)
- Modify: `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` (pop-detector block)
- Modify: `validate_pop_repair.m` (`stageConfig` lines ~166-172 and the call ~186-192)
- Test: `tests/test_RELAX_detect_steps.m` (new sections 1f-1h; a `legacy` option list for
  sections that test bracketing mechanics)
- Test: any other test file that fails under the new defaults — **only** by passing `legacy{:}`,
  see Step 5

**Interfaces:**
- Produces:
  - `info.params.derived.lag` = `max(1, round(TransitionMs * srate / 1000))` exactly;
    `derived.gPre = floor(lag/2)`, `derived.gPost = lag - gPre`; `derived.g` is gone.
  - `Dslope(t) = mc(t + gPost) - mc(t - gPre)`; `D(t) = mPost(t + gPost) - mPre(t - gPre)`. The
    level windows are `x(t-gPre-W+1 .. t-gPre)` and `x(t+gPost .. t+gPost+W-1)`, leaving
    `lag - 1` samples between them. For even lags this is identical to Task 2; for `lag = 1` it
    is the uncentred `mc(t+1) - mc(t)` with adjacent windows.
  - `steps.rawZ` (reported, never gating): `|amplitude| / hypot(robustSd(pre), robustSd(post))`
    over the two raw level windows, the superseded level gate's denominator, so 4B can price a
    raw-noise gate by filtering scored events.
  - `info.floorBoundFraction`, `info.rawScale`, `info.diffScale` (nbchan x 1). In order:
    - the fraction of valid samples where `s < SlopeNoiseFloorMicrovolts`;
    - `1.4826 * MAD` of the raw channel after subtracting each segment's own median;
    - `1.4826 * MAD` of `diff(x)`.

    All are taken over valid samples only.
  - `SlopeSpreadQuantile` (default 0.85, was 0.9): `s(t)` and `channelScaleSlope` are
    `quantile(|v|, q) / (sqrt(2) * erfinv(q))` of `v = Dslope - b` (local, over the
    `SlopeBaselineMs` window, evaluated every `N/20` samples and held between) and of
    `v = Dq - median(Dq)` (channel, quiet samples). `q = 0.5` reproduces `1.4826 * MAD`.
    Implemented in the working tree on 2026-09-23 (`quantileSd`, `movingQuantileSd`, and
    `segmentMedians` for `rawScale`); this task's agent tests and documents it rather than
    re-implementing it.
  - Defaults (function = wrapper = entry script = `validate_pop_repair`): `WindowMs` 300,
    `TransitionMs` 4, `ZThreshold` 12, `ChannelZThreshold` 10, `MinAmplitudeMicrovolts` 15,
    `RefractoryMs` 10, `MaxTransitionMs` 500, `RefineTransition` false,
    `MaxPlateauDriftFraction` Inf, `HoldMs` 500, `MaxHoldDriftFraction` Inf,
    `CommonModeChannelFraction` 0.15, `SlopeBaselineMs` 1000, `SlopeNoiseFloorMicrovolts` 0.1,
    `SlopeSpreadQuantile` 0.85, `StepsPerMinuteToRejectChannel` 6,
    `MaxTransitionFractionToRejectChannel` 0.02. `WindowMs` 300, `TransitionMs` 4, `ZThreshold`
    12 and `ChannelZThreshold` 10 are **provisional** until Task 4C. They were picked from the
    synthetic check in job 19904839. At lag 4, with (12, 10) and with (20, 35), every W in
    {51, 151, 301} gave 0 steps on the 2 Hz and 1 Hz waves and the random walk, and exactly
    1 step on each pop shape and on the step riding drift.
  - New config fields, passed through by the wrapper and `validate_pop_repair`:
    `PopDetect_RefineTransition`, `PopDetect_SlopeBaselineMs`,
    `PopDetect_SlopeNoiseFloorMicrovolts`, `PopDetect_SlopeSpreadQuantile`,
    `PopDetect_MaxHoldDriftFraction`, `PopDetect_HoldMs`.

- [ ] **Step 1: Suite baseline (Agent)**

```bash
mkdir -p logs
JOB=$(sbatch --wait --parsable run_all_tests.sh); tail -60 logs/run_all_tests_${JOB}.log
```

Report the `N/M test files passed` line and the names of any failures. That count, less the known
`test_RELAX_bridged_channels` failure, is the target.

- [ ] **Step 2: Write the failing tests (Agent)**

Near the top of `tests/test_RELAX_detect_steps.m`, after `numericCols`, add `'rawZ'` to
`numericCols` and define the pre-2026-09-23 operating point that the bracketing sections exercise:

```matlab
    % The operating point every section below 1f was written against. Sections that test the
    % bracketing, plateau/hold and floor mechanics pass it explicitly; those mechanics still exist
    % as options, they are just no longer the defaults (2026-09-23).
    legacy = {'WindowMs', 150, 'TransitionMs', 100, 'MinAmplitudeMicrovolts', 50, ...
              'RefractoryMs', 20, 'MaxTransitionMs', 300, 'RefineTransition', true, ...
              'MaxPlateauDriftFraction', 0.3, 'SlopeBaselineMs', 2000, ...
              'CommonModeChannelFraction', 0.08};
```

Add after section 1e:

```matlab
    %% ---- 1f. TransitionMs is Dslope's exact lag, down to an uncentred one-sample difference ----
    rng(7);
    q = randn(1, 20000) * 5;
    for tr = [1 2 4 6 10 100]
        [~, info] = RELAX_detect_steps(mkEEG(q), 'TransitionMs', tr);
        d = info.params.derived;
        assert(d.lag == tr, 'TransitionMs %d at 1 kHz must be a %d-sample lag, got %d', tr, tr, d.lag);
        assert(d.gPre == floor(tr/2) && d.gPost == tr - floor(tr/2), ...
            'lag %d must split floor/ceil, got %d + %d', tr, d.gPre, d.gPost);
    end
    assert(~isfield(info.params.derived, 'g'), 'the symmetric half-gap g was replaced by gPre/gPost');

    % an ideal step is found once at every lag in the grid, at its edge, at full amplitude, and
    % with refinement off its interval is the gate's own run: it contains the edge and stays local
    x = q + 200 * ((0:19999) >= 10000);
    for tr = [1 2 4 6 10]
        steps = RELAX_detect_steps(mkEEG(x), 'TransitionMs', tr);
        assert(numel(steps.sample) == 1, 'lag %d: one step expected, got %d', tr, numel(steps.sample));
        assert(abs(steps.sample - 10000) <= tr, 'lag %d: step at %d, expected within %d of 10000', ...
            tr, steps.sample, tr);
        assert(abs(steps.amplitude - 200) < 10, 'lag %d: amplitude %g, expected ~200', tr, steps.amplitude);
        assert(steps.onset <= 10000 && steps.offset >= 9999, ...
            'lag %d: interval [%d %d] must contain the edge', tr, steps.onset, steps.offset);
        assert(steps.offset - steps.onset < 50, ...
            'lag %d: the gate run must stay local, got %d samples', tr, steps.offset - steps.onset);
        assert(abs(steps.rawZ - 200 / (5 * sqrt(2))) < 5, ...
            'lag %d: rawZ %g, expected ~%g', tr, steps.rawZ, 200 / (5 * sqrt(2)));
    end

    % the lag-1 windows are indexed from the segment start: a step 600 samples in (two windows)
    % is found at its edge, not shifted and not an index error. Z 20, not the provisional 12: at
    % W 301 / lag 1 a median-bend shoulder 12 samples past this edge reaches z 14 for this noise
    % draw (1 of 40 seeds; none at Z 20); how often shoulders occur is measured by the survey.
    x = q; x(1:600) = x(1:600) + 100;              % a step 600 samples in: inside the valid span
    steps = RELAX_detect_steps(mkEEG(x), 'TransitionMs', 1, 'ZThreshold', 20);
    assert(numel(steps.sample) == 1 && abs(steps.sample - 600) <= 1, ...
        'lag 1: a step 600 samples in must be found at its edge');

    %% ---- 1g. diagnostics for the noise floor and the channel's raw noise ----
    [~, iq]  = RELAX_detect_steps(mkEEG(q), 'TransitionMs', 1);
    [~, iq5] = RELAX_detect_steps(mkEEG(q * 100), 'TransitionMs', 1);
    for f = {'floorBoundFraction', 'rawScale', 'diffScale'}
        assert(isequal(size(iq.(f{1})), [1 1]), 'info.%s must be nbchan x 1', f{1});
    end
    assert(iq.floorBoundFraction >= 0 && iq.floorBoundFraction <= 1);
    assert(iq5.floorBoundFraction < iq.floorBoundFraction, ...
        'a fixed floor binds less on a noisier channel: %g vs %g', iq5.floorBoundFraction, iq.floorBoundFraction);
    assert(abs(iq.rawScale - 5) < 0.5, 'rawScale %g of white noise sd 5', iq.rawScale);
    assert(abs(iq.diffScale - 5 * sqrt(2)) < 0.7, 'diffScale %g of white noise sd 5', iq.diffScale);

    %% ---- 1h. the defaults are the detection-only operating point of the entry script ----
    [~, info] = RELAX_detect_steps(mkEEG(q));
    P = info.params;
    assert(isequal(P.RefineTransition, false), 'refinement is off by default (2026-09-23)');
    assert(P.WindowMs == 300 && P.TransitionMs == 4 && P.ZThreshold == 12 && ...
           P.ChannelZThreshold == 10 && P.MinAmplitudeMicrovolts == 15 && P.RefractoryMs == 10 && ...
           P.MaxTransitionMs == 500 && isinf(P.MaxPlateauDriftFraction) && P.HoldMs == 500 && ...
           isinf(P.MaxHoldDriftFraction) && P.CommonModeChannelFraction == 0.15 && ...
           P.SlopeBaselineMs == 1000 && P.SlopeNoiseFloorMicrovolts == 0.1 && ...
           P.SlopeSpreadQuantile == 0.85 && P.CommonModeGapMs == 10 && ...
           ~isfield(P, 'CommonModeWindowMs'), ...
        'function defaults must match RELAX_SET_PARAMETERS_AND_RUN_slurm.m');

    %% ---- 1i. the spread is a high quantile, so it survives a mostly-zero short-lag Dslope ----
    % Neighbouring centred medians are equal whenever the entering and leaving samples fall on the
    % same side of the median, so at lag 1 most of Dslope is exactly 0 and its MAD collapses.
    [~, iq9]   = RELAX_detect_steps(mkEEG(q), 'TransitionMs', 1);
    [~, iqMad] = RELAX_detect_steps(mkEEG(q), 'TransitionMs', 1, 'SlopeSpreadQuantile', 0.5);
    assert(iq9.channelScaleSlope > 0.01, ...
        'the default-q channel scale must not collapse at lag 1, got %g', iq9.channelScaleSlope);
    assert(iqMad.channelScaleSlope < iq9.channelScaleSlope / 10, ...
        'q 0.5 is the MAD, which does collapse at lag 1: %g vs %g', ...
        iqMad.channelScaleSlope, iq9.channelScaleSlope);

    % a DC offset between inter-boundary segments is not channel noise
    x = q; x(10001:end) = x(10001:end) + 1000;
    [~, iseg] = RELAX_detect_steps(mkEEG(x, 10000.5), 'TransitionMs', 1);
    assert(abs(iseg.rawScale - 5) < 0.5, ...
        'rawScale %g must measure noise within segments, not the 1000 uV offset between them', ...
        iseg.rawScale);
```

- [ ] **Step 3: Run it to verify it fails (Agent)**

```bash
JOB=$(sbatch --wait --parsable run_one_test.sh test_RELAX_detect_steps); tail -40 logs/one_test_${JOB}.log
```

Expected: FAIL in 1f on `derived.lag` (the field does not exist). Report the log verbatim.

- [ ] **Step 4: Implement (Agent)**

`derivedSamples` — replace the `d.g` line:

```matlab
    d.lag   = ms2n(opt.TransitionMs);      % Dslope's lag in samples, exactly
    d.gPre  = floor(d.lag / 2);
    d.gPost = d.lag - d.gPre;              % a one-sample lag reads mc(t+1) - mc(t)
```

`minSegLen = 2*d.W + d.lag - 1;` and `valid = validMask(segments, EEG.pnts, d.W, d.gPre, d.gPost);`

```matlab
function valid = validMask(segments, n, W, gPre, gPost)
% true where both level windows fit inside the segment holding t
    valid = false(1, n);
    for s = 1:size(segments,1)
        lo = segments(s,1) + gPre + W - 1;
        hi = segments(s,2) - gPost - W + 1;
        if hi >= lo; valid(lo:hi) = true; end
    end
end
```

`levelStatistic(x, W, gPre, gPost, valid, baselineN, noiseFloor) -> [D, Dslope, z, s]`, with
`D(t) = mPost(t + gPost) - mPre(t - gPre)` and `Dslope(t) = mc(t + gPost) - mc(t - gPre)`;
everything else in it is unchanged. Returning `s` lets the channel loop fill
`info.floorBoundFraction(ch) = mean(s(valid) < opt.SlopeNoiseFloorMicrovolts)`. In the same loop,
fill `info.rawScale(ch)` and `info.diffScale(ch)` with `robustSd` of `x(valid)` and of
`diff(x)` restricted to pairs whose samples are both valid.

Window positions, everywhere `d.g` appeared: the pre-window ends at `t - d.gPre` and the
post-window starts at `t + d.gPost`, so in `refineCandidate` `searchLo = max(a, rs - (d.gPre + d.h))`,
`searchHi = min(b, re + (d.gPost + d.h))`, the two level windows at lines ~428-429 and ~452-453,
and in `rawCandidate` `preWin`/`postWin` follow that rule. Both `rawCandidate` and
`refineCandidate` add `'rawZ', abs(A) / max(hypot(robustSd(pre), robustSd(post)), eps)` to the
row; add `rawZ` to `emptySteps` and to the column assembly after `varianceRatio`.

Defaults: set the `inputParser` defaults, the wrapper's `PopDetect_*` defaults, the entry
script's values and `validate_pop_repair`'s `stageConfig` to the list under Interfaces. Add the
five new config fields to the wrapper default struct, to both wrapper detector calls, to the
entry script (next to their siblings, one-line comments in the existing style) and to
`validate_pop_repair`. In the entry script, `PopDetect_TransitionMs` goes from 100 to 4 and
`PopDetect_WindowMs` stays 300.

Help block: `TransitionMs` becomes "the exact lag of Dslope in samples (min 1), split
floor/ceil either side of t; 1 is the uncentred difference `mc(t+1) - mc(t)` with adjacent level
windows"; the rise-attenuation note becomes `L/T` (`L` = lag) in place of `2g/T`; the
placement-margin paragraph becomes `W + gPre` / `W + gPost`; document `rawZ` and the three new
`info` fields; `RefineTransition`'s paragraph states the default is false.

- [ ] **Step 5: Run the suite; pin, never loosen (Agent, Evaluator rules)**

```bash
JOB=$(sbatch --wait --parsable run_one_test.sh test_RELAX_detect_steps); tail -60 logs/one_test_${JOB}.log
JOB=$(sbatch --wait --parsable run_all_tests.sh); tail -60 logs/run_all_tests_${JOB}.log
```

Changing the defaults will break sections written against the old operating point. For each
failure, the agent reports it and **stops**; the evaluator rules on it:
- **Pin** — the section tests a mechanism that still exists as an option (bracketing,
  plateau/hold rejection, the 50 uV floor, `MaxTransitionMs`, event grouping at 0.08, run merging
  at 20 ms). The only permitted edit is appending `legacy{:}` (or, in another test file, an
  equivalent option list defined in that file) to the detector call. Assertions do not change.
- **Finding** — the section asserts something the new defaults are meant to keep true (section
  1c's drift/pop claims, 2a's detection and location, boundaries in section 4). Those run at the
  new defaults, and a failure there is reported to the user, not pinned.

Repeat until `ONE_TEST_PASS` and the suite is at the Step 1 count. Report every pinned call site.

- [ ] **Step 6: Four-way default check (Agent)**

```bash
for k in WindowMs TransitionMs ZThreshold ChannelZThreshold MinAmplitudeMicrovolts RefractoryMs \
         MaxTransitionMs RefineTransition MaxPlateauDriftFraction HoldMs MaxHoldDriftFraction \
         CommonModeChannelFraction SlopeBaselineMs SlopeNoiseFloorMicrovolts; do
  echo "== $k"; grep -n "'$k'\|PopDetect_$k\b" RELAX_detect_steps.m RELAX_Wrapper_slurm.m \
       RELAX_SET_PARAMETERS_AND_RUN_slurm.m validate_pop_repair.m | tr -d '\r'
done
grep -rn "derived\.g\b\|d\.g\b" --include=*.m . | grep -v "^./.superpowers\|^./slope_gate_survey\|^./variance_ratio_survey"
```

Report verbatim. Expected: one value per parameter across all four files, and no `d.g` outside
the superseded exploratory scripts under `slope_gate_survey/` (see Task 3B).

- [ ] **Step 7: Checkpoint (Evaluator)**

Diff every file listed above. Diff each test file against this task's text: the only changes
to existing sections are `legacy{:}` appends at call sites the evaluator ruled "pin". Record
the pinned sites, the suite count and any findings. **Do not commit.**

### Decisions and findings, 2026-09-24 (user rulings)

- **`SlopeSpreadQuantile` 0.9 -> 0.85** (function, wrapper, entry script, `validate_pop_repair`).
  At q 0.9 a linear 100-150 ms ramp masks itself: the ramp raises many `Dslope` samples at once,
  the 0.9-quantile spread absorbs them, and the local z never clears the gate. Legacy-operating-
  point runs of `test_RELAX_detect_steps` sections 2d-3 differing only in q: q 0.9 found no step
  on the 2d, 2e and section-3 ramp fixtures; q 0.5 passed them all (scratch jobs 19908068,
  19908072). q 0.9 had been chosen for shoulder suppression (1/40 seeds at W 301 / lag 1 / 12/10,
  none at Z >= 20); shoulders were not re-measured at 0.85.
- **Section 3 of `test_RELAX_detect_steps` passes `legacy{:}`**: it tests the bracketing that
  the defaults switch off. **Open decision:** once the defaults
  are settled (Task 4C), decide whether that machinery is removed from `RELAX_detect_steps`
  entirely, and those sections with it.
- **Section 4's "step well inside a segment" fixture** moved from 160 to 350 samples in: a step
  must sit about W samples in, and W is 300 at the defaults.
- Section 3's drift channel recovered only 4 of 6 pops at the legacy operating point with q 0.5
  (the test wants >= 5); at q 0.85 it passes. The expectation stands.
- **Findings at q 0.85 (user kept 0.85):** ramps that fill much of the local window still mask
  themselves. The 2d fixture (300 uV over 300 ms, then a jump) and the 2e 100 ms rise give no
  candidate at the legacy operating point (both found at q 0.5), and section 6's 16 x 150 ms ramps
  give transition fraction 0.8% (< 2%) at the defaults. Section 3's 150 ms ramp is found.
- **Sections 2d and 2e parked** (user): the plateau- and hold-drift gates are not tested until it
  is decided whether they stay; the sections are in `tests/parked/`. **Section 6's transition-
  fraction case** now uses transitions no longer than 15 ms (90 alternating 15 ms ramps/min).
- **`test_RELAX_detect_reference_pops` and `test_RELAX_bridged_channels` parked** in
  `tests/parked/` (user): those features are in limbo; restore the tests when work resumes.
- **Event grouping** (`CommonModeGapMs` 10, replacing 50 ms `CommonModeWindowMs` bins): steps
  chain into one event while each is within 10 ms of the previous one, across channels. Fixed bins
  joined steps up to ~100 ms apart depending on bin phase. 10 ms covers the worst placement spread
  of one edge across channels at every grid lag (L - 1 samples, L <= 10).
- **Common-mode events are flagged, not excised, by default** (`PopRepair_ExciseCommonMode`
  false): they are neither cut, filled nor partition points; pass 1 still detrends every channel.
- **Z 20 pins in `test_RELAX_detect_steps`** (user): section 3's drift channel (whose 150 ms ramp
  was removed; >= 4 of 5 pops), section 5, both 5b calls and section 6's transition-fraction case.
  At q 0.85 / Z 12 a 300 uV edge can carry a median-bend shoulder at z 12-14 13-20 samples away
  (4 among ~45 large edges in sections 5/5b; 93 vs 90 steps in section 6), and 1/f drift gives
  z 13-18 false steps (job 19909316). **Finding for 4A**: shoulders are back at the provisional
  default Z 12. `validate_pop_repair.m` now also has `PopRepair_ExciseCommonMode = false`.
- **Task 3A checkpoint:** full suite 41/41 (job 19909436; 43 files less the two parked);
  four-way default check clean after the `validate_pop_repair` fix.
- **Grid `WindowMs` {7, 151, 301}**: W 7 replaces W 51. Partial W 7 check (job 19905531, cancelled):
  on synthetic signals W 7 gave no false steps from noise, 2/10/100 Hz waves or the walk and found
  the 200 uV pop, but missed the 30 uV step and the slow sigmoid, and caught the faster sigmoid
  and exponential tail only at some settings. On real data at 12/10 lag 1 it gave ~200 (m583)
  and ~520 (m159) steps/min.
- **`'ThresholdPairs'`** (N x 2 `[Z CZ]`): the per-channel statistics are computed once and gated
  for each pair; each pair's result is identical to a single-pair call. Channels run in a `parfor`
  when a pool is already open (the detector never opens one). m159: 10-pair grid 231 s serial,
  85 s on 8 workers, against ~116-190 s per single call before. Task 3B's `detect_grid` should
  call it once per (W, lag) with all ten pairs, on a pool.
- **Correction to the survival curves** (`zhist_summarise`): they ran each (W, lag) at Z = CZ = 4
  and filtered by threshold afterwards. A loose run merges into one run what splits into several
  at stricter thresholds, so the curves under-count candidates at high z and flatten there. Grid
  choices drawn from their high-z tail are weaker than stated; Task 3B detects at each pair and is
  unaffected.

---

### Task 3B: Detection grid, candidate pool, blind scoring sample

**Files:**
- Create: `slope_gate_survey/detect_grid.m` (one array task: one file × one W × one lag × one baseline,
  all 10 z pairs in one `'ThresholdPairs'` call)
- Modify: `RELAX_detect_steps.m` (adds the `info.zeroSlopeFraction` diagnostic only)
- Create: `slope_gate_survey/pool_candidates.m`, `slope_gate_survey/score_configs.m`
- Create: `slope_gate_survey/draw_scoring_sample.m` (pools, stratifies, samples, renders)
- Create: `run_slope_gate_grid.sh` (array `1-120`), `run_slope_gate_pool.sh`
- Test: `tests/test_slope_gate_pool.m`
- Superseded, left untouched and not run: `slope_gate_survey/survey_slope_gate.m`,
  `sample_blind_steps.m`, `rerender_blind_steps.m`, `sample_config_candidates.m`,
  `sweep_step_configs.m`, `compare_bracket_modes.m`, `legacy/`. Some read `derived.g`, which
  Task 3A removed. `matchSteps.m` and its test stay; the suite still runs them. Whether to delete
  or park the superseded files is decided at Task 6's hand-back.

**Interfaces:**
- `pool_candidates(stepsByConfig, tolSamples) -> [ev, member, count]`. `stepsByConfig` is a
  cell of detector `steps` tables for one file. Common-mode steps are dropped. On each channel,
  steps whose `[onset - tol, offset + tol]` intervals overlap, transitively, form one event.
  `ev` is a table (`channel`, `label`, `onset`, `offset`, `sample` = median member sample,
  `nMembers`); `member(e, c)` is true when config `c` put at least one step in event `e`;
  `count(e, c)` is how many.
- `score_configs(w, score, member) -> T`, with `T` one row per config: `precision`, `recall`,
  `nScored`, `unclearShare`. `score` is a char array of `'p'` (pop, to be repaired here), `'d'`
  (handled downstream: drift, slow wave, blink, muscle, other), `'u'` (unclear). `w` holds the
  sampling weights. Precision is `sum(w.*isP.*m) / sum(w.*(isP|isD).*m)`. Recall is
  `sum(w.*isP.*m) / sum(w.*isP)`, **relative to the grid's union**, not absolute. Unclear events
  are left out of both and reported as `unclearShare`.

- [ ] **Step 1: Write the failing test (Agent)**

`tests/test_slope_gate_pool.m`:

```matlab
function test_slope_gate_pool()
% pool_candidates merges one step found by several configurations into one event and keeps
% other channels, distant steps and common-mode steps out of it; score_configs is the weighted
% precision / relative-recall estimator the scoring sample is read with.
    rng(7);
    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'slope_gate_survey'));
    mk = @(ch, on, off, cm) struct('channel', ch(:), ...
        'label', {arrayfun(@(k) sprintf('E%d', k), ch(:), 'UniformOutput', false)}, ...
        'onset', on(:), 'offset', off(:), 'sample', round((on(:) + off(:)) / 2), ...
        'commonMode', logical(cm(:)));

    A = mk([1 1 2], [1000 5000 1000], [1004 5002 1003], [0 0 0]);
    B = mk([1 1],   [1008 7000],      [1030 7010],      [0 1]);    % 7000 is common-mode
    C = mk([1 1 1], [1001 1012 5003], [1002 1020 5004], [0 0 0]);  % a staircase: two steps, one event
    [ev, member, count] = pool_candidates({A, B, C}, 5);
    assert(height(ev) == 3, 'events E1@1000, E1@5000, E2@1000 expected, got %d', height(ev));
    e1 = find(ev.channel == 1 & ev.onset == 1000);
    assert(isscalar(e1) && ev.offset(e1) == 1030, 'an event spans the union of its members');
    assert(isequal(member(e1, :), [true true true]));
    assert(isequal(count(e1, :), [1 1 2]), 'the staircase config contributes two steps to one event');
    e2 = find(ev.channel == 1 & ev.onset == 5000);
    assert(isequal(member(e2, :), [true false true]));
    e3 = find(ev.channel == 2);
    assert(isequal(member(e3, :), [true false false]), 'the same samples on another channel are another event');
    assert(~any(ev.onset == 7000), 'common-mode steps are not pooled');

    % the tolerance is exact: padded intervals that share a sample join, one sample apart do not
    ev2 = pool_candidates({mk(1, 100, 110, 0), mk(1, 120, 130, 0)}, 5);
    assert(height(ev2) == 1, 'gap of 2*tol joins');
    ev3 = pool_candidates({mk(1, 100, 110, 0), mk(1, 121, 131, 0)}, 5);
    assert(height(ev3) == 2, 'gap of 2*tol + 1 does not');

    % an empty config contributes nothing and keeps its column
    [ev4, m4] = pool_candidates({A, mk([], [], [], [])}, 5);
    assert(height(ev4) == 3 && size(m4, 2) == 2 && ~any(m4(:, 2)));

    % estimator: config 1 holds events 1-3, config 2 holds events 3-4
    w      = [1; 1; 2; 4];
    score  = ['p'; 'd'; 'p'; 'u'];
    member = logical([1 0; 1 0; 1 1; 0 1]);
    T = score_configs(w, score, member);
    assert(abs(T.precision(1) - 3/4) < 1e-12, 'config 1: (1+2)/(1+1+2)');
    assert(abs(T.precision(2) - 1)   < 1e-12, 'config 2: the unclear event leaves the denominator');
    assert(abs(T.recall(1) - 1)      < 1e-12, 'config 1 holds every scored pop');
    assert(abs(T.recall(2) - 2/3)    < 1e-12, 'config 2 holds weight 2 of 3');
    assert(abs(T.unclearShare(2) - 4/6) < 1e-12);
    assert(isequal(T.nScored(:)', [3 2]));

    fprintf('test_slope_gate_pool: PASS\n');
end
```

- [ ] **Step 2: Run it to verify it fails (Agent)** — `run_one_test.sh test_slope_gate_pool`;
  expected `Undefined function 'pool_candidates'`. Report verbatim.

- [ ] **Step 3: Implement `pool_candidates` and `score_configs` (Agent)**, then re-run the test to
  `ONE_TEST_PASS`. `pool_candidates` sorts one channel's padded intervals by start and sweeps,
  opening a new event whenever `start > runningEnd`. Configs with no steps still get their
  `member`/`count` column.

- [ ] **Step 4: `detect_grid` (Agent)**

`detect_grid(taskId, outDir)`: `taskId` 1..120 decodes to file (4) × `WindowMs` {7,151,301} ×
lag {1,2,4,6,10} × baseline {1000,500} (file slowest). For that file it:
1. Loads the file with `pop_loadset`.
2. Restricts it to `/home/nga4004/TMS_response_variability/montage_204_egi_indices.txt` plus the
   `E257` drop, as `survey_slope_gate.m`'s `restrictMontage` does.
3. Applies the impedance screen exactly as stage 2 of `RELAX_Wrapper_slurm.m` does
   (`RELAX_impedance_bad_channels(EEG.chanlocs, 500, 0.25, ...)`, then `pop_select`), so no
   candidate lands on a channel production has already removed. The bridge screen is skipped
   (it needs `ebridge`); say so in the log.
4. Opens a local parpool of `SLURM_CPUS_PER_TASK` workers with its own `JobStorageLocation`
   under `$TMPDIR`, exactly as `RELAX_Wrapper_slurm.m` does (the detector never opens one).
5. Calls `RELAX_detect_steps` **once**, with `'ThresholdPairs'` holding all ten
   `(ZThreshold, ChannelZThreshold)` rows of `{7,12,20,35,50} × {10,35}` (Z slowest), and every
   other option at the Task 3A defaults, passed **explicitly** from one `base` cell printed at
   the top of the log (`WindowMs`, `TransitionMs`, `SlopeBaselineMs` from the task id).
   Element k of the returned `steps`/`info` arrays is identical to a single-pair call
   (Task 3A), so nothing is lost by the single call.

**Revised 2026-09-24** (steps 4-5 above replaced ten single-pair calls): one call per
(file, W, lag, baseline) on a pool is ~3x faster than ten serial calls (m159: 85 s on 8 workers).

It saves `<outDir>/grid_<fileId>_W<W>_L<lag>_B<baseline>.mat` holding:
- `cfg`: a 10-row table of `W`, `Z`, `CZ`, `lag`, `baseline`;
- `steps`: a 10-cell array of detector `steps` structs (`num2cell` of the returned array);
- `summ`: a 10-row table with `nSteps`, `nNonCM`, `nCMEvents`, `nChannelsBad`, `badChannels`,
  `badReasons`, `medianRiseMs`, `secondsDetected`, `elapsedSeconds` (the one call's wall time,
  repeated);
- `chan`: `labels`, `rawScale`, `diffScale`, `floorBoundFraction`, `zeroSlopeFraction` (these do
  not depend on the thresholds; take them from `info(1)`) and `channelScaleSlope` as an
  nChannels × 10 matrix (it depends on the thresholds through its quiet mask);
- the impedance-rejected labels.

It prints one line per z pair, and a floor line: the median and 90th percentile of
`floorBoundFraction` over channels, the number of channels whose `channelScaleSlope` fell back
to `eps` (at any pair), and the median over channels of `zeroSlopeFraction`.

**Detector diagnostic (Agent, part of this step):** add `info.zeroSlopeFraction` (nbchan × 1,
NaN for degenerate channels, like `floorBoundFraction`): the fraction of valid samples where
`Dslope` is exactly 0. It is computed in the statistics stage next to `floorBoundFraction`,
copied into every pair's `info`, initialised in `makeInfo`, and listed in the help text's
`info` block. It changes no detection result. `test_RELAX_detect_steps` and the full suite must
still pass unchanged.

`run_slope_gate_grid.sh`: copy `run_survey_reference_pops.sh`; `--array=1-120%40`, `--mem=32G`,
`--time=03:00:00`, `--cpus-per-task=8`, log `logs/slope_gate_grid_%A_%a.log`, call
`detect_grid(${SLURM_ARRAY_TASK_ID}, 'slope_gate_survey/grid')`.

- [ ] **Step 5: `draw_scoring_sample` (Agent)**

`draw_scoring_sample(gridDir, figDir)`:

1. **Pool.** Per file, pool all 300 configs (3 W × 5 lags × 2 baselines × 10 z pairs) with
   `pool_candidates(..., tol = 10 samples)`. Record, per event:
   - `strict`: the largest `ZThreshold` of any member config, one of {7,12,20,35,50};
   - `fastLag`: whether any member config has lag ≤ 2;
   - `chained`: whether the event is longer than 1000 ms. Chained events are counted and
     reported, and not sampled.
2. **Strata.** `strict` (5 levels) × `fastLag` (2 levels) gives 10 strata, with `k = min(20, U)`
   units drawn from each.
3. **Sample, two stages per stratum.** A unit is a (file, channel) pair with at least one event
   in the stratum. Draw `k` units uniformly without replacement (`U` = units in the
   stratum), then one event uniformly within each drawn unit. An event's weight is
   `w = (U / k) * nEventsOnUnitInStratum`. This stops one telegraph channel filling the sample
   while keeping the estimate unbiased. Seed `rng(20260923)`.
4. **Hidden repeats.** Add 20 of the drawn events a second time under new numbers, to measure how
   consistent the scorer is with themself.
5. **Shuffle and number.** Shuffle everything and number it `cand001..candNNN`.
6. **Render** each candidate with the three-panel layout of `rerender_blind_steps.m`'s plotting
   code (copy it; don't call the superseded file):
   - 10 s of context;
   - a 300 ms zoom on the event;
   - the 5 nearest neighbours over ±1 s.

   Shade the event's pooled `[onset, offset]`. Title: `candNNN   <file id>   <label>   <sample>`.
   **Nothing** on the figure or in its filename says which config, lag, z, stratum or weight
   produced it.
7. **Write:**
   - `<figDir>/scores.csv` with columns `png,score,notes`, the last two left empty for the user;
   - `<figDir>/key.mat`: per candidate the event, file, stratum, weight, member row, count row,
     repeat-of, and the event's `amplitude`, `riseMs`, `plateauDrift`, `holdDrift`, `rawZ`, `z`
     and `channelZ` taken from its loosest member config;
   - `<figDir>/grid_summary.csv`: every file × config row of `summ`;
   - `<figDir>/pool_summary.csv`: events per stratum per file, and the chained count.
8. **Heatmaps.** Render `sweep_L<lag>_B<baseline>.png`: nSteps per minute (non-common-mode) over
   the Z × CZ grid, one panel per file. This is the unscored view of how the grid behaves.

`figDir` = `docs/superpowers/reports/slope_gate_pool/`. `run_slope_gate_pool.sh`: `--mem=48G`,
`--time=04:00:00`, submitted with `--dependency=afterok:<array job id>`.

- [ ] **Step 6: Run (Agent submits, Evaluator reads)**

```bash
mkdir -p logs slope_gate_survey/grid
AJ=$(sbatch --parsable run_slope_gate_grid.sh)
PJ=$(sbatch --parsable --dependency=afterok:${AJ} run_slope_gate_pool.sh)
echo "grid ${AJ}  pool ${PJ}"
```

The evaluator reads:
- all 120 grid logs' floor lines and z-pair lines;
- the pool log;
- `pool_summary.csv`;
- the heatmaps;
- **four rendered candidates, for rendering quality only**, without forming or recording a
  score.

Three conditions are findings for the user **before** scoring starts:
- If `floorBoundFraction` has a median above 0.5 at any lag, the floor, not the local spread,
  is setting `z` at that lag, and the user decides whether to continue.
- If channels fall back to `channelScaleSlope = eps` at any lag, the channel gate is vacuous
  there: every candidate passes it, and the CZ axis of the grid means nothing at that lag. This
  was expected at lag 1 before `SlopeSpreadQuantile` (a q-quantile, not the MAD, now sets the
  channel scale), so any fallback now is a finding. `zeroSlopeFraction` shows how much of
  `Dslope` is exactly 0 at each lag.
- If any stratum has fewer than 10 units, the user decides whether to merge strata.

- [ ] **Step 7: Checkpoint (Evaluator)** — diff the new files and the test; suite count; record the
  job ids and log paths. **Do not commit.**

**3B run record (2026-09-24).** Steps 1-5: Red job 19909560, green 19909616; suite 42/42 (job
19911455). Grid array 19911744 (120/120 COMPLETED, 90-160 s per detector call on 8 workers),
logs `logs/slope_gate_grid_19911744_<n>.log`. Pool job 19911989 (9 min), log
`logs/slope_gate_pool_19911989.log`: 118,673 / 86,895 / 190,602 / 166,834 pooled events
(m159 / m208 / m497 / m583), 0 chained, 220 candidates (200 + 20 hidden repeats) in
`docs/superpowers/reports/slope_gate_pool/`. Profiling (m159, 300 configs): 232 s, 64 % of it
figure printing; `pool_candidates` 16 s. Not optimised. After the suite run, the evaluator added a
dashed marker at the candidate's time on every panel and relabelled the zoom panel `±300 ms`.

Pre-scoring checks:
- **Floor sets z at W 301 / lag 1:** median `floorBoundFraction` 0.78-0.94 in all four files.
  It is lower but material at W 151 / lag 1 (0.14-0.41) and W 301 / lag 2 (0.10-0.48). **User
  (2026-09-24): continue, and keep W 301 / lag 1 in the comparison.** The decision rule and its
  tie-breakers are in `docs/superpowers/reports/2026-09-23-slope-gate-survey.md`.
- No channel fell back to `channelScaleSlope = eps` anywhere.
- No stratum has fewer than 10 units.
- Counts stay high at the strictest pair: at W >= 151, Z 50 still gives 17-320 non-common-mode
  steps/min.
- Stricter CZ can *raise* the step count: at m583 W 151 / lag 1 / Z 7, CZ 35 gives 1016/min
  against 935/min at CZ 10. Fewer candidate samples split one run into several steps; this is the
  same effect as the survival-curve correction above.

---

### Task 3C: Blind scoring (User)

- [x] **Step 1: Write down the decision rule before scoring (User).** The rule goes at the top of
  `docs/superpowers/reports/2026-09-23-slope-gate-survey.md`, under "Decision rule (written
  before scoring)". It says how the operating point will be picked from precision and recall.
  A suggested form: "the config with the highest relative recall among those whose lower 95 %
  bound on precision is at least X; ties go to the larger lag, then the higher thresholds". The
  user chooses X and the form.
- [x] **Step 2: Score (User).** Open `docs/superpowers/reports/slope_gate_pool/cand*.png` in
  order and fill `scores.csv`: `p` = pop to repair here, `d` = something a later stage handles
  (drift, slow wave, blink, muscle, other), `u` = unclear, plus optional notes. Do not open
  `key.mat` or the grid/pool summaries until scoring is finished.
- [x] **Step 3: Validate the sheet (Evaluator).** Every row carries one of `p`/`d`/`u`. Then
  unblind. (2026-09-24: 220/220 scored, 49 p, 110 d, 61 u, 4 with notes.)

---

### Task 3D: Redraw and rescore (amendment, 2026-09-24, after the first round was unblinded)

**Why.** The first round (3B/3C, 200 first scores) found no config eligible under the
pre-registered rule. That first pass led to changes in how detections are counted (report,
"Estimation definitions (decided after scoring)"). Three problems made the sample unfit to
re-estimate from:
- the pool had dropped common-mode steps, so about two thirds of the real detections were
  invisible;
- 130 of the 200 draws were eye events under the new eye rule (version B), so they carry no
  information;
- the weight (U/k) x nEventsOnUnit let one channel's pop (`cand016`) carry 62 % of the weighted
  pop mass.

The user also asked to widen the grid to Z {5, 7, 12, 20, 35, 50} x CZ {5, 10, 35}. The first
round's scores stay in the report as a pilot and are **not** pooled with the second round.

**Files:**
- Modify: `slope_gate_survey/detect_grid.m`, `run_slope_gate_grid.sh` (done: ZLIST/CZLIST gained
  5; output dir from `GRID_DIR`; run 19918057 writes `slope_gate_survey/grid_ext/`)
- Modify: `slope_gate_survey/pool_candidates.m`: an opt-in name-value `'KeepCommonMode'`
  (default false, so `tests/test_slope_gate_pool.m` is unchanged)
- Create: `slope_gate_survey/draw_scoring_sample_v2.m`, `run_slope_gate_pool_v2.sh`
- Output: `docs/superpowers/reports/slope_gate_pool_v2/`. The first round's
  `slope_gate_pool/` is left as it is.

- [x] **Step 1: Reproducibility (Evaluator).** For every `grid_ext` file, the 10 original
  (Z, CZ) pairs match `grid/` step for step. (2026-09-24: grid job 19918057, 120/120 completed;
  1200 pairs checked, 0 mismatches, `isequaln` on the steps structs.)
- [x] **Step 2: Pool and flag (Agent).** Pool all 540 configs per file with
  `'KeepCommonMode', true` and a tolerance of 10 samples. Build the `configs` table canonically,
  as before. For every config, run `RELAX_flag_eye_steps` on its steps (steps need `.label`), with:
  - the blink mask from `RELAX_pop_eye_mask` (HEOG off) on the restricted, impedance-screened
    montage;
  - eye labels = BlinkElectrodes union the resolved HEOG pair;
  - `MaxAmplitudeMicrovolts` 500.

  A pooled event is an **eye event** if any member step, in any config, is flagged (version B).
  It is **chained** if longer than 1 s, as before. Eye and chained events are not sampled. Record
  the counts per file and stratum.
- [x] **Step 3: Stratify and draw (Agent).**
  - Strata are strict (the maximum member Z: 5, 7, 12, 20, 35, 50) x fastLag (any member
    lag <= 2), giving 12.
  - Within each stratum, draw k = min(20, N) events uniformly without replacement, with at most
    2 events per unit (file, channel). Where the cap binds, redraw from the remaining units.
  - Record each event's stratum weight N/k for the sensitivity analysis only.
  - Add hidden repeats of about 10 % of the draws (24), shuffle and number `cand001...`.
  - `rng(20260925)`.
  - Render as in 3B: the three panels, the shaded [onset, offset] and the red dashed `xline(0)`.
  - Write `scores.csv` (png, score, notes) and `key.mat` (`key` with `member`/`count` over 540
    configs, `configs`, `strata`), plus `pool_summary.csv` and `grid_summary.csv`.
  - Put each event's eye and chained flags in `key` and the summaries.
  *Run record (2026-09-24):* pool job 19918488, 1 h 38 min, MaxRSS 36.5 GB. `tests/test_slope_gate_pool.m`
  still passes (19920626). Pooled events per file, common-mode included: m159 569,572; m208
  445,666; m497 719,218; m583 738,741. Of these, 906,794 were eye events (version B); none were
  chained. Four fastLag=0 strata are nearly or entirely empty once eye events are dropped: strict
  12 has N=15, strict 20 has N=1, strict 35 and 50 have N=0. That gives 176 draws plus 18 hidden
  repeats (10 %), so 194 candidates. The eye mask uses `RespectBoundaries=1`, as the SLURM
  parameters do.
- [x] **Step 4: Score (User).** Use `score_candidates('docs/superpowers/reports/slope_gate_pool_v2')`,
  with the same p/d/u definitions. Do not open `key.mat` or the summaries until scoring is done.
- [x] **Step 5: Validate (Evaluator).** Every row carries one of p/d/u. Then unblind.
  (2026-09-24: 194/194 scored, 19 p, 149 d, 26 u, no notes.)

**Task 4A changes (from 3D).** `analyse_scores` runs on `slope_gate_pool_v2` and `grid_ext`:
- A config detects a candidate if it has any step (common-mode included) on the candidate's
  channel within 10 samples of its sample. Rejected channels are ignored.
- The **primary estimates are unweighted**: each first score counts once. Stratum-weighted (N/k)
  estimates are a sensitivity check only.
- Precision excludes `u`. Recall is reported with `u` excluded from and added to the denominator.
- Bootstrap: resample events within strata.

The pre-registered rule (F1, the 0.8 floors, the precision walk, the tie-breakers) is unchanged;
it is applied to the unweighted estimates. The report says that the recall choice (u excluded or
included) was made after the first round.

---

### Task 4A: Precision and recall per config; choose the detection operating point

**Files:**
- Create: `slope_gate_survey/analyse_scores.m`
- Create: `docs/superpowers/reports/2026-09-23-slope-gate-survey.md` (sections 1-6 below)

- [x] **Step 1: Analysis (Agent writes, Evaluator reads).** `analyse_scores(figDir)` joins
  `scores.csv` to `key.mat`. It uses first scores only for estimates; the hidden repeats go to
  the consistency check. Outputs:
  - **Self-agreement** on the repeats: percent agreement and Cohen's kappa.
  - **Per config** (300 rows): `score_configs` precision, relative recall and unclear share, with
    a 95 % bootstrap interval (2000 resamples of units within strata), plus nSteps/min, events,
    steps per event, nChannelsBad and the median `floorBoundFraction` at that lag. Written to
    `config_table.csv`.
  - **Figures:**
    - precision against recall, one panel per (W, lag, baseline), with points labelled `Z/CZ`;
    - precision and recall heatmaps over Z × CZ per lag;
    - the pop share per stratum.
- [x] **Step 2: Apply the pre-registered rule (Evaluator proposes, User decides).** Name the
  config the rule picks and the runners-up. The user confirms or overrides; record which.
  `SlopeBaselineMs` 1000 vs 500 is read off the same table.
- [x] **Step 3: Report sections 1-6 (Evaluator):**
  1. what changed and why (the revision note above, condensed);
  2. the decision rule, verbatim, and its date;
  3. the pool summary;
  4. self-agreement;
  5. the full config table and the figures;
  6. the decision, with what the chosen point costs in recall against its neighbours and what
     the floor diagnostic said at its lag.

  **Standard: a reader with only this file can check the choice and disagree on the evidence.**
  Embed every scored `p` and `u` candidate and a representative third of the `d` candidates
  (say that is what was done), each with its per-step columns from `key.mat`.

---

**4A outcome (2026-09-24).** The round-2 analysis ran as job 19921144, cross-checked by the
evaluator. The rule, applied to the unweighted estimates with `u` excluded from recall, reached
eligibility only at Pmin 0.30, with a single config: W 151 / lag 6 / baseline 1000 / Z 7 / CZ 5
(P 0.30, R 0.81). **The user overrode it** in favour of the highest F1, with the rule's
tie-breakers: **W 151 / lag 4 / baseline 1000 / Z 20 / CZ 35** (P 1.00 on 10 of 10, R 0.62
[0.42, 0.85], F1 0.77). The user reviewed the 6 missed pops and judged them borderline or
acceptable. Report sections 1-6 hold the evidence. The provisional values in Task 3A (Z 12 /
CZ 10, window 300 ms) change in 4C.

---

### Task 4B: Downstream screens against the scored events

Detection is now fixed at the 4A point. Each screen is judged only by which **scored** events it
removes, at that point.

**Files:** `slope_gate_survey/analyse_screens.m`; report section 7; `RELAX_detect_steps.m` and
its test only if the raw-noise gate is adopted (Step 3).

- [x] **Step 1: Filter screens — DROPPED by the user (2026-09-25).** At the 4A point the scored detections are 10 p, 0 d, 0 u, so a screen could only remove pops; the user kept the config as chosen, with no screen (MinAmplitude stays 15, plateau/hold drift stay Inf, no raw-noise gate). Original text kept for the record: **Filter screens (Agent writes, Evaluator reads).** These screens only need the
  scored events' columns, so nothing is re-detected. Each value is applied alone to the events
  in the chosen config, reporting the weighted pops removed, non-pops removed and the ratio:
  - `MinAmplitudeMicrovolts` {15, 25, 50} on `|amplitude|`;
  - `MaxPlateauDriftFraction` {0.3, 0.5, 1} on `plateauDrift`;
  - `MaxHoldDriftFraction` {0.3, 0.5, 1} on `holdDrift`;
  - a raw-noise gate `RawZThreshold` {2, 3, 4, 6} on `rawZ`.

  Add one row for the per-channel noise: the pop share among events on channels in the top decile
  of `rawScale`.
- [x] **Step 2: Refinement (Agent runs, Evaluator reads).** Done 2026-09-25 (jobs 19921562, 19921583); see report §7. Run the chosen config once per file
  with `RefineTransition` `'once'` and once with `true`, keeping MaxTransitionMs / plateau / hold
  at their current values. Pool each run against the scored events (tolerance 10), then report:
  - which scored pops and non-pops each mode **rejects**;
  - how far the refined intervals move against the gate-run intervals (median and 90th
    percentile width change).

  Refinement's value to the **repair** is judged in Task 5, not here.
- [x] **Step 3: Adopt or drop each screen (User decides).** 2026-09-25: no filter screen, no raw-noise gate (Step 1 dropped); `RefineTransition` stays false. The user ruled that the bracketing code stays until Task 5 decides whether refinement helps the repair. For each screen, the evaluator
  proposes adopt (with value) or drop, from the Step 1-2 numbers. The user decides. If the
  raw-noise gate is adopted, it goes in by TDD:
  - an `inputParser` option `RawZThreshold` (default `0` = off, then set to the adopted value in
    4C), gating in `rawCandidate`/`refineCandidate` after the amplitude floor with rejection
    reason `rawNoise` in `info.nRejected`;
  - a test that a 30 uV step on a channel of raw sd 40 is rejected at `RawZThreshold` 2 and
    found at 0, written verbatim by the evaluator into this step **before** the agent runs.
- [x] **Step 3b: Common-mode check (Evaluator, User).** 2026-09-25: 16 events rendered (the files have 3, 5, 3 and 34). The user sorted them into eye artifacts and possible reference-channel pops, none a per-channel pop. `PopRepair_ExciseCommonMode` stays false. Reference pops are handled separately by the reference-pop pass, which the user still has to finish. At the chosen point, render 20
  common-mode events the way `slope_gate_survey/zhist_and_events.m` does, for the user to look
  at. At Z = CZ = 4 none of them was a pop (2026-09-24). If that still holds, excision is
  removing eye movement, drift and oscillation that later stages handle, and the user decides
  whether to keep `PopRepair_ExciseCommonMode`.
- [x] **Step 4: Channel rejection re-check (Evaluator).** 2026-09-25: 0 channels rejected on all four files after eye flagging. The detector's own pre-eye lists (m497 11, m583 13) were all frontal/periocular and were over the limit only because of blink-coincident steps. The user kept `StepsPerMinute` 6 and `TransitionFraction` 0.02 unchanged; these files do not exercise them. At the chosen point, list per file the
  channels `RELAX_rank_pop_channels` rejects and their reasons. With refinement off,
  `transitionFraction` is gate-run time and no longer bracket time, so for each rejected channel
  also report the scores of its sampled events. If the rejections no longer track channels the
  user scored as bad contacts, that is a finding and the user decides whether to re-price
  `StepsPerMinute` / `TransitionFraction`.
- [x] **Step 5: Report section 7 (Evaluator):** written 2026-09-25. the screen table, the refinement comparison, the
  decisions and who made them.

---

### Task 4C: Ship the defaults

**Files:** `RELAX_detect_steps.m`, `RELAX_Wrapper_slurm.m`, `RELAX_SET_PARAMETERS_AND_RUN_slurm.m`,
`validate_pop_repair.m`, `tests/test_RELAX_detect_steps.m` section 1h, report sections 8-10.

- [x] **Step 1: Failing test (Agent).** Red 2026-09-25, job 19921697. The evaluator rewrites section 1h's assertion with the
  4A/4B values, **verbatim in this step before the agent runs**. Its message becomes `'the
  shipped defaults are the ones the 2026-09-23 survey chose'`. Run it red. If every value
  matches the provisional one, it passes immediately; record that in the report as confirmed
  rather than inventing a change.

  **Written by the evaluator 2026-09-25.** 4A sets W 151 / Z 20 / CZ 35; lag 4 and baseline 1000
  are unchanged. 4B changes nothing else. Section 1h's second assertion becomes:

  ```matlab
    assert(P.WindowMs == 151 && P.TransitionMs == 4 && P.ZThreshold == 20 && ...
             P.ChannelZThreshold == 35 && P.MinAmplitudeMicrovolts == 15 && P.RefractoryMs == 10 && ...
             P.MaxTransitionMs == 500 && isinf(P.MaxPlateauDriftFraction) && P.HoldMs == 500 && ...
             isinf(P.MaxHoldDriftFraction) && P.CommonModeChannelFraction == 0.15 && ...
             P.SlopeBaselineMs == 1000 && P.SlopeNoiseFloorMicrovolts == 0.1 && ...
             P.SlopeSpreadQuantile == 0.85 && P.CommonModeGapMs == 10 && ...
             ~isfield(P, 'CommonModeWindowMs'), ...
          'the shipped defaults are the ones the 2026-09-23 survey chose');
  ```
- [x] **Step 2: Set the values in all four places (Agent)**. Done 2026-09-25. The four-way grep shows one value per parameter. Three failures under the new values were each ruled **pin** by the evaluator, with no assertion changed:
  - `legacy` gained Z 12 / CZ 10: section 3's ramp, since the list omitted the two values when they equalled the defaults;
  - section 2c pins Z 12: z on a silent channel caps near amplitude / noise floor, 19.6 there;
  - `test_RELAX_detect_steps_grid` section 2 reads the default pair from `info.params` instead of hard-coding `[12 10]`.

  Suite 42/42 (job 19921780). Original step text:, re-run the task test, the full
  suite, and Task 3A Step 6's four-way grep. A section-1c or 1f assertion failing under the chosen
  values is a finding, never a loosened test.
- [x] **Step 3: Report sections 8-10 (Evaluator):** written 2026-09-25.
  8. **Threats**:
     - one scorer;
     - relative recall only, bounded by the grid's union;
     - four files, one dataset, one montage;
     - the bridge screen was skipped in the survey;
     - any stratum or axis that was thin.
  9. **Shipped values**: a table of old value, provisional value, shipped value and source for
     each.
  10. **Reproduction**: the sbatch commands, job ids, log paths and the `.mat`/`.csv` paths.
- [x] **Step 4: Checkpoint (Evaluator).** Diff everything this task touched. **Do not commit.** Done 2026-09-25, nothing committed.

---

### Task 4D: Remove the unused screen machinery (added 2026-09-25, user request)

**Order (user, 2026-09-25): Task 5 runs before 4D.** The removal happens after stage validation,
so 4D can also act on Task 5's refinement verdict.

The user dropped the downstream filter screens (4B Step 1) and wants the machinery for screens
the shipped detector no longer uses removed, not left as dormant options. Runs after 4C, so the
shipped values are settled first. **The bracketing path is not removed here**: the user ruled on
2026-09-25 that it stays until Task 5 judges refinement's value to the repair, and only then is
it decided whether the bracketing path (`RefineTransition` true / `'once'`, `MaxTransitionMs`,
`RefineIterations`, `BoundaryFraction`, `BoundaryNoiseMultiple`, `SmoothMs`) is also removed.
That folds in the open decision under Task 3A about section 3's bracketing tests.

**Candidates** (grep on 2026-09-25):
- `MaxPlateauDriftFraction`, `HoldMs` and `MaxHoldDriftFraction`, with their rejection code and
  `info.nRejected` reasons. Inert with `RefineTransition` false. Referenced in
  `RELAX_detect_steps.m`, `RELAX_Wrapper_slurm.m`, `RELAX_SET_PARAMETERS_AND_RUN_slurm.m`,
  `validate_pop_repair.m`, `tests/test_RELAX_detect_steps.m`,
  `tests/test_RELAX_detect_steps_grid.m` and `tests/parked/detect_steps_plateau_hold_sections.m`.
- **Caveat:** `RELAX_detect_reference_pops.m` (the reference-pop pass, off by default) and its
  parked test also pass these options. Removal must either keep that pass working or be
  coordinated with the open item on parking or removing the reference-pop tests.
- Whether the `steps.plateauDrift` / `holdDrift` / `rawZ` / `varianceRatio` report columns stay
  as diagnostics is part of the list.
- Superseded survey drivers (`survey_hold_drift.m`) go into Task 6's list of files to decide on,
  not here.

- [x] **Step 1: Removal list (Evaluator proposes, User approves).** Every option, code path,
  `RELAX_cfg` field, test section and doc passage to remove, with its references.
  **Approved 2026-09-25 (user):**
  - Remove `MaxPlateauDriftFraction`, `HoldMs`, `MaxHoldDriftFraction` and their rejection code,
    `PopDetect_*` fields (wrapper defaults, both detector calls, entry script,
    `validate_pop_repair.m`, `survey_reference_pops.m` default `detectOpts`).
  - Remove the full iterative bracketing: `MaxTransitionMs`, `RefineIterations`, the iteration
    loop and the `noCrossing`/`unbounded` reasons. **`RefineTransition` becomes a logical:
    `true` (the default) is today's `'once'` behaviour, `false` the unbracketed diagnostic.**
    `'once'` is no longer accepted. `SmoothMs`, `BoundaryFraction`, `BoundaryNoiseMultiple`
    stay. `info.nRejected` keeps `amplitude` and `edge`.
  - Drop the `plateauDrift`, `holdDrift`, `rawZ`, `varianceRatio` columns from `steps` and from
    `RELAX_detect_reference_pops`' candidate table.
  - Tests (evaluator edits; deletions approved): drop the four columns from `numericCols`, the
    removed options from `legacy`, the `rawZ` assertion in 1f, the removed-option checks in 1h
    (1h's refinement assert becomes `RefineTransition == true`), section 3's
    plateauDrift/varianceRatio assertion; section 3 runs at `RefineTransition true` (the old
    `'once'`) with its other assertions unchanged - any failure goes back to the user. Grid option
    set 3 loses `MaxPlateauDriftFraction` and switches to `RefineTransition false`, so the grid
    test covers both modes. Delete `tests/parked/detect_steps_plateau_hold_sections.m`; delete the
    plateauDrift/holdDrift assertions in `tests/parked/test_RELAX_detect_reference_pops.m`.
  - Move the superseded drivers that use removed options/columns into `slope_gate_survey/legacy/`:
    the `slope_gate_survey/*.m` drivers that reference them, `survey_hold_drift.m`, and their
    `run_*.sh` submit scripts (`run_bracket_modes.sh`, `run_raw_run_survey.sh`,
    `run_survey_hold_drift.sh`).
  - Before-baseline: detector-only run at the defaults on the four files (job 19929524,
    `~/tmp_slope_char/equiv4d/`), taken before any edit.
- [x] **Step 2: Remove (Agent).** Code and config first. Deleting a test section that exercises
  only removed machinery is allowed once the user has approved it in Step 1; editing any other
  assertion is not. Full suite green, and a grep shows no remaining references outside
  `slope_gate_survey/legacy/`, `.superpowers/` and the reports.
  **Done 2026-09-25.** The implementer's suite run (job 19929841) was 41/42: section 3's 200 ms
  ramp gave onset 4037 (limit 4035) under single-pass bracketing. A 20-seed probe (job 19930490)
  showed ramp onsets 4018-4045 (3/20 over the limit) and one pre-ramp offset 2 samples off (limit
  1); every other section 3 shape passed on all seeds. **User decision:** widen those two
  limits to 4050 and 2, with the reason in a comment (the iterative bracketing they were written
  for is removed). Suite 42/42, job 19930557. The grep leaves only the report-record scripts
  under `docs/`.
- [x] **Step 3: Equivalence check (Evaluator).** At the shipped config on the four survey files,
  the steps table is identical before and after the removal (the removed paths were inert).
  **Done.** The detector-only baseline (job 19929524) was lost when `~/tmp_slope_char` was cleared.
  Used Task 5's full-stage `pop_validation_slopegate_once.mat` (pre-4D, `'once'`) as "before"
  instead, and re-ran `validate_pop_repair` after the removal without the ICA and legacy arms (job
  19930489). Every surviving steps column, including `eyeMovement`, is identical on all four
  files (1594 / 4843 / 1523 / 7957 steps), and so are the rejected channels and every shared
  repair-info field (`~/tmp_slope_char/equiv4d/compare4d.m`).
- [x] **Step 4: Checkpoint (Evaluator).** Diff. **Do not commit.** Not committed.

---

### Task 5: Stage-level validation and documentation

**Files:**
- Modify: `CLAUDE.md` (the stage-2 paragraph describing the detector; the fresh-output-dir note)
- Modify: `docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md` (addendum)
- Modify: `docs/superpowers/reports/2026-09-23-slope-gate-survey.md` (stage results appended)

**Interfaces:**
- Consumes: Task 4C's shipped defaults.
- Produces: no code interface; the branch is ready to review.

- [x] **Step 1: Re-run the stage validation (Agent submits, Evaluator reads)**

**Done 2026-09-25** (jobs 19921823 off, 19921824 'once'; `pop_validation_slopegate*.mat`). Load fell
about tenfold in samples touched: worst channel 0.04-0.07 % vs 0.51-0.81 %, no channels rejected,
no excision. The post-highpass residual is not worse: the fill removes 64-75 % of the excess over
control against 55-79 % on 2026-09-17 (the higher ratios come with sharper sites, raw 2.15-3.20x).
Refinement: ringing a wash, leave-one-out slightly worse, m159 ICA eye-class 0.84 -> 0.98 (7
components). **User decision: ship `RefineTransition 'once'`.** In 'once' mode MaxTransitionMs and
both drift fractions are inert, so Task 4D can still remove them; the single-pass bracketing stays.
Test 1h now asserts 'once' and checks that an exponential approach is bracketed at its mid-level
crossing.


First point `validate_pop_repair.m` at the four Sep26 files (Revision note). Run it at the
shipped defaults, and once more with `RefineTransition` switched: `'once'` if refinement shipped
off, `false` if it shipped on. Refinement changes only the `[onset, offset]` that repair fills,
so this residual comparison is where its value is decided. If refinement beats off here,
report it to the user; do not flip the default silently.

```bash
mkdir -p logs && sbatch run_validate_pop_repair.sh
```

Expected: the post-1 Hz-highpass residual at filled transitions, the leave-one-out spatial-fill
check and the ICA agreement numbers, now measured at the new detection sites. Compare against
`docs/superpowers/reports/2026-09-17-pop-repair-validation.md`. **Detection load should fall
(fewer, better-justified fills); residual quality at the sites that remain should not get
worse.** A worse residual means the gate dropped repairs that were doing real work — a finding,
not a footnote.

- [x] **Step 2: Re-check the reference-pop path (Agent submits, Evaluator reads)**

**Done 2026-09-25.** `survey_reference_pops.m` hard-coded the pre-3A operating point; it now takes
`opts.detectOpts`, and `run_survey_reference_pops_slopegate.sh` passes what the wrapper passes
(job 19922508 at refine off; rerun at 'once' after the default changed). Candidate counts fall
sharply (m497's 2026-09-21 gated set was slow-wave flanks), and no gated candidate survives
subtraction on any file. Recall of the reference path at per-channel thresholds is unmeasured: left
to the user's reference-pop work.


```bash
mkdir -p logs && sbatch run_survey_reference_pops.sh
```

`RELAX_detect_reference_pops` runs the same detector on the cross-channel median at
`RefPop_MinAmplitudeMicrovolts = 10`; the median averages out per-channel noise, so its
`channelScaleSlope` is smaller and the channel gate behaves differently there. Compare the
candidate and acceptance counts against
`docs/superpowers/reports/2026-09-21-reference-pop-survey.md`. Append the comparison to the
Task 4 report.

- [x] **Step 3: Update `CLAUDE.md` (Agent)**

In the stage-2 paragraph, replace the sentence describing the detector — "compares the medians
of two `WindowMs` windows separated by a `TransitionMs` allowance, gated by a local z (pooled
within-window robust sd), a channel-wide z ... and an absolute microvolt floor" — with the
median-slope description: the same centred median filter read across `TransitionMs` alone
(an exact lag in samples; 1 is the uncentred `mc(t+1) - mc(t)`), refinement off by default with
the gate run as the transition interval, `steps.rawZ` reported,
scored as an excess over the local baseline slope (`SlopeBaselineMs`) in local robust sd, plus
the channel-wide z on the same statistic and the unchanged floor on the long-lag level
difference. Delete the sentence describing `MaxUnstableFractionToRejectChannel` from the
channel-rejection list, note that `PopNChannelsRejectedByUnstable` is gone and that a run
against an existing `RELAXProcessed/` therefore needs a fresh output directory, and point at
`docs/superpowers/specs/2026-09-22-median-slope-gate-design.md`.

- [x] **Step 4: Add the spec addendum (Agent)**

At the end of `docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md`, add an
`## Addendum 2026-09-22 - median-slope gate` section: what changed, the shipped thresholds,
and a pointer to the new spec and the survey report. Update that file's top-of-document
"Superseded in part" banner to name 2026-09-22 as well.

- [x] **Step 5: Full suite, one last time (Agent)**

```bash
JOB=$(sbatch --wait --parsable run_all_tests.sh); tail -60 logs/run_all_tests_${JOB}.log
```

- [x] **Step 6: Checkpoint (Evaluator)**

**Done 2026-09-25.** CLAUDE.md stage-2 detection and ranking paragraphs describe the median-slope
gate with `'once'` (not "refinement off": the user's decision postdates this step's text), the
fresh-output-dir note, the new pointers, and the four-file validation set. Spec addendum
2026-09-22 added and banner updated. The entry script's stale two-window and unstable-criterion
comments were rewritten too (CRLF/ASCII kept). Report sections 9-11 updated. Suite 42/42 (job
19929369). Deferred: reference-path recall at per-channel thresholds (user's reference-pop work);
slow-pop recall on large drift (report section 8). No commit.

Read the whole diff this task produced:

```bash
git --no-pager diff -- CLAUDE.md docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md docs/superpowers/reports/2026-09-23-slope-gate-survey.md
git --no-pager status --short -- CLAUDE.md docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md docs/superpowers/reports/2026-09-23-slope-gate-survey.md
```

(`--no-pager` always: a pager called from inside MATLAB's `system()` hangs silently.)
Record in the task notes what changed, the suite count, and anything deferred. **Do not commit.** This plan produces no commits — the whole electrode-pop detection/repair line of work lands as one commit the user makes at the end.

---

### Task 6: Final review

- [x] **Step 1: Dispatch `complex-code-reviewer` (Evaluator)** — done 2026-09-25 (opus; probes in `~/tmp_slope_char/review6/`, jobs 19936549/647/651/652/681).

One agent, on the whole line of work. The user has chosen `complex-code-reviewer` for this
review and has approved its use of opus — this is the standing exception to the sonnet-only
rule for review agents. The dispatch prompt must carry, because the agent starts cold:

1. **What to review:** the working-tree diff of `RELAX_detect_steps.m`,
   `RELAX_rank_pop_channels.m`, `RELAX_Wrapper_slurm.m`,
   `RELAX_SET_PARAMETERS_AND_RUN_slurm.m`, `validate_pop_repair.m`, `slope_gate_survey/`,
   `tests/test_RELAX_detect_steps.m`, `tests/test_RELAX_rank_pop_channels.m`,
   `tests/test_RELAX_protect_eye_channels.m`, `tests/test_survey_slope_gate.m`,
   `tests/assert_smoke_outcome.m`, `CLAUDE.md` and the two docs files. **Nothing is
   committed**, so the review target is the working tree, not a commit range — say so
   explicitly and tell it to use `git --no-pager diff` and to expect unrelated uncommitted
   hold-drift and reference-pop work in the same tree, which is **not** under review.
2. **What the code is supposed to do and where that intent is recorded:**
   `docs/superpowers/specs/2026-09-22-median-slope-gate-design.md` (the spec), this plan, and
   `docs/superpowers/reports/2026-09-23-slope-gate-survey.md` (what the real data said). Name
   the published methods it claims to follow: the median filter's edge-preserving property;
   `1.4826 * MAD` as a robust sd; de Cheveigné & Arzounian (2018) for the `varianceRatio`
   column it must leave untouched; `eeg_interp 'spherical'` for the spline weights the repair
   reuses.
3. **How to run things here:** `module load matlab`; EEGLAB at `~/tools/eeglab2026.0.0`;
   MATLAB only through `sbatch` on `scu-cpu`, **never** on a login node and **never** via the
   MATLAB MCP; `sbatch --wait run_one_test.sh <test>` and `sbatch --wait run_all_tests.sh`;
   real inputs at
   `~/TMS_response_variability/store/mdd_dlpfc/preprocessed_Aug26/tms_artifact_removed`; saved
   real-data output at `slope_gate_survey.mat` and `pop_validation.mat`.
4. **Must not write to:** the input data directory, `RELAXProcessed/`, and any file outside
   the list in (1) — in particular it must not `git commit`, `git stash` or `git checkout --`
   anything.
5. **Known symptoms, as symptoms only** (no diagnosis attached, nothing off-limits):
   `tests/test_RELAX_bridged_channels.m` fails on this branch and did so before this work; any
   file where the stage rerun's residual got worse rather than better; any axis of the sweep
   where the knee was ambiguous.
6. **A floor, not a scope** — things to be sure it looks at: the index arithmetic of
   `mc(t) = mPre(t+h)` and `Dslope(t) = mc(t+g) - mc(t-g)` against the `validMask` margins and
   the segment edges; whether the baseline span can be contaminated by the step it is scoring;
   whether `loc`'s new definition still selects sensible quiet samples for
   `channelScaleSlope`; the `gPre`/`gPost` split at odd lags (above all lag 1) against
   `validMask` and both candidate paths' window positions; whether `channelScaleSlope` and `s`
   degenerate to `eps` / the floor at lag 1 because more than half of `Dslope` is exactly zero;
   whether the precision/recall estimator's weights match the two-stage draw; whether
   `steps.z`'s change from a magnitude to a signed excess is
   handled everywhere it is read; and whether the retired instability criterion left anything
   dangling in the covariate path.

- [x] **Step 2: Triage the findings (Evaluator)** — done 2026-09-25:
  - **Duplicate steps (fixed, user chose merge):** separate same-sign gate runs bracketed to the
    same transition (99-572 same-sample pairs per file). `RELAX_detect_steps` now merges a
    channel's overlapping same-sign brackets through `RELAX_merge_overlapping_steps` when
    `RefineTransition` is true. New test section 2d (failed at 12 overlaps before the fix, job
    19936891). Real data (job 19936897 vs `equiv4d`): steps 1594/4843/1523/7957 ->
    1465/4283/1401/7202, 0 overlaps left, fills and interpolated fraction unchanged, no change in
    rejected channels. Suite 42/42 (job 19936898).
  - **Self-masking of slow rises (accepted by the user, documented):** the 0.85-quantile spread
    breaks down at 15 % of the 1000 ms baseline, so linear rises >= ~150 ms are never found.
    Recorded in the help text, CLAUDE.md, the 2026-09-22 spec status and report §8.
  - **Unweighted recall (documented):** report §8 now states the weighted recall (0.38 vs 0.62),
    the tilt toward strict configs, and the corrected one-sided bound (0.74).
  - **Early-segment steps (open, rendered for the user):**
    `docs/superpowers/reports/early_segment_steps/` (job 19936894).
  - **Doc/comment fixes:** detector help (sign from Dslope, q 0.85, quantile scale, rejections
    left, raw-mode riseMs -1), CLAUDE.md (`CommonModeGapMs` is not a cfg field), spec addendum,
    stale test comments, and the `validate_pop_repair` runStage docstring (no impedance screen, no
    reference-pop pass).
  - **Not changed:** b/s windows cross segment boundaries (5-9 of 417-1598 runs differ on real
    data); parameter validation gaps (`SlopeSpreadQuantile` 0/1, `TransitionMs` 0);
    `PopRepair_MaxExciseMs` 500 in PopDefaults vs 300 in SET (inert while excision is off).
  - **Out of scope, for the user:** SET sets `OverwriteExistingCleanedFiles = 1`;
    `run_validate_pop_repair.sh` with no args overwrites `pop_validation.mat`; the wrapper's
    stale-option guard does not name the options removed here.

Fix what is real; record what is not, and why. A finding that the synthetic tests are weaker
than the real-data survey is expected and is answered by the survey, not by new synthetics. A
finding about the uncommitted hold-drift or reference-pop work is out of scope — note it for
the user and move on.

- [ ] **Step 3: Hand back to the user (Evaluator)**

Report: the suite count, the survey's retained/lost/gained tallies, the adjudication result,
the stage-validation comparison, the shipped thresholds, the review's findings and their
disposition, and the full list of paths this work touched — the user makes the single commit
covering the whole electrode-pop line of work, and needs that list to write it.

---

## Self-review notes (written with the plan)

- **Spec coverage.** Statistic and gate → Task 2. `SlopeNoiseFloorMicrovolts` and
  `SlopeBaselineMs` → Task 2 (options) and Task 4 (values). Instability removal and its seven
  call sites → Task 1. Covariate removal and the fresh-output-dir consequence → Task 1 (code)
  and Task 5 (documentation). Pinned legacy detector, step-table diff, blind adjudication,
  sweep, channel-rejection check, reference-pop re-check, stage rerun → Tasks 3-5. Every
  spec-listed unit test appears in Task 1 Step 2 or Task 2 Step 2.
- **The "rise longer than `2g`" test is written as `TransitionMs` 100 rejects / 400 accepts.**
  That asserts the lag semantics rather than a particular amplitude, so it stays true if Task 4
  moves `TransitionMs` — unless Task 4 chooses 400, in which case the test's two arms must be
  re-scaled, not deleted.
- **Task 2 Step 1's characterisation exists because a red test must actually be red.** If the
  current detector already rejects the 2 s ramp, that assertion proves nothing and the
  evaluator substitutes a shape from the characterisation log that it does accept.
- **One agent per task instead of a tester/implementer pair** (user decision 2026-09-22). The
  guard that split provided — an implementer cannot loosen a test it did not write — is
  replaced by three named requirements in the execution protocol: verbatim test code in the
  plan, a reported Red log per task, and the evaluator diffing each test file against this
  plan before accepting. Task 2 runs at effort high because its index arithmetic fails
  silently rather than loudly.
- **Nothing is committed** (user decision 2026-09-22): the whole electrode-pop line of work is
  one commit, made by the user at the end. The consequence inside the plan is Task 1 Step 7 —
  Task 3's comparison arm has to be snapshotted from the working tree at the end of Task 1,
  because there will be no commit to `git show` it from.
- **Known weak point.** `describeSteps` recomputes `b` and `s` with a copy of four lines from
  `levelStatistic`. That is deliberate duplication in survey-only code, to keep the detector's
  interface from growing a diagnostic output it does not otherwise need; if a third caller ever
  wants them, export them from `levelStatistic` instead.
- **2026-09-23 revision.** Tasks 3 and 4 were replaced by 3A-3C and 4A-4C (see the Revision
  note). There are three structural changes:
  - Detection is scored before any screen is judged, so no screen shadows the detection
    parameters.
  - One set of labels prices all 300 grid points, because candidates are pooled into events
    across configs.
  - Precision and recall replace old-vs-new counts.

  The predicted risk is lag 1. More than half of `Dslope` can be exactly zero there, which would
  make both the local spread and the channel scale degenerate. Task 3B Step 6 checks this before
  any scoring effort is spent.
