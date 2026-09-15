# Audit: electrode-pop detect-and-repair feature branch (RELAX-GEDAI)

## Header

- **Commit reviewed:** `9ca263f78b9826990090c2c1d903fa67d0e30560` (branch `tmseeg`), diff range
  `7962ea7..9ca263f` (14 commits), against the prepared package at
  `.superpowers/sdd/2026-09-14-pop-artifact-repair/review-7962ea7..9ca263f.diff`.
  `git status --short` shows only `CLAUDE.md` and `run_relax-gedai_slurm.sh` modified in the
  working tree (both explicitly out of scope per the dispatcher's constraints — these are the
  user's own uncommitted edits, including the Step-6 documentation task for this same feature,
  which is drafted but not committed). No other file under review has uncommitted changes; I
  read the working tree directly for everything else and confirmed it is byte-identical to
  `9ca263f` for the files that matter (spot-checked `RELAX_SET_PARAMETERS_AND_RUN_slurm.m`,
  which an earlier snapshot of git status had shown modified but is clean against HEAD now).
- **Components reviewed:** `vendor/noisetools/*`, `RELAX_add_vendored_noisetools.m`,
  `RELAX_channel_budget.m`, `RELAX_impedance_bad_channels.m` (4th-arg addition),
  `RELAX_detect_steps.m`, `RELAX_repair_steps.m`, the pop-repair wiring block in
  `RELAX_Wrapper_slurm.m`, the new `RELAX_cfg` fields in `RELAX_SET_PARAMETERS_AND_RUN_slurm.m`,
  `validate_pop_repair.m`, all four new/changed test files, the design spec, the plan, and the
  validation report.
- **Intent established from:** the design spec (`docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md`,
  rank-1 authority), the implementation plan, docstrings/comments in the code itself, the unit
  tests (numerically assertive, not just "ran without erroring"), and the validation report's own
  measured data (cross-checked against its saved `.mat` artifacts, see below).
- **Executed:** full test suite (`run_all_tests`, login-node-safe per CLAUDE.md), 6 purpose-built
  MATLAB probe scripts against the vendored source (each under a few seconds), and independent
  re-derivation of several of the validation report's headline numbers directly from its saved
  `detect_*.mat`/`repair_*.mat`/`sweep_m159.mat` artifacts (no new compute jobs submitted; nothing
  run on a login node; no file in the repo or in `~/TMS_response_variability` was written to).

## Hunting list

Researched from the design spec's own cited mechanism (NoiseTools `nt_regw`/`nt_detrend`),
first-principles reasoning about zero-phase-filter step distortion, and the dispatcher's
"where to focus" list. No web access was used for this review — the method (`nt_detrend`'s
weighted-mean-preserving regression) is fully specified by its own source and was verified
directly against that source and numerically, which is the strongest form of evidence
available here.

1. **Zero-weighted-mean claim (`nt_regw`/`nt_detrend`)** — checked-clean, see Finding 0 below
   (not a defect; reported per instructions as a positive confirmation).
2. **Boundary correctness: segments never cached, never spanned by a fit** — checked-clean.
   `RELAX_get_segments(EEG)` is called fresh inside both `RELAX_detect_steps` and
   `RELAX_repair_steps` on every invocation; `RELAX_get_segments` partitions `[1,EEG.pnts]`
   into contiguous, non-overlapping, gapless segments; the detector never computes a
   derivative across a segment boundary (`d(a:b-1)` per segment, transition at a boundary
   stays `NaN`); the repairer's per-channel cut-merge logic keeps step-derived cuts at least
   `MergeSamples` away from any boundary, and detection's own `SustainMs` requirement keeps
   detected steps at least `sustain` samples from any segment edge, so the two guarantees
   compose safely at shipped defaults. `test_RELAX_repair_steps_boundaries` proves
   bit-identical output to processing each side of a boundary standalone. Confirmed by
   reading + passing tests.
3. **Units: `wsize` in samples, not seconds** — checked-clean.
   `RELAX_repair_steps.m:97`: `wsize = 2*floor(round(opt.MaxWindowSeconds*EEG.srate)/2)`, matches
   `nt_detrend.m:54`'s own sample-based convention. No missed conversion.
4. **Shared budget's worst-first prefix-truncation contract** — checked-clean.
   `RELAX_detect_steps.badChannels` is explicitly sorted worst-first by step rate
   (`RELAX_detect_steps.m:213-216`, and pinned by test 8b in `test_RELAX_detect_steps.m`, which
   is specifically designed to catch a regression to channel-index order). The impedance
   screen is wired through the documented bypass (4-argument `maxRejectableCount` form), so
   its own index-order return value never reaches `'spend'`'s truncation logic — confirmed by
   reading `RELAX_Wrapper_slurm.m:486-489` and `RELAX_impedance_bad_channels.m:125-150`.
5. **Aggregate-table safety** — checked-clean for crash risk, minor semantic inconsistency found
   (Finding 3). `RELAX_append_table_row` (pre-existing, unmodified by this branch) already
   reconciles schema drift between rows (backfills missing columns with NaN/`{}` on either
   side, matches by name not position) — this is exactly the mechanism that makes the
   pop-stage's fields being present on some files/absent on others (`Perform_pop_artifact_repair`
   off by default, `PopScreenInfo.available`/`PopRepairInfo.available` both independently
   gate which fields get set) safe rather than a table-building crash. All Pop* fields are
   scalars or cell arrays of labels — none is a shape that would defeat `struct2table(...,
   'AsArray',true)`.
6. **Validation report's honesty** — checked-clean. See "Validation report" section below.
7. **Vendored NoiseTools closure correctness / dependency shadowing** — Finding 1
   (Important, confirmed by execution).
8. **`chooseOrder`/`useOverlapAdd` consistency when `MaxPolyOrder` is reconfigured** — Finding 2
   (Important, confirmed by execution).
9. **Entry-point disabling of `Perform_pop_artifact_repair`** — checked-clean.
   `RELAX_SET_PARAMETERS_AND_RUN_slurm.m:340` sets `=0`; `RELAX_Wrapper_slurm.m:255-257`'s
   absent-means-off guard is a second, independent safety net; grep across the whole repo
   found no other place that sets this field to `1`; `RELAX_Wrapper.m`/`RELAX_Wrapper_bart.m`
   never reference it at all (matches the design's "`RELAX_Wrapper_slurm.m` only" scope).
10. **Population A / channel-index remapping after `pop_select`** — checked-clean.
    `RELAX_Wrapper_slurm.m:298-303` re-derives `PopSteps.channel` from `.label` via `ismember`
    against the post-`pop_select` `EEG.chanlocs`, rather than reusing stale indices — exactly
    the index-invalidation pattern this kind of review exists to catch, and it is handled
    correctly here.

## Finding 0 — the zero-weighted-mean claim: CONFIRMED CORRECT

This is the claim the whole design rests on, so I verified it three ways.

**By reading.** `nt_regw.m:29-32`: `mn = y - nt_demean(y,w)` (the weighted mean of the data being
fit). The returned trend estimate `z` is built as `nt_demean(x,w)*topcs*b + mn` (lines 43/58/76) —
a linear combination of *demeaned* regressor columns, plus `mn` added back. Since a weighted mean
is a linear functional and `nt_demean(x,w)` has zero weighted mean by construction, the regression
term contributes zero weighted mean to `z`, so `z`'s weighted mean equals `mn`, i.e. the *data's*
weighted mean. Therefore `x - z` (what `nt_detrend_helper` returns as the detrended residual,
`nt_detrend.m:261`) has zero weighted mean over the fit interval, for any polynomial order and any
weight vector, exactly as the design's docstring/comments in `RELAX_repair_steps.m:11-15` claim.

**Numerically, against `nt_regw` directly** (`probe_zero_mean.m`, saved at
`/home/nga4004/.claude/jobs/91008025/tmp/probe_zero_mean.m`): unweighted residual mean
`-5.7e-17`; weighted residual mean `1.9e-15` against a weight vector with real zeros and partial
weights — both at floating-point noise floor, while the *plain* (unweighted) mean of the same
weighted-fit residual is `-2.5e-2`, confirming the guarantee is specifically about the *weighted*
mean, as claimed, not a coincidence of unweighted symmetry.

**Numerically, end to end, on a genuine step** (same probe): a synthetic channel with a real
+400 µV step at sample 1001, split at the step sample into two pieces (never touching the
amplitude), each independently run through `nt_detrend(..., thresh=3, niter=3)` — the exact call
RELAX_repair_steps makes. Result: `median(pre-repaired) = -0.126`, `median(post-repaired) =
-0.115`, **difference 0.010 µV**, not 400 µV. The step offset disappears as a side effect of
per-interval detrending alone, with no explicit subtraction anywhere in the call. This is a direct
demonstration of the design's central mechanism working as claimed against the actual vendored
source shipped in this branch.

One minor nuance surfaced by the same probe, not a bug: with the default `niter=3` (iterative
reweighting), the *plain* residual mean is not exactly machine-zero (`-6.9e-3` on a ~3-8 µV-scale
signal) because a few samples get down-weighted as outliers across iterations, and the guarantee is
about the weighted mean under whatever `w` produced that specific fit, not the plain mean. The
residual is three orders of magnitude below the step sizes this stage targets and does not
materially affect anything.

## Findings

### Finding 1 — `RELAX_add_vendored_noisetools()`'s self-check does not cover its own dependency closure; a second, incompatible vendored copy of 5 of those 9 files exists in this same repo and can silently shadow them

- **Severity:** Important (latent — confirmed reachable by execution, confirmed currently
  numerically benign for the actual calling pattern; would become a live correctness bug the
  moment the calling pattern changes without anyone noticing, since nothing detects it).
- **Confidence:** Confirmed by execution (`probe_shadow.m`, saved at
  `/home/nga4004/.claude/jobs/91008025/tmp/probe_shadow.m`) plus the `run_all_tests` run's own
  console output.
- **Location:** `RELAX_add_vendored_noisetools.m:18-24` (the self-check); the shadowing source
  is `vendor/zapline-plus/nt_demean.m`, `nt_fold.m`, `nt_unfold.m`, `nt_greetings.m`,
  `nt_version.m` (all pre-existing, not touched by this branch) versus
  `vendor/noisetools/nt_demean.m`, `nt_fold.m`, `nt_unfold.m`, `nt_greetings.m`, `nt_version.m`
  (new in this branch).
- **What's wrong.** `vendor/zapline-plus/` (unrelated, pre-existing vendored tree) bundles its
  own copies of `nt_demean.m`, `nt_fold.m`, `nt_unfold.m`, `nt_vecadd.m`, `nt_vecmult.m`,
  `nt_greetings.m`, `nt_version.m` — five of which (`nt_demean`, `nt_fold`, `nt_unfold`,
  `nt_greetings`) are **not byte-identical** to the new `vendor/noisetools/` copies (confirmed
  by direct diff; e.g. zapline's `nt_demean.m` calls `nt_unfold(x)` unconditionally where
  noisetools' vendored copy only does so when `ndims(x)==3`; zapline's `nt_fold.m` is a
  different implementation entirely, down to the internal function name (`fold` vs `fold`,
  different truncation-warning behavior on non-integer epoch counts)). `RELAX_add_vendored_noisetools()`'s
  self-check only verifies `which('nt_detrend')` resolves into `vendor/noisetools/` — it never
  checks `nt_regw`, `nt_demean`, `nt_fold`, `nt_unfold`, or `nt_greetings`, even though
  `nt_detrend` calls all of them transitively. Since `nt_detrend.m`'s internal calls to
  `nt_regw`/`nt_demean`/`nt_fold` are resolved by MATLAB's normal path search at call time (not
  pinned to the calling file's own directory), if `vendor/zapline-plus/` ever ends up ahead of
  `vendor/noisetools/` on the path (confirmed to happen mechanically the moment
  `RELAX_add_vendored_zapline()` is called again after `RELAX_add_vendored_noisetools()` — which
  is exactly what several of the branch's own new/pre-existing test files do in sequence within
  `run_all_tests`, e.g. `test_pwelch_segmented.m` calls `RELAX_add_vendored_zapline()` after
  `test_RELAX_repair_steps.m` has already called `RELAX_add_vendored_noisetools()`), then
  `nt_demean`/`nt_fold`/`nt_unfold`/`nt_greetings` silently resolve to the zapline-plus copies for
  the rest of the MATLAB session, and `RELAX_add_vendored_noisetools()`'s own self-check would
  report success regardless, because it never looks at those functions.
- **Confirmed by execution.** `probe_shadow.m` reproduces the shadowing mechanically
  (`which('nt_demean')` and `which('nt_fold')` switch from `vendor/noisetools/` to
  `vendor/zapline-plus/` after a second, out-of-order `addpath` call exactly like the one the
  test suite performs; `nt_detrend` itself keeps resolving correctly since zapline has no
  `nt_detrend.m`/`nt_regw.m` of its own) — and this is directly visible as a side effect in the
  plain `run_all_tests` output, which prints **two different NoiseTools citation banners**
  (`version 29-Apr-2023` from the vendored copy, then later `version 25-Oct-2020` from
  zapline-plus's bundled copy) during a single test run, proving two different `nt_greetings.m`
  files actually executed in the same session.
- **Currently benign, but only by luck of the calling pattern, not by design.** The same probe
  compared `nt_detrend`'s numerical output on identical input under the correct path order versus
  the shadowed order: **bit-identical** (`isequal` true, max abs diff `0`). This is because
  `RELAX_repair_steps` only ever calls `nt_detrend` on a single-channel `Tx1` column vector
  (never a genuinely 3-D/multi-epoch array), and at `ndims==2`/`nepochs==1` both the vendored and
  the zapline-plus versions of `nt_demean`/`nt_fold`/`nt_unfold` happen to take their respective
  no-op branches. If that calling pattern ever changes (e.g. a future batched multi-channel call,
  or a future NoiseTools upstream sync that changes 2-D behavior in either tree independently),
  this would silently diverge with no error and no test catching it — the existing
  `test_RELAX_add_vendored_noisetools` test only repeats the same `nt_detrend`-only check.
- **Scientific impact.** None currently (confirmed numerically benign for the actual usage, and
  the pop-repair stage is disabled by default besides). This is a robustness/reproducibility gap,
  not an active bias.
- **Fix (pick one):** (a) extend `RELAX_add_vendored_noisetools()`'s self-check to verify
  `which('nt_demean')`, `which('nt_fold')`, `which('nt_unfold')`, `which('nt_greetings')` all
  resolve into `vendor/noisetools/`, not just `nt_detrend`; or (b) have
  `RELAX_repair_steps.m` call `RELAX_add_vendored_noisetools()` itself immediately before its
  first `nt_detrend` call (defensive re-assertion at point of use, catching any path
  perturbation that happened between setup and the stage actually running); (a) is cheaper and
  catches the failure mode directly at its source.

### Finding 2 — `RELAX_repair_steps`: raising `PopRepair_MaxPolyOrder` above 3 silently disables the overlap-add windowing for every long sub-segment, with no error or warning

- **Severity:** Important (dormant at the shipped default of 3; silently and completely breaks
  a stated design invariant the moment this exposed, documented config knob is changed to
  anything above 3 — which is exactly the kind of change the validation report's own
  recommended "rederivation pass" would plausibly touch, since it explicitly calls out
  retuning several sibling parameters together).
- **Confidence:** Confirmed by execution (`probe_maxpolyorder.m`, saved at
  `/home/nga4004/.claude/jobs/91008025/tmp/probe_maxpolyorder.m`).
- **Location:** `RELAX_repair_steps.m:141` (`order = chooseOrder(...)`) and `:149`
  (`useOverlapAdd = (order == opt.MaxPolyOrder) && ...`), interacting with `chooseOrder`'s
  hardcoded `C = [1.69 1.24 1.09]` (3 entries) at `:79` and `:226` (`for n =
  1:min(maxOrder, numel(C))`).
- **What's wrong.** `chooseOrder` can never return an order above `numel(C) = 3`, regardless of
  `opt.MaxPolyOrder`. `useOverlapAdd` is gated on `order == opt.MaxPolyOrder` (the raw configured
  value), not on `order` having reached its *actual* achievable ceiling. At the shipped default
  (`MaxPolyOrder = 3`), these coincide and everything works as documented (verified: a 10 s
  no-step channel correctly lands on `order3oa`). The moment `MaxPolyOrder` is configured to 4 or
  5, `chooseOrder` still only ever returns at most 3, so `order == opt.MaxPolyOrder` (`3 == 4` or
  `3 == 5`) is **never true**, and `useOverlapAdd` is false for every sub-segment regardless of
  length. On the same 10-second no-step test channel, `MaxPolyOrder=4` and `MaxPolyOrder=5` both
  silently fall back to a **single global order-3 polynomial fit over the entire 10-second
  segment** (`treatmentCounts.order3oa=0, order3=1`) instead of the intended windowed
  overlap-add — directly violating design decision 5's stated invariant ("a maximum effective
  window so the cutoff never falls below 0.5 Hz"): a single global order-3 fit over a segment many
  times longer than `MaxWindowSeconds` has an effective cutoff far below 0.5 Hz and is a much
  poorer high-pass than what the ladder is documented to guarantee. Lowering `MaxPolyOrder` below
  3 (e.g. to 2) does **not** trigger the bug — confirmed the equality condition self-corrects
  correctly in that direction (`order` saturates at the lower `opt.MaxPolyOrder` itself, so the
  equality holds) — the failure is one-directional, only above 3.
- **A related, purely cosmetic consequence of the same root cause**: when `useOverlapAdd` *is*
  correctly triggered at a non-default `MaxPolyOrder` (e.g. `MaxPolyOrder=2`, confirmed via the
  same probe), the counter that gets incremented is still hardcoded `counts.order3oa` regardless
  of the actual order used — so `treatmentCounts.order3oa` would misreport an order-2 overlap-add
  as an order-3 one if this parameter is ever changed. Not independently scientifically
  consequential, but worth fixing alongside the main issue since it's the same code path.
- **Scientific impact if triggered:** on real ~8 s inter-boundary TMS-EEG segments (the actual
  data this stage targets, per the design's own problem statement), essentially every sub-segment
  without steps would fall into this regime once `MaxPolyOrder>3`, silently defeating the
  frequency-response guarantee the spectral-fidelity test (`test_RELAX_repair_steps.m` test 2)
  exists to certify — a global low-order fit over a long window leaves substantially more slow
  drift unremoved than the documented ladder promises, in the direction of *under*-cleaning
  low-frequency artifact (not over-cleaning neural signal), but silently and without any
  diagnostic signal distinguishing it from correct operation.
- **Not currently live**: the shipped `RELAX_cfg.PopRepair_MaxPolyOrder = 3` matches `numel(C)`
  exactly, and the stage is disabled by default regardless. This is a landmine for the
  "rederive thresholds" pass the validation report explicitly recommends as the next step,
  not a currently-active defect.
- **Fix (pick one):** (a) simplest and most robust — change the `useOverlapAdd` condition to
  compare against the *achievable* ceiling rather than the raw config value:
  `useOverlapAdd = (order == min(opt.MaxPolyOrder, numel(C))) && ...`; (b) extend `C` to cover
  higher orders if higher-order fits are ever actually wanted; (c) at minimum, have
  `RELAX_repair_steps` warn or error if `opt.MaxPolyOrder > numel(C)` so the mismatch cannot pass
  silently. (a) is a one-line fix and preserves all current behavior at the default.

### Finding 3 — `PopMedianStepAmplitude` covariate is computed after Population-A channel rejection filters the step list, while `PopNSteps`/`PopNChannelsAffected` are computed before it

- **Severity:** Minor (no current scientific impact — the stage is disabled by default — but
  directly touches the covariate-integrity concern CLAUDE.md calls out by name for this exact
  kind of stage, so worth recording for whoever re-enables it).
- **Confidence:** Confirmed by reading (inferred intent — the design spec's "Metrics and
  covariates" section names the covariate but does not specify whether it should reflect the
  full detected population or only the retained/repaired one, so this is a genuine ambiguity in
  the spec, not a directly contradicted contract).
- **Location:** `RELAX_Wrapper_slurm.m:293-334`.
- **What's wrong.** `PopChannelToReject` (Population A, chronically-popping channels rejected
  outright) is computed at line 293; if non-empty, `PopSteps` is filtered at lines 298-299 to
  drop every step belonging to a rejected channel, *before* `PopMedianStepAmplitude =
  median(abs(PopSteps.amplitude))` is computed at line 334. Meanwhile `PopNSteps` (line 332) and
  `PopNChannelsAffected` (line 333) are taken directly from `PopScreenInfo`, which reflects the
  detector's *original*, unfiltered output. So on a file with heavy Population-A activity (the
  design's own survey found Population-A step amplitudes running far larger than Population B's —
  m155's dominant channel had median step size 1365 µV versus Population B's 50-1500 µV,
  generally at the low end), `PopNSteps`/`PopNChannelsAffected` measure "how poppy was this
  session overall" while `PopMedianStepAmplitude` measures "how big were the steps on the
  channels that were *kept and repaired*" — two different populations answering two different
  questions, filed under covariate names that both read as unqualified session-level summaries.
- **Scientific impact if re-enabled**: a downstream model using these as covariates (exactly the
  use case CLAUDE.md's "Saving covariates for downstream ML classification tasks" section
  anticipates) could see `PopMedianStepAmplitude` *drop* on a session with heavier Population-A
  contamination (because the largest-amplitude population gets filtered out of the median before
  it's computed), which is the opposite of what the covariate's name suggests it should do,
  potentially masking rather than exposing the session-quality signal it exists to capture.
- **Fix (pick one):** (a) compute `PopMedianStepAmplitude` from the *unfiltered* `PopScreenInfo`-era
  step list (before the Population-A filter), matching `PopNSteps`; or (b) keep the current
  post-filter semantics but rename the field (e.g. `PopMedianRepairedStepAmplitude`) and add a
  second covariate for the pre-filter population if both numbers are useful; either is a
  one-line/one-field change. No fix is urgent while the stage stays disabled.

## Validation report: honesty check

Cross-checked the report's headline numbers directly against its own saved artifacts in
`/home/nga4004/.claude/jobs/91008025/tmp/` (not re-run — used the existing `detect_*.mat`,
`repair_*.mat`, `sweep_m159.mat`):

- Section 1's per-file `nSteps`/`nChannelsAffected` table: reproduced exactly from
  `detect_*.mat` for all 5 files (e.g. m155: 1,390,837 steps / 248 channels; m159: 177,999 /
  184; m208: 505,514 / 198; m497: 151,638 / 237; m583: 145,066 / 250).
- Section 2's budget table: reproduced exactly (`nBadWanted`/`nBadFit`/`remainingBudgetForImpedance`
  fields in the saved structs match the report's "pop screen wants" / "budget allows" / "left for
  impedance" columns row for row, including m155's 248/64/0).
- Section 3's m159 parameter sweep: reproduced exactly from `sweep_m159.mat` for all 10 rows
  (defaults through `MinRiseFraction=0.9`), including the `ZThreshold=200` row (16,522 steps /
  157 channels / 22 bad-contact channels) the report calls out as the best single-knob match.
- Section 4's ringing table: reproduced exactly from `repair_*.mat` for all 5 files (raw/repaired/
  control medians and residual ratios, e.g. m155: 225.5/223.9/185.2 → 1.21×; m159:
  18.4/14.9/12.2 → 1.22×), including `treatmentCounts.fallback = 0` on every file (matching the
  report's Section 5 table) and repair elapsed times matching to the tenth of a second.

No discrepancy found between any number quoted in the report's prose/tables and the underlying
saved data. The report's conclusions are supported by its own evidence: it correctly identifies
`ZThreshold` as the dominant cause of over-detection without overclaiming a fix (explicitly
flags the sweep as single-file, and flags the `MinAmplitudeMicrovolts` deficit as an unresolved
second problem); it correctly declines to read the 1.09-1.22× ringing ratios as reassurance,
giving the specific mechanistic reason (dilution by false positives) and the specific tell
(m155's 0.7% effect on an essentially-all-false-positive file); and it correctly disables the
stage at the entry point rather than merely recommending it be disabled. I found nothing in the
report that overstates what was measured, and the one deliberately-left-open question (guard-band
ringing, design decision 7) is stated as open, not resolved either way, matching what the
underlying (contaminated) measurement actually supports.

## Unverified concerns

- **`RELAX_detect_steps`'s "computed per channel over the WHOLE recording" claim for the
  median/MAD scale** (docstring line 63) is, on inspection, actually computed only from samples
  inside inter-boundary segments long enough to hold a full sustain window on both sides (a
  `minSegLen` gate filters out too-short segments before the median-filtered-derivative array is
  built). This is inconsequential in practice for the real data (61 segments of ~8 s each,
  `minSegLen` at defaults is ~0.2 s), and arguably the *safer* behavior, but the docstring is
  technically imprecise. Not raised as a formal finding since I could not construct a real
  scenario where it changes behavior meaningfully — would need a file with a very short
  inter-boundary segment mixed among normal-length ones to see whether the scale estimate shifts
  in a way that matters.
- **`RELAX_repair_steps`'s defensive sub-segment "clip" logic** (`RELAX_repair_steps.m:124-128`,
  comment "a cut may fall inside a later boundary segment") appears, by construction, to be
  unreachable given the guarantees `RELAX_detect_steps`'s `SustainMs` requirement and
  `RELAX_repair_steps`'s own `MergeSamples` merge logic jointly provide at any config where
  `MergeSamples <= 2*refractory` roughly holds (true at shipped defaults). I could not construct
  a config combination that makes it actually fire and tested none directly (it would require
  deliberately setting `RELAX_detect_steps`'s `SustainMs`/`RefractoryMs` and
  `RELAX_repair_steps`'s `MergeSamples` to values that violate that relationship, which no
  shipped config does). If it never fires, it's dead defensive code; if it does fire under some
  unexamined config combination, verify whether it can silently drop a strip of samples from
  `repaired(ch, a:b)` coverage (i.e. leave them at their raw, undetrended value) rather than
  reassigning them to an adjacent sub-segment.

## Coverage and limits

Read in full: `RELAX_detect_steps.m`, `RELAX_repair_steps.m`, `RELAX_channel_budget.m`,
`RELAX_add_vendored_noisetools.m`, `RELAX_impedance_bad_channels.m`, all of
`vendor/noisetools/nt_detrend.m`/`nt_regw.m`/`nt_demean.m` plus enough of the rest of the closure
(`nt_fold.m`, `nt_unfold.m`) to assess the dependency-shadowing question, `vendor/noisetools/VENDORING.md`,
the pop-repair block and its call sites in `RELAX_Wrapper_slurm.m`, the new config block in
`RELAX_SET_PARAMETERS_AND_RUN_slurm.m`, `RELAX_get_segments.m` (for the boundary-partition
contract), `RELAX_append_table_row.m` (pre-existing, to assess aggregate-table safety), the design
spec (all 493 lines), the validation report (all 357 lines), all four new/changed test files, and
`validate_pop_repair.m`. Skimmed (not deeply audited, since out of the branch's algorithmic core
and already covered by the report's own honesty check): the implementation plan
(`docs/superpowers/plans/2026-09-14-pop-artifact-repair.md`, 1945 lines) and the task
briefs/reports under `.superpowers/sdd/`.

Ran: the full test suite (`run_all_tests`, 24/25 passing, the one failure being the documented
pre-existing non-hermetic `eeg_emptyset` issue — no regression from this branch); 6 targeted MATLAB
probes against the vendored source, all on synthetic data, each a few seconds of compute, all in
`/home/nga4004/.claude/jobs/91008025/tmp/` (`probe_zero_mean.m`, `probe_shadow.m`,
`probe_maxpolyorder.m`, and three small `.mat`-artifact-inspection scripts); no probe modified
anything outside scratch, no file under `~/TMS_response_variability` was touched, nothing was
submitted to a job queue, nothing ran on a login node.

Did not run: `validate_pop_repair.m` itself (per instructions — it's a ~50-minute job and the
saved artifacts from the prior run already answer every question I had; I independently
re-derived several of the report's tables from those artifacts instead of re-running it) and did
not attempt to run the full `RELAX_Wrapper_slurm.m` end to end on real data (out of scope for a
component-level audit and would require write access I was told not to use).

No findings at Critical severity. Two Important findings, both currently dormant under shipped
defaults (one a robustness/self-check gap with demonstrated-benign current numerical impact, one
a silent-degradation bug reachable only by a config change no one has made yet) and one Minor
covariate-semantics finding, plus a positive, numerically-confirmed verification of the design's
central mathematical claim.
