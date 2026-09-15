# Electrode-pop repair: real-file validation and the ringing question

Date: 2026-09-15
Scope: `RELAX_detect_steps.m` / `RELAX_repair_steps.m` (Tasks 1-5), exercised directly (not
through `RELAX_Wrapper_slurm.m`) on the five files the design's exploratory survey characterised.
Companion to `docs/superpowers/specs/2026-09-14-pop-artifact-repair-design.md`.

**Headline: the stage fails validation at its shipped defaults. Recommend
`RELAX_cfg.Perform_pop_artifact_repair` stay OFF until the detector's thresholds are
rederived against real data.** Every number below is measured, not modelled.

**This is not hypothetical: the stage WAS enabled in the shipped entry point.**
`RELAX_SET_PARAMETERS_AND_RUN_slurm.m:336` read `RELAX_cfg.Perform_pop_artifact_repair =
1` on HEAD before this report — set that way by the Task 5 wiring commit — so running the
pipeline as shipped executed the stage at exactly the defaults this report proves broken.
This report's commit changes that line to `= 0` with a comment pointing back here; every
`PopDetect_*`/`PopRepair_*` parameter line stays in place unchanged, so the stage is one
flag away from being re-enabled once its thresholds are rederived (Section 3). **If you
maintain your own copy of `RELAX_cfg`** (a script derived from this one, a notebook that
sets the struct by hand, etc.), check it explicitly: omitting the field is safe — the
wrapper's own absent-means-off guard at `RELAX_Wrapper_slurm.m:255-257` covers that case —
but a copy that carries the old `= 1` forward is not, and this report's fix does not reach
that copy.

## Reproducibility

- SLURM job `19639307`, partition `sackler-cpu`, node `c6`, `COMPLETED`, exit `0:0`,
  elapsed `00:50:54` (`sacct -j 19639307`).
- Driver: `/home/nga4004/.claude/jobs/91008025/tmp/run_pop_validation.sbatch` (scratch,
  not part of this repo) — `sbatch -p sackler-cpu -c 4 --mem=64G --time=8:00:00`, log at
  `/home/nga4004/.claude/jobs/91008025/tmp/pop_validation_19639307.log`.
- The sbatch script runs three phases, each file/combination as its **own** `matlab
  -batch` process under a shell `timeout`, so a slow or crashing file cannot take the
  others down and every result is saved to its own `.mat` immediately on completion:
  1. detection only, all 5 files, shipped defaults (`pop_detect_one.m`), plus a
     simulation of the shared channel-rejection budget exactly as
     `RELAX_Wrapper_slurm.m` spends it for the pop screen;
  2. a one-parameter-at-a-time sweep on `m159` (`pop_sweep_m159.m`);
  3. repair + the ringing measurement, one file per process
     (`pop_repair_one.m`), logic identical to `validate_pop_repair.m` below.
- `validate_pop_repair.m` (repo root) is the single-process reference script this task
  specifies; it implements the same detection/repair/ringing logic used in phases 1 and
  3 above. An earlier attempt to run it directly inside an interactive `srun` session
  died silently when the session's shell exited mid-run (no SLURM accounting record, no
  error in its log) partway through `m159`'s repair — the batch job above is the actual
  source of every number in this report, run as an independent, monitorable SLURM job
  rather than tied to an interactive shell. This is worth recording as an operational
  lesson for future validation runs on this cluster: always `sbatch`, never rely on an
  `srun --pty` session surviving unattended.
- EEGLAB `~/tools/eeglab2026.0.0`, MATLAB R2023a (`module load matlab/R2023a`).
- Baseline test suite immediately before this run: `24/25 test files passed` via
  `matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests'));
  run_all_tests"`, the one pre-existing failure being the documented
  `test_RELAX_excluding_channels_and_epoching_boundary` / `eeg_emptyset` non-hermetic-suite
  issue. No regression.

## 1. Detection at shipped defaults vs. the exploratory survey

All five files loaded and ran to completion; nothing crashed. `RELAX_detect_steps` was
called with no arguments (i.e. every shipped default:
`MedianFilterOrder=9, ZThreshold=20, MinAmplitudeMicrovolts=10, RefractoryMs=20,
SustainMs=100, MinRiseFraction=0.5, StepsPerMinuteToRejectChannel=30`).

| file | steps (production) | steps (survey) | ratio | channels affected (prod/257) | channels affected (survey) | flagged bad-contact (prod) |
| --- | ---: | ---: | ---: | --- | --- | ---: |
| m159 | 177,999 | 15,045 | **11.8x** | 184 | 169 | 166 |
| m497 | 151,638 | 8,387 | **18.1x** | 237 | 203 | 138 |
| m208 | 505,514 | 5,269 | **95.9x** | 198 | 251 | 193 |
| m583 | 145,066 | 2,110 | **68.8x** | 250 | 215 | 171 |
| m155 | 1,390,837 | 635 | **2190.3x** | 248 | **2** | 248 |

Every file over-detects, by an order of magnitude to over three orders of magnitude. The
design explicitly predicted the opposite: "The production detector adds the rise-time
criterion the exploratory one lacked... the production counts should be lower; if they
are higher, something is wrong" (design spec, and repeated in the task brief). They are
higher on all five files. Something is wrong — see Section 3 for the diagnosis.

**m155 is the decisive case.** The survey's whole point in including m155 was as the
"pure Population A" example: one chronically bad contact (E91, 624/635 = 98% of all
steps) on an otherwise clean 256-channel net (2/257 channels affected total). At shipped
defaults, the production detector instead flags **248 of 257 channels** (96.5% of the
montage) as step-affected, and **all 248** exceed the 30 steps/min bad-contact threshold.
The survey's central structural claim — two well-separated populations, a handful of bad
contacts vs. a large clean majority — is not reproduced on this file at these settings.
This is not a borderline disagreement; it inverts the survey's headline result on its own
cleanest example.

## 2. Shared channel-rejection budget: the pop screen starves impedance on every file

`RELAX_SET_PARAMETERS_AND_RUN_slurm.m` sets
`RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted = 0.25`, so
`maxRejectable = floor(0.25 * 257) = 64` for all five files (same montage). Simulated by
calling `RELAX_channel_budget(EEG,'init',0.25)` then
`RELAX_channel_budget(EEG,'spend','pop',info.badChannels)`, exactly as the wrapper does
before the impedance screen would run:

| file | pop screen wants to reject | budget allows (of 64) | budget exhausted | left for impedance |
| --- | ---: | ---: | :---: | ---: |
| m159 | 166 (64.6% of montage) | 64 | yes | **0** |
| m497 | 138 (53.7%) | 64 | yes | **0** |
| m208 | 193 (75.1%) | 64 | yes | **0** |
| m583 | 171 (66.5%) | 64 | yes | **0** |
| m155 | 248 (96.5%) | 64 | yes | **0** |

On every one of the five files, the pop screen alone wants more channels than the entire
shared budget permits, so it consumes the budget completely (`exhausted=true`) and leaves
the impedance screen **zero** allowance downstream, regardless of how many
impedance-failing channels that screen would otherwise flag. `RELAX_channel_budget`
itself behaves exactly as designed (worst-first truncation, `exhausted` flag set,
warning fired) — this is the budget doing its job on inputs the pop screen should never
have generated at this volume.

## 3. Root cause: `ZThreshold`, and it is fixable — but only demonstrated on one file

**The shipped defaults, stated explicitly, since naming them is what makes the fix
obvious:** `RELAX_cfg.PopDetect_ZThreshold = 20` and
`RELAX_cfg.PopDetect_MinAmplitudeMicrovolts = 10`. Both are FieldTrip `ft_artifact_jump`'s
own defaults, carried over by design decision 2 as "well established and trivially
defensible to reviewers" — a reasonable starting point for a detector never before run
against this data, but not validated against it before shipping. The sweep below shows
the gate needs to move roughly 10x, to `ZThreshold≈200`, before this data's own noise
floor stops dominating the count.

One-parameter-at-a-time sweep from shipped defaults, m159 only (`pop_sweep_m159.m`;
survey for this file: 15,045 steps / 169 channels):

| setting | steps | channels affected | bad-contact channels |
| --- | ---: | ---: | ---: |
| defaults | 177,999 | 184 | 166 |
| `RefractoryMs=100` | 92,910 | 184 | 156 |
| `RefractoryMs=1000` | 24,999 | 184 | 0 |
| `MinAmplitudeMicrovolts=50` | 2,606 | 112 | 0 |
| `MinAmplitudeMicrovolts=92` (= survey's own median step for this file) | 666 | 36 | 0 |
| `MinAmplitudeMicrovolts=150` | 265 | 21 | 0 |
| `ZThreshold=100` | 35,963 | 183 | 55 |
| **`ZThreshold=200`** | **16,522** | **157** | 22 |
| `RefractoryMs=100 + MinAmplitudeMicrovolts=92` | 409 | 36 | 0 |
| `MinRiseFraction=0.9` | 114,807 | 183 | 132 |

`ZThreshold=200` — ten times the shipped default of 20 — reproduces the survey on **both**
marginals at once: 16,522 vs. 15,045 steps (10% high), 157 vs. 169 channels (7% low). No
other single knob tested matches on both axes simultaneously.

Mechanistically: `ZThreshold=20` is FieldTrip's `ft_artifact_jump` default, calibrated by
convention where z=5-10 already flags an extreme outlier. On this data the per-channel
MAD of the median-filtered derivative is tiny (a spot check on m155's channel E1, run
ad hoc outside the three sbatch phases — `diag_m155.m`, not one of the saved-artifact
figures elsewhere in this report, though its arithmetic is exactly
`RELAX_detect_steps.m`'s own `med`/`mad`/`z` definitions applied by hand: median
and MAD of `abs(diff(movmedian(x,9)))` both ≈2.06 uV, i.e. the derivative distribution is
sharply peaked near a small, non-zero value — plausibly quantisation or continuous
low-amplitude EMG/instrument noise riding under the EEG). At z=20 that peaked scale puts
the absolute threshold at only ≈63 uV for a single-sample jump, well inside the range of
this data's ordinary noise floor and any residual muscle activity in these post-TMS-train
segments (61 inter-boundary segments on m155, matching the survey's segment count exactly
— boundary handling is not the issue here). Against a scale this fine, genuine 50-1500 uV
pops sit hundreds of MADs out, so nothing is lost by pushing the threshold far higher;
noise sits much closer in, so a 10x threshold increase removes it disproportionately.

**The amplitude floor tells the same story from a different angle, and it disagrees with
the z-threshold story about *which* events survive.** At `MinAmplitudeMicrovolts=92` — the
survey's own reported *median* |step| for this file, so roughly half its 15,045 steps
(~7,500) should clear that floor — the detector returns only **666**. That is an ~11x
*deficit* at the same time `ZThreshold=20` alone produces a ~12x *excess*. The two
instruments are not disagreeing only about how many events cross a count threshold; they
disagree about what amplitude the real population of pops actually has. Reconciling this
requires rederiving `MinAmplitudeMicrovolts` from real-data amplitude distributions, not
just retuning `ZThreshold` in isolation.

**`RefractoryMs` is entangled with `StepsPerMinuteToRejectChannel`, so they cannot be set
independently.** A refractory period caps the maximum achievable detection rate: a 1 s
refractory implies at most 60 events/min can ever be reported on one channel, so a
30/min bad-contact threshold with `RefractoryMs=20` (the shipped default, 3000/min
ceiling) is not testing the same thing as 30/min with `RefractoryMs=1000`. Any retuning
pass must move `RefractoryMs`, `ZThreshold`, `MinAmplitudeMicrovolts` and
`StepsPerMinuteToRejectChannel` together and re-validate the pair, not treat
`StepsPerMinuteToRejectChannel=30` as fixed while the other three move.

### Does `StepsPerMinuteToRejectChannel = 30` still sit in an empty gap on real data?

The design justified 30/min by arguing the two populations are separated by "more than
two orders of magnitude in per-channel rate," so the threshold "is not finely balanced."
That claim was checked directly against the shipped defaults' own output: for each file,
`stepRatePerMinute` was pulled from the saved `detect_*.mat`, restricted to channels with
at least one detected step, and its distribution computed (percentiles by linear
interpolation over the sorted per-channel rates —
`/home/nga4004/.claude/jobs/91008025/tmp/rate_gap_check.m`, reading the same Phase 1
artifacts as Sections 1-2, no new detection run):

| file | channels w/ steps | min | p10 | p25 | median | p75 | max | channels in [20,50]/min |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| m159 | 184 | 0.12 | 34.37 | 57.55 | 89.67 | 167.44 | 358.91 | 24 (13%) |
| m497 | 237 | 0.11 | 1.32 | 10.14 | 49.95 | 108.71 | 687.43 | 34 (14%) |
| m208 | 198 | 5.56 | 102.29 | 168.56 | 248.27 | 282.82 | 546.00 | 10 (5%) |
| m583 | 250 | 0.11 | 8.91 | 21.76 | 45.41 | 89.93 | 283.95 | 78 (31%) |
| m155 | 248 | 153.21 | 606.61 | 624.58 | 647.70 | 683.86 | 889.80 | 0 (0%) |

**No file shows an empty gap anywhere near 30/min.** On four of the five files (m159,
m497, m208, m583) the [20,50]/min band is densely populated — from 5% up to 31% of every
step-affected channel's rate falls inside it, and m583 alone puts nearly a third of its
affected channels there. m155 is the exception, but not in the design's favour: every
one of its 248 affected channels sits at 153/min or higher, so 30/min classifies literally
all of them as bad contacts with no gap to speak of either side of the threshold — there
is no low-rate population on this file at all at shipped defaults, which is the same
over-detection problem from Section 1 restated as a rate histogram rather than a raw
count. Across all five files, the two-population, wide-empty-gap structure the design
built `StepsPerMinuteToRejectChannel=30` around is not present in this data at shipped
detector settings. (This is consistent with, not separate from, Section 1: a detector
that is 10x too sensitive turns what should be a handful of genuinely fast-popping
channels into a continuum running from near-zero to several hundred per minute, and a
single fixed-rate cutoff drawn through a continuum has no non-arbitrary place to sit.)
This does not by itself say where 30 belongs once `ZThreshold` is corrected — that
distribution has not been recomputed at the corrected settings — but it does mean 30
cannot currently be defended as "sitting in an empty gap," and the coupling with
`RefractoryMs` above means it should be re-derived alongside `ZThreshold`, not audited on
its own.

**Caveat, stated plainly: this sweep was run on m159 only.** It was not the most extreme
file (m155 was, at 2190x) and it has not been re-run for m155, m208, m497, or m583.
`ZThreshold≈200` is a promising starting point for a rederivation pass, not a validated
fix — it is not yet known whether it also resolves m155's 248-channel result, and it
almost certainly needs `MinAmplitudeMicrovolts` retuned alongside it given the deficit
above.

## 4. The deferred guard-band ringing question: still open, not answered

Design decision 7 deferred the choice of whether to interpolate the few samples spanning
each step transition, pending a measurement of how much filter ringing survives at a
repaired step site after the pipeline's 1 Hz highpass. That measurement was taken:

| file | raw median (uV) | repaired median (uV) | step-free control median (uV) | residual ratio | repair time |
| --- | ---: | ---: | ---: | ---: | ---: |
| m155 | 225.5 | 223.9 | 185.2 | 1.21x | 54.6 s |
| m159 | 18.4 | 14.9 | 12.2 | 1.22x | 154.2 s |
| m208 | 20.5 | 16.6 | 15.0 | 1.10x | 151.7 s |
| m497 | 19.2 | 16.4 | 15.0 | 1.09x | 167.5 s |
| m583 | 18.5 | 15.4 | 13.3 | 1.16x | 182.8 s |

At face value a residual ratio of 1.09-1.22x looks reassuring — repaired step sites read
close to ordinary step-free data after the highpass. **Do not read it that way.** The
measurement is diluted by the detector's own false-positive rate from Section 1
(11.8x-2190.3x too many "step" sites), and most of the windows being averaged into
`ringingRawMedian`/`ringingRepairedMedian` are centred on ordinary data, not real pops —
that pulls both the raw and repaired medians toward the step-free control and compresses
the ratio toward 1.0 regardless of how well genuine pops are being handled.

**Limitation:** the "step-free control" window itself is not guaranteed step-free.
`validate_pop_repair.m` / `pop_repair_one.m` pick it as the same channel offset by
roughly half the recording length, with no check that the offset window avoids *other*
detected step sites on that channel. On a channel carrying hundreds of steps/min (routine
at these settings per Section 1) the control window can itself straddle a detected step.
This is a limitation of the measurement, not a caveat that weakens the "unresolved"
conclusion above — a control window contaminated by steps would itself show elevated
ringing, pulling `controlMedian` up and the ratio further toward 1.0, so if anything this
reinforces rather than undermines the case that the measured ratios cannot be trusted as
reassurance yet.

**m155 is the tell.** It has essentially one real bad channel and the rest of its
"1,390,837 steps" are near-certainly the detector's own noise (Section 1). Repair moved
its ringing median from 225.5 to 223.9 uV — a **0.7%** reduction — against a control
already at 185.2 uV. If repair were doing meaningful work at genuine pop sites on this
file, the reduction would look like m159/m208's ~19%, not 0.7%. The 0.7% figure is
measuring "what happens when you run a drift-detrending stage on 1.39 million site that
are mostly not pops," not "how much ringing survives at a genuine repaired pop."

**Committed answer: the guard-band ringing question cannot be answered from this run.**
The instrument used to select the measurement windows (the step detector) is not
trustworthy at its shipped defaults, so the windows being measured are not reliably real
pops. Recommendation: after Section 3's threshold rederivation lands (and is validated
across all five files, not just m159), rerun exactly this measurement
(`validate_pop_repair.m` / `pop_repair_one.m`, unchanged logic) at the corrected
defaults. If residual ratios stay in the 1.1-1.2x range on a detector whose false-positive
rate has been brought down to roughly match the survey, that is real evidence guard-band
interpolation is unnecessary. Until then, no claim either way is supportable — the
guard-band-left-as-recorded default should remain unless and until that rerun says
otherwise.

## 5. Repair runtime as a diagnostic: fast is suspicious, not reassuring

Ladder-rung distribution per file (`treatmentCounts`: `order3oa` = order-3 windowed
overlap-add, `order3/order2/order1` = single-window polynomial fits of that order,
`mean` = weighted-mean-only, `fallback` = the ladder stepped down after an `nt_detrend`
failure):

| file | sub-segments | order3oa | order3 | order2 | order1 | mean-only | fallback | mean-only % | median sub-seg | repair time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| m155 | 1,406,514 | 1,027 | 241 | 824 | 2,487 | 1,401,935 | 0 | 99.7% | 0.050 s | 54.6 s |
| m159 | 193,676 | 14,837 | 3,095 | 5,584 | 9,162 | 160,998 | 0 | 83.1% | 0.152 s | 154.2 s |
| m208 | 521,191 | 11,579 | 3,667 | 7,087 | 13,892 | 484,966 | 0 | 93.0% | 0.088 s | 151.7 s |
| m497 | 167,315 | 15,037 | 3,128 | 5,568 | 9,940 | 133,642 | 0 | 79.9% | 0.236 s | 167.5 s |
| m583 | 160,743 | 18,229 | 4,023 | 6,417 | 10,468 | 121,606 | 0 | 75.7% | 0.259 s | 182.8 s |

(`fractionSamplesMeanOnly` from `RELAX_repair_steps` — a sample-weighted fraction, not the
plain rung-count fraction shown here — reads slightly lower per file: 91.7% / 24.9% /
42.5% / 23.2% / 20.3% for m155/m159/m208/m497/m583 respectively. Both tell the same
story.)

Counter-intuitively, **m155 (1.39 million steps) repaired in 54.6 s, faster than every
other file**, and repair time does not track step count at all: m208 (505,514 steps, most
of the five) took 151.7 s, close to m497's 167.5 s (151,638 steps, fewest of the four
"normal" files) and less than m583's 182.8 s (145,066 steps, the fewest overall). The
driver is **treatment**, not count: more cuts produce shorter sub-segments, and a
sub-segment under 0.845 s degrades straight to weighted-mean subtraction (`nt_detrend` is
never called), which is far cheaper than an iterative robust polynomial fit. m155's median
sub-segment is 0.050 s — a fortieth of a second — so 99.7% of its sub-segments hit the
mean-only floor immediately; almost none of `RELAX_repair_steps`'s actual polynomial
machinery ever runs on this file.

This makes runtime (or equivalently `fractionSamplesMeanOnly`) a usable, already-recorded
diagnostic for exactly the failure mode in this report: **a suspiciously fast repair pass
combined with a high mean-only fraction indicates the detector over-fired and the stage
silently degenerated into piecewise-constant baseline chopping rather than the polynomial
detrend the design specifies.** No new instrumentation is needed — `PopFractionSamplesMeanOnly`
and the per-rung `treatmentCounts` columns specified in the design's "Metrics and
covariates" section already capture this per file.

## Conclusions and recommendations, in order

1. **The stage fails validation at its shipped defaults.** m155 alone is sufficient: the
   survey's 2-bad-channel example returns 248/257 channels flagged. `RELAX_cfg.Perform_pop_artifact_repair`
   was `= 1` in the shipped `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` before this report —
   i.e. real runs following that entry point unmodified were already executing the stage
   at exactly the defaults measured broken above. This report's commit flips that line to
   `= 0` with a comment pointing here; leave it off until thresholds are rederived and
   re-validated on all five files. `RELAX_Wrapper_slurm.m`'s own absent-means-off guard
   (`:255-257`) only protects a `RELAX_cfg` that never sets the field at all — it does not
   help a copy of the entry point that still carries the old `= 1` forward, so check any
   such copy explicitly rather than relying on the wrapper's default.
2. **The cause is `ZThreshold`, and it is fixable, pending broader validation.**
   `ZThreshold≈200` reproduces the survey on m159 on both step count and channel count;
   `MinAmplitudeMicrovolts` needs independent rederivation (the 92 uV probe showed an
   11x deficit, not a match); `RefractoryMs` and `StepsPerMinuteToRejectChannel` are
   coupled and must be retuned together, not independently. Rerun the Section 3 sweep on
   all five files, not just m159, before committing to new defaults.
3. **The guard-band ringing question (design decision 7) is unresolved.** The 1.09-1.22x
   residual ratios measured here are contaminated by the detector's own false-positive
   rate and cannot be trusted as evidence either for or against guard-band
   interpolation — m155's 0.7% repair effect on an essentially all-false-positive file is
   the proof. Rerun `validate_pop_repair.m` after Section 3's fix lands; only then does
   the residual-ratio number mean what decision 7 needs it to mean.
4. **Repair runtime and `fractionSamplesMeanOnly` are a live diagnostic for this exact
   failure mode**, already wired into the per-file covariates the design specifies — a
   fast repair pass with a high mean-only fraction is a red flag, not a clean bill of
   health.

**None of this is a defect in the code built in Tasks 1-5.** `RELAX_detect_steps` and
`RELAX_repair_steps` behave exactly as their design and tests specify — the unit tests
(including the spectral-fidelity test and the numerically-verified ladder boundaries)
pass, `RELAX_channel_budget` truncates correctly and warns as designed, and the repair
stage's degradation to mean-only on short sub-segments is documented, intended behaviour,
not a bug triggered here. The defect is in the **design spec's shipped default
thresholds**, which were derived from an exploratory survey (`pop_detect2.m` /
`pop_char.m`, not present in this repo, not recoverable during this task) whose detection
procedure evidently differs materially from the production detector's — most clearly in
effective sensitivity, since matching the survey's absolute-microvolt sense of "a pop"
required moving `ZThreshold` an order of magnitude from FieldTrip's convention-based
default.
