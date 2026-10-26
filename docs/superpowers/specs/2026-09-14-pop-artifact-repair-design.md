# Electrode-pop detection and repair for the RELAX-GEDAI SLURM pipeline

> **Superseded in part on 2026-09-17 and 2026-09-22** - see the addenda at the end: the 2026-09-17
> detector (a two-window level comparison with transition bracketing) is itself superseded by the
> 2026-09-22 median-slope gate; transitions are filled (linear chord + spherical-spline detail),
> common-mode events are flagged (excision off by default), and eye-movement clusters are skipped.

Date: 2026-09-14
Status: approved design, not yet implemented
Scope: `RELAX_Wrapper_slurm.m` and new functions it calls; `vendor/noisetools/`

## Problem

The input files
(`~/TMS_response_variability/store/mdd_dlpfc/preprocessed_Aug26/tms_artifact_removed/`) are
saline-base EGI GSN MR-compatible recordings. Most files contain **electrode pops**: large,
near-instantaneous shifts in the DC baseline of a single channel, caused by a momentary change in the
electrode-electrolyte contact. A pop is a *step*, not a transient - the baseline moves and stays
moved.

The pipeline currently has no stage that recognises a step. The 1 Hz acausal Butterworth highpass at
`RELAX_Wrapper_slurm.m:268` is the first stage that sees one, and a zero-phase filter converts a step
into a **bipolar excursion extending both forwards and backwards in time** for the length of the
filter's transient - on the order of 1-2 s in each direction at 1 Hz - with peak amplitude on the
order of the step itself. On this data that means a several-hundred-microvolt artifact smeared across
a multi-second window, an order of magnitude above resting EEG, straddling the acausal boundary so it
contaminates data *before* the pop occurred.

Everything downstream inherits the damage. Extreme-value rejection deletes the rung-out window rather
than the millisecond that was actually bad; MWF and ICA see a high-variance structured component that
is neither brain nor a stationary artifact; and because pops are channel-specific while the
pipeline's thresholds are computed from within-recording statistics, a poppy session gets cleaned
differently from a quiet one in a way that tracks that session's electrode quality.

### Measured properties of the input data

Five files surveyed with an exploratory median-filtered-derivative detector (scratchpad
`pop_detect2.m` / `pop_char.m`; 257 channels, 1000 Hz, 61 inter-boundary segments each).

| File | channels affected | confirmed steps | rate/min | median \|step\| | top-1 chan share |
| --- | --- | --- | --- | --- | --- |
| m159_dlpfc_day1 | 169/257 | 15045 | 1759 | 92 uV | 11% |
| m497_dlpfc_day1 | 203/257 | 8387 | 925 | 121 uV | 6% |
| m208_dlpfc_day1 | 251/257 | 5269 | 473 | 150 uV | 4% |
| m583_dlpfc_day1 | 215/257 | 2110 | 227 | 254 uV | 19% |
| m155_dlpfc_day1 | **2**/257 | 635 | 73 | 1365 uV | 98% |

The exploratory detector has **no rise-time criterion** and therefore also flags smooth ramps (e.g.
E166 in m583). The counts above are upper bounds on true step count; they are used here to establish
*structure*, not absolute prevalence. The production detector adds the rise-time criterion the
exploratory one lacked.

**The structure is two populations, and they need opposite treatment:**

- **Population A - chronically bad contacts.** A small number of channels produce hundreds to
  thousands of steps each. m155 is the pure case: 98% of all detected steps come from one channel
  (E91, 624 of 635), with median step size 1365 uV. A channel stepping at that rate has no stable
  baseline to restore. Repair is meaningless; the channel is a bad contact and belongs in the
  rejection list.
- **Population B - sporadic isolated pops.** Moderate (50-1500 uV), isolated to 1-3 channels, at
  **1.3 to 16.8 per minute** across the five files, on 1 to 54 channels per file. These are genuine
  single events on otherwise usable electrodes. This is the repair target, and it is the population
  that is currently being converted into multi-second filter ringing.

The `offset retained at +1s / +0.25s` ratio (median 0.91-1.96 across files, where 1.0 is a pure DC
step) confirms the events are sustained level changes rather than spikes.

### Signal model that justifies repair rather than interpolation

Two components sit in the same low-frequency band but differ in spatial structure:

- **Drift** - slow, largely common-mode, present on every channel.
- **Pops** - abrupt, channel-specific, uncorrelated across electrodes.

Repair exploits the difference: fit and remove a smooth trend *within each interval between steps*,
independently per channel. The step is never fitted - it is excluded by construction, because it sits
on a partition boundary rather than inside a fit window. Interpolation, by contrast, discards the
channel's real neural signal for the whole affected period and replaces it with a spatial estimate;
on a 257-channel net with pops on up to 54 channels per file that is a large amount of real data
thrown away for an artifact lasting milliseconds.

### Why per-sub-segment detrending subsumes step removal

The obvious two-stage design - robust detrend, then subtract the residual between-step offset
(`nt_destep`) - is redundant. Reading `nt_regw.m` (the local NoiseTools tree,
`~/TMS_response_variability/NoiseTools`):

```matlab
mn = y - nt_demean(y,w);                              % line 31: weighted mean
...
if nargout>1; z = nt_demean(x,w)*topcs*b; z = z+mn; end   % lines 43, 58: added back into the trend
```

The fitted trend `z` returned to `nt_detrend` **carries the weighted mean**. Therefore `x - z` has
zero weighted mean over the fit interval. Detrending each inter-step sub-segment independently
already removes that sub-segment's offset; a subsequent `nt_destep` pass would find nothing to do.
This holds in all three local variants of `nt_regw` examined. **No `nt_destep` in this design.**

### Measured frequency response of the detrend operator

The `n/(2T)` rule of thumb for a polynomial detrend's -3 dB point is wrong by up to 70%, so the ladder
below is built from measurement rather than the heuristic (scratchpad `detrend_response.m`,
`wsize_response.m`).

Single-window, measured -3 dB frequency (Hz):

| order \ T | 0.5 s | 1.0 s | 1.5 s | 2.0 s | 4.0 s | 8.0 s | ratio to n/(2T) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1.686 | 0.848 | 0.565 | 0.427 | 0.208 | 0.105 | **1.69** |
| 2 | 2.453 | 1.234 | 0.822 | 0.621 | 0.312 | 0.157 | **1.24** |
| 3 | 3.249 | 1.634 | 1.089 | 0.822 | 0.414 | 0.202 | **1.09** |
| 5 | 4.726 | 2.377 | 1.584 | 1.196 | 0.602 | 0.303 | 0.96 |
| 8 | 6.875 | 3.458 | 2.304 | 1.740 | 0.875 | 0.440 | 0.87 |

The ratio is constant across duration, so `f_3dB(n,T) = c(n) * n / (2T)` with `c = [1.69 1.24 1.09]`
for orders 1-3 is exact enough to derive thresholds from.

`nt_detrend`'s windowed overlap-add path (order 3, `wsize` in seconds):

| wsize | 1.0 | 1.5 | 2.0 | 2.5 | 3.0 | 4.0 |
| --- | --- | --- | --- | --- | --- | --- |
| f_3dB (Hz) | 1.767 | 1.132 | **0.841** | 0.673 | 0.538 | 0.431 |

A `[0.5, 1.0] Hz` bound on the effective cutoff pins `wsize` to roughly 1.7-3.0 s.

## Decisions taken

1. **Repair, not interpolation.** Per the signal model above.
2. **Detect with a median-filtered derivative and a robust z-threshold** - the FieldTrip
   `ft_artifact_jump` recipe (median filter order 9, z threshold 20), plus an absolute microvolt
   floor. Well established and trivially defensible to reviewers. Spatial model-mismatch (STAR-style)
   is the documented second option, deferrable behind the same interface.
3. **Run before any filtering**, at native 1000 Hz, between `RELAX_Wrapper_slurm.m:237` and `:239`.
   Repairing after the highpass is not possible: the ringing has already spread the artifact across
   seconds, and no step remains to find.
4. **Split detection from repair behind a stable `steps` struct.** Either half can be replaced
   without touching the other.
5. **Adaptive fit order per sub-segment, capped at 3**, with a maximum effective window so the
   cutoff never falls below 0.5 Hz and never rises above 1.0 Hz. Derived arithmetically from the
   measured `c(n)`, not hardcoded.
6. **Sub-segments shorter than the order-1 threshold get weighted-mean subtraction only.** A
   polynomial fit on a sliver eats real signal.
7. **Guard-band samples are left as recorded.** The few samples spanning the step transition are not
   interpolated in this version. The seam is built so interpolation can be dropped in later, and
   **residual ringing severity is to be measured and reported after implementation** so the decision
   can be made on evidence.
8. **Population A escalates inside the stage** - marked bad, repair skipped, removed immediately via
   `pop_select`, exactly as the impedance and bridge screens do.
9. **A shared running channel-rejection budget** across the pop and impedance screens, so the two
   together cannot exceed `MaxProportionOfElectrodesThatCanBeDeleted` of the original montage. PREP
   does not currently spend from this budget (see below).
10. **Reference pops measured, not specially handled.** Steps appearing simultaneously across a large
    fraction of channels are counted and reported; whether they need common-mode treatment is decided
    from that evidence later.
11. **Auxiliary electrodes untouched** - requires no code, they are split into `Non_cleaned_electrodes`
    at line 231, before the insertion point.
12. **Vendor the runtime closure of the local NoiseTools tree** (9 files), not the 138-file tree.
13. **`RELAX_Wrapper_slurm.m` only.** The sibling wrappers are not modified; the gap is documented in
    `CLAUDE.md` per the naming-convention rule.

## Design

### Vendored dependency

`vendor/noisetools/` receives the runtime closure of `nt_detrend`, copied from
`~/TMS_response_variability/NoiseTools`:

```
nt_detrend.m  nt_regw.m  nt_demean.m  nt_demean2.m
nt_fold.m  nt_unfold.m  nt_vecadd.m  nt_vecmult.m  nt_greetings.m
```

`nt_detrend_helper` is a local subfunction of `nt_detrend.m` (line 163), not a separate file.
`nt_detrend2`, `nt_read_data` and `nt_imagescc` appear only inside the `if 0` demo block at the
bottom of `nt_detrend.m` and are not needed.

`vendor/noisetools/VENDORING.md` follows `vendor/zapline-plus/VENDORING.md`: records the source tree
path, the copy date, a per-file line count and SHA-256 so later drift from the source is detectable,
the reason the local tree was chosen over the public mirror (the public `nt_regw` has its
`error('weights all zero')` commented out, silently returning `b=0` and propagating NaN; the local
tree raises), the closure analysis above, and the procedure for adding a further `nt_*` function.

`RELAX_add_vendored_noisetools()` mirrors `RELAX_add_vendored_zapline()`: prepends the directory to
the path, errors with a named identifier if the directory is missing, and is called from
`RELAX_SET_PARAMETERS_AND_RUN_slurm.m` after `eeglab()`. `nt_greetings` prints its citation banner
once per session at first call, not per invocation.

### `RELAX_detect_steps(EEG, 'Name', value, ...)` -> `[steps, info]`

Detection, per channel, on the raw unfiltered data:

1. Median-filter the channel at `MedianFilterOrder` (default 9) to suppress single-sample spikes
   without smearing a genuine step edge.
2. Take `abs(diff(...))` of the filtered trace.
3. Robust z-transform using **median and MAD computed per channel over the whole recording**, not
   per segment. An 8 s inter-boundary segment is too short for a stable scale estimate, and MAD is
   insensitive to the pops themselves.
4. Threshold at `ZThreshold` (default 20) **and** require the robust level change across the
   candidate to exceed `MinAmplitudeMicrovolts` (default 10). The absolute floor prevents a very
   quiet channel from generating detections from its own noise floor.
5. Require the level change to be **sustained** - the rise-time criterion the exploratory detector
   lacked - and enforce `RefractoryMs` (default 20) minimum spacing between steps on one channel.

`SustainMs` (default 100) sets the window either side used to measure the pre- and post-step level -
the level change must persist across it, which is what makes this a step detector rather than a spike
detector. `MinRiseFraction` (default 0.5) is the rise-time criterion itself: the fraction of the total
level change that must occur in the single largest derivative sample. A genuine pop delivers most of
its change in one sample; a ramp spreads it over thousands and fails.

`steps` is a struct of parallel arrays, one entry per detected step, sorted by channel then sample.
An empty detection returns `0x1` / `{}` fields, never a bare `[]`, so `numel(steps.sample)` is always
valid:

| field | type | meaning |
| --- | --- | --- |
| `.channel` | `[N x 1]` | row index into `EEG.data` |
| `.label` | `{N x 1}` | channel label, so the struct survives a `pop_select` |
| `.sample` | `[N x 1]` | first sample of the **new** level |
| `.amplitude` | `[N x 1]` | signed uV, robust post-level minus pre-level |
| `.z` | `[N x 1]` | detection statistic at that sample |

`info` leads with `.available` / `.status` to match `RELAX_bridged_channels` and
`RELAX_impedance_bad_channels`, then: `.params` (echo of every effective setting, so the saved struct
is self-describing), `.nSteps`, `.nChannelsAffected`, `.stepsPerChannel` `[nbchan x 1]`,
`.stepRatePerMinute` `[nbchan x 1]`, `.badChannels` + `.badChannelReason`, `.budgetExceeded` +
`.badChannelsBeforeBudget`, `.degenerateChannels`, `.nCommonModeSteps`, `.commonModeSamples`.

A **common-mode (reference) step** is defined for reporting as a 25 ms bin in which steps are detected
on more than 25% of surviving channels. That criterion is measurement-only in this version - such
steps are still partitioned and repaired per channel exactly like any other, per decision 10. The 25 ms
bin width matches the co-occurrence analysis in the survey above.

Population A criterion: `stepRatePerMinute > StepsPerMinuteToRejectChannel` (default 30). Candidates
are ordered worst-first by step rate before being passed to the budget.

Unlike the detector settings, this default is not borrowed from FieldTrip and is derived from the
survey above. The two populations are separated by more than two orders of magnitude in per-channel
rate: the clearest Population A channel (E91 in m155, 624 steps over ~8.5 min) runs at roughly
72 steps/min, while Population B channels average well under 1 step/min (m583: 156 sporadic steps
spread over 54 channels in ~9.3 min). A threshold of 30/min sits in the empty gap between them, so
the value is not finely balanced - but it is the one default that should be re-checked against the
Step 5 report on real files.

### `RELAX_repair_steps(EEG, steps, 'Name', value, ...)` -> `[EEG, info]`

**Per-channel partition.** For channel `c`:

```
cuts_c = unique( segment_starts  U  steps.sample(steps.channel == c) )
```

`segment_starts` is re-derived from the **current** `EEG` via `RELAX_get_segments`, per the design
rule in the boundary spec. Both kinds of cut are **hard**: no padding, no overlap, each sub-segment
fitted in complete isolation. That isolation is what makes the fit absorb the offset, per the
`nt_regw` analysis above.

`RELAX_apply_by_segments` is deliberately **not** reused. It asserts a single partition tiling the
whole record and shared by all channels, and it reflection-pads. Here the partition differs per
channel, and padding across a step would reintroduce exactly the discontinuity the partition exists
to isolate.

**Treatment ladder.** For a sub-segment of duration `T`, choose the largest order `n <= MaxPolyOrder`
satisfying `c(n) * n / (2T) <= MaxCutoffHz`, with `c = [1.69 1.24 1.09]` held as a constant in the
function. At the defaults (`MaxPolyOrder = 3`, `MaxCutoffHz = 1.0`, `MaxWindowSeconds = 2.0`) this
yields:

| sub-segment duration T | treatment | resulting cutoff |
| --- | --- | --- |
| T >= 2.000 s | `nt_detrend` order 3, `wsize = 2.0 s` (overlap-add) | 0.84 Hz |
| 1.635 <= T < 2.000 s | order 3, single window | 0.82-1.00 Hz |
| 1.240 <= T < 1.635 s | order 2, single window | 0.76-1.00 Hz |
| 0.845 <= T < 1.240 s | order 1, single window | 0.68-1.00 Hz |
| T < 0.845 s | weighted mean subtraction only | DC only |

`nt_detrend` takes `wsize` in **samples** (`nt_detrend.m:54`, `wsize = 2*floor(wsize/2)`), so
`MaxWindowSeconds` is converted with `round(MaxWindowSeconds * EEG.srate)` at the call site and
forced even. The ladder is continuous at T = 2.0 s (0.822 Hz single-window vs 0.841 Hz overlap-add). Generating it
from two config knobs rather than five hardcoded durations means the invariant "this stage never
attenuates above `MaxCutoffHz`" is stated once and enforced arithmetically; raising `MaxPolyOrder`
moves the thresholds automatically instead of leaving them stale.

**Hard vs soft subdivision.** Length-limited subdivision of a long sub-segment uses `nt_detrend`'s
`wsize` overlap-add (50% overlap, triangular weighting, continuous reconstructed trend) and **not**
further hard cuts. Hard-cutting a smooth stretch every 2 s would manufacture new step discontinuities
at 2 s intervals - precisely the input this stage exists to eliminate.

Channels with no detected steps are still partitioned by boundaries alone and still detrended. Every
channel therefore receives the same drift treatment, so the stage introduces no systematic difference
between repaired and unrepaired channels.

Fitting is per channel (`nt_detrend` receives a `[T x 1]` column) because the partition is
per channel. Initial weights `w0` are all ones; `nt_detrend`'s own iterative reweighting
(`DetrendThresh = 3` sd, `DetrendIterations = 3`) handles outliers. No prior artifact mask exists at
this point in the pipeline - extreme-value marking is stage 5, far downstream.

`info`: `.available`, `.status`, `.nChannelsRepaired`, `.nSubSegments`, `.subSegmentsPerChannel`,
`.treatmentCounts` (one count per rung plus `meanOnly`, `fallback`, `skipped`),
`.minSubSegmentDuration`, `.medianSubSegmentDuration`, `.lowFreqVarianceRemoved` `[nbchan x 1]`,
`.params`.

### `RELAX_channel_budget(EEG, op, ...)`

```matlab
EEG            = RELAX_channel_budget(EEG, 'init', maxProportion);
[fits, EEG]    = RELAX_channel_budget(EEG, 'spend', stageName, worstFirstLabels);
```

State on `EEG.RELAX.ChannelBudget`: `.nOriginalChannels`, `.maxRejectable`, `.nSpent`, `.spentBy`
(struct keyed by stage name), `.exhausted`. `spend` takes a **worst-first ordered** candidate list and
returns only what fits, so truncation always sacrifices the least-bad candidates.

Initialised at `RELAX_Wrapper_slurm.m:238`, immediately after `EEG.allchan = EEG.chanlocs;` - that
line records the reference montage the proportion is taken against, so the budget must be created
against the same snapshot.

Participating screens: **pop, impedance**. PREP does **not** currently spend from this budget -
`RELAX_Wrapper_slurm.m` has no `RELAX_channel_budget(...,'spend','PREP',...)` call (the only `'PREP'`
spend in the repo is a scenario in `tests/test_RELAX_channel_budget.m`). Consequence: if
`RELAX_cfg.Perform_PREP_bad_channel_rejection` is enabled (it is 0 by default in the SLURM entry
point), total deletions across all three screens can still exceed
`MaxProportionOfElectrodesThatCanBeDeleted`, since PREP applies only its own independent cap. Wiring
PREP into the shared budget is out of scope for this plan. Two further deliberate exclusions:

- The **bridge screen stays uncapped**, as it is today. A bridge is a physical short; a retained
  cluster member holds a smeared mixture while its `chanlocs` entry claims a precise location, so
  capping it would preserve channels known to be invalid.
- **`RELAX_excluding_channels_and_epoching`** keeps its existing independent later cap. Bringing it in
  means touching three sibling variants (`.m`, `_iTBS.m`, `_boundary.m`) and is a separate change.

`RELAX_impedance_bad_channels` gains an optional 4th argument `maxRejectableCount` which, when
supplied, overrides the internal `floor(maxProportion * nChannels)` at line 117. Its existing
worst-first sort, `capApplied` flag and warning are reused unchanged, and the 3-argument signature
keeps working for the sibling wrappers.

The comment at `RELAX_Wrapper_slurm.m:364` claiming the impedance cap is already shared with the
epoching stage is corrected: today those are two independent caps at the same proportion, not one
budget.

### Configuration

All new fields, with absent treated as **off** for `Perform_pop_artifact_repair` (new behaviour must
not switch itself on for other datasets) and absent treated as the stated default for the rest:

```matlab
RELAX_cfg.Perform_pop_artifact_repair             = 1;
% detector - FieldTrip ft_artifact_jump defaults
RELAX_cfg.PopDetect_MedianFilterOrder             = 9;
RELAX_cfg.PopDetect_ZThreshold                    = 20;
RELAX_cfg.PopDetect_MinAmplitudeMicrovolts        = 10;
RELAX_cfg.PopDetect_RefractoryMs                  = 20;
RELAX_cfg.PopDetect_SustainMs                     = 100;
RELAX_cfg.PopDetect_MinRiseFraction               = 0.5;
RELAX_cfg.PopDetect_StepsPerMinuteToRejectChannel = 30;
% repair
RELAX_cfg.PopRepair_MaxPolyOrder                  = 3;
RELAX_cfg.PopRepair_MaxCutoffHz                   = 1.0;
RELAX_cfg.PopRepair_MaxWindowSeconds              = 2.0;
RELAX_cfg.PopRepair_DetrendThresh                 = 3;
RELAX_cfg.PopRepair_DetrendIterations             = 3;
```

`RELAX_check_boundary_support` gains **no** new entry: the stage is boundary-aware by construction,
because every partition begins from `RELAX_get_segments`.

### Wrapper wiring

Inserted between `RELAX_Wrapper_slurm.m:237` (`EEG.allchan = EEG.chanlocs;`) and `:239` (the notch
block). The ordering matters: `allchan` must record the full montage *before* this screen removes
anything, as it already does for the impedance and bridge screens.

The block follows the bridge screen's shape verbatim - output field initialised before the `if`,
`warning` on `~available`, `fprintf` for the no-op case, `pop_select` + `eeg_checkset` only when there
is something to remove, `else` branch printing the skip message naming the config field.

## Error handling

| hazard | handling |
| --- | --- |
| `nt_detrend` windowed path contains `if ~sum(w(start:stop)); error('!'); end` - a bare error with no identifier | each sub-segment fit is wrapped; on failure step **down** the ladder (order 3 -> 1 -> mean-only) and increment `treatmentCounts.fallback`. One pathological window must never kill a 257-channel SLURM task hours into a run. |
| flat or dead channel (MAD = 0) | detector emits no steps, records the label in `info.degenerateChannels`; repair applies mean-only. No `0/0` z-scores. |
| NaN samples inside a sub-segment | mean-only, counted. `nt_regw` would otherwise propagate NaN through the entire fit. |
| step within `RefractoryMs` of a boundary | merged into the boundary cut rather than creating a 1-sample sub-segment. |
| `vendor/noisetools/` missing | `RELAX_add_vendored_noisetools()` errors with a named identifier, mirroring the zapline helper. |
| budget exhausted mid-stage | excess candidates recorded in `info.badChannelsBeforeBudget`, `info.budgetExceeded = true`, warning naming the file. Nothing silently disappears. |

## Metrics and covariates

`CLAUDE.md`'s "Saving covariates for downstream ML classification tasks" applies directly: this
stage's aggressiveness scales with how poppy a session was, so a model predicting treatment response
from session-to-session EEG change could learn the cleaning instead of the neural signal. The
aggressiveness must therefore be recorded per file, not just the outcome.

```matlab
EEG.RELAXProcessingExtremeRejections.PopBasedChannelToReject = {...};  % matches the sibling fields
EEG.RELAX.PopScreen = detectInfo;
EEG.RELAX.PopRepair = repairInfo;
```

Folded into the per-participant aggregate row and CSV: `PopNSteps`, `PopNChannelsAffected`,
`PopMedianStepAmplitude`, `PopNChannelsRejected`, `PopMedianSubSegmentDuration`,
`PopFractionSamplesMeanOnly`, `PopMedianLowFreqVarianceRemoved`, `PopNCommonModeSteps`, and one column
per ladder rung from `treatmentCounts`.

`PopFractionSamplesMeanOnly` is the covariate to watch hardest: it measures directly how much of the
record was too fragmented to detrend properly, and it will correlate with session noise.

Initialising `PopBasedChannelToReject` to `{}` before the `if` keeps the `struct2table` aggregation
shape identical whether the stage runs, finds nothing, or is switched off - the same property the
boundary spec relied on for `PREPBasedChannelToReject`.

## Implementation steps

### Step 0 - Vendor NoiseTools
Copy the 9-file closure, write `VENDORING.md` with checksums, add
`RELAX_add_vendored_noisetools.m`, call it from `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` after
`eeglab()`. Verify `nt_detrend` resolves to the vendored copy with `which -all`.

### Step 1 - `RELAX_channel_budget` and the impedance change
New helper plus the optional 4th argument on `RELAX_impedance_bad_channels`. Tests first. The
impedance screen's behaviour must be **byte-identical** when the new argument is absent.

### Step 2 - `RELAX_detect_steps`
Detector plus its tests. Built and validated against synthetic signals before any real file.

### Step 3 - `RELAX_repair_steps`
Partition construction, ladder, fallback handling, plus tests - including the spectral-fidelity test.

### Step 4 - Wrapper wiring and config
Insert the block, add the `RELAX_cfg` fields, wire `EEG.RELAX.*` and the aggregate-row columns,
initialise the budget at line 238, correct the stale comment at line 364.

### Step 5 - Real-file validation and the deferred ringing report
Run on the five surveyed files. Report: step counts before/after, channels escalated, ladder rung
distribution, `PopFractionSamplesMeanOnly`, and **the residual ringing severity at repaired step
sites after the 1 Hz highpass** - the measurement promised by decision 7, which determines whether
guard-band interpolation is added.

### Step 6 - Documentation
`CLAUDE.md`: new stage in the pipeline architecture list **between stage 1 (load & prep) and
stage 2 (notch filter)** - it runs before any filtering, so it is not adjacent to the other
channel-rejection screens at stage 4 even though it can reject channels. The `vendor/noisetools/`
entry alongside `vendor/zapline-plus/`, the new `RELAX_cfg` fields, the shared-budget behaviour, and
an explicit note that the sibling wrappers do **not** have this stage.

## Testing

Following the `tests/test_<function>.m` convention; run with
`matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"`.

**`test_RELAX_detect_steps.m`**
- recovers a known step's sample index and signed amplitude on a synthetic channel
- **no detection on a smooth ramp** - the exact false-positive mode the exploratory detector had
- no detection on pink noise at the default z threshold
- `RefractoryMs` enforced: two steps 5 ms apart yield one detection
- flat channel (MAD = 0) yields no steps and is listed in `degenerateChannels`, no NaN
- step placed one sample after a boundary is merged into the boundary cut
- `MinAmplitudeMicrovolts` floor suppresses a high-z but tiny step on a very quiet channel

**`test_RELAX_repair_steps.m`**
- a 500 uV step is reduced below tolerance
- **a 2 Hz sinusoid survives with < 5% attenuation while 0.3 Hz is strongly attenuated** - this is
  the test that certifies the design claim, not just the implementation; if the measured `c(n)`
  constants were wrong or the ladder drifted out of sync with `MaxCutoffHz`, nothing else notices
- each duration threshold selects the expected ladder rung (boundary cases at 0.845, 1.240, 1.635,
  2.000 s)
- a channel with zero steps is still detrended by boundary partition alone
- the all-zero-weight path is forced and lands on mean-only with `treatmentCounts.fallback == 1`
- a sub-segment containing NaN takes the mean-only path without propagating NaN

**`test_RELAX_repair_steps_boundaries.m`**
- synthetic `EEG` with boundary events: no fit window spans a boundary, and data either side is
  bit-identical to processing that side alone as a standalone record

**`test_RELAX_channel_budget.m`**
- two screens exhausting one budget; worst-first truncation; `exhausted` flag; `spentBy` attribution
- `RELAX_impedance_bad_channels` with 3 arguments is unchanged from current behaviour

## Files touched

New:
- `RELAX_detect_steps.m`
- `RELAX_repair_steps.m`
- `RELAX_channel_budget.m`
- `RELAX_add_vendored_noisetools.m`
- `vendor/noisetools/` (9 files + `VENDORING.md` + `LICENSE` if the source carries one)
- `tests/test_RELAX_detect_steps.m`
- `tests/test_RELAX_repair_steps.m`
- `tests/test_RELAX_repair_steps_boundaries.m`
- `tests/test_RELAX_channel_budget.m`

Modified:
- `RELAX_Wrapper_slurm.m` (stage block at 238, budget init, impedance call, comment at 364, aggregate row)
- `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` (new `RELAX_cfg` fields, vendored-path call)
- `RELAX_impedance_bad_channels.m` (optional 4th argument)
- `CLAUDE.md`

Not modified, deliberately: `RELAX_Wrapper.m`, `RELAX_Wrapper_bart.m`,
`RELAX_excluding_channels_and_epoching*.m`, `RELAX_bridged_channels.m`.

---

## Addendum 2026-09-17 - two-window detector, transition fill, common-mode excision

The survey of the GSN-256 saline-net TMS-EEG data (`variance_ratio_survey/`, characterization in
`variance_ratio_survey/pop_characterization/`) showed that the artifact this pipeline meets is
not the single-sample jump the design above assumed:

- 10-90 % rise median 9.5 ms (IQR 3-20 ms, max ~112 ms); staircases (2-3 rungs over ~80 ms),
  sigmoids (~150 ms), exponential approaches (tau ~40 ms), occasional ~300 ms pre-ramps.
- |A| 50 uV - 13 mV (median ~190 uV); telegraph toggling (plateaus 0.1-5 s) and 100-200 ms square
  pulses on face/neck-ring channels.
- Four populations: (A) single-channel large bistable pops; (B) local same-sign clusters of 2-12
  neighbours; (C) rim-wide mixed-sign shifts on 30-250 channels (movement); (D) uniform small
  all-channel shifts (reference; removed by the average re-reference).

A slow rise spreads the derivative over many samples so it never reaches z > 20; staircases fail
the sustain / variance-ratio windows; and the single-sample cut leaves the ramp inside a
sub-segment, which the later 1 Hz zero-phase highpass turns into a multi-second bipolar
excursion. Real-data validation of the shipped stage (`docs/superpowers/reports/
2026-09-14-pop-repair-validation.md`) failed. Both halves were therefore rewritten behind the same
interfaces; this addendum records what changed and which decisions above it supersedes.

### Detection - `RELAX_detect_steps.m` (rewrite; supersedes decision 2)

Two level windows of `WindowMs` (150) separated by a `TransitionMs` (100) allowance:
`D(t) = median(post) - median(pre)`, `S(t) = hypot(robust sd pre, robust sd post)`, `z = D / S`.
Three gates: local `|z| > ZThreshold` (4); an absolute floor `|D| >= MinAmplitudeMicrovolts` (50);
and a **channel-wide gate** `|D| / (1.4826 MAD of D over the channel's quiet samples) >
ChannelZThreshold` (4). The channel gate is the finding that made the two-window statistic
usable: on a channel carrying 1/f^2 drift at 150-200 uV rms, the local z alone produces ~80 false
candidates per minute because the within-window robust sd underestimates the drift's
level-to-level variation; the channel-wide scale of D absorbs it (6/6 synthetic pops found, 0
false). Same-sign super-threshold runs closer than `RefractoryMs` (20) are one candidate, located
at the run midpoint (the |z| plateau is flat over ~W + 2g samples). The transition is then
bracketed on a `SmoothMs` (11) moving median: onset = last sample before the midpoint within a
settle band of the pre level, offset = first sample after within the band of the post level, band
= max(`BoundaryFraction` (0.10) |A|, `BoundaryNoiseMultiple` (2) x robust sd), levels re-measured
from windows adjacent to the bracket, `RefineIterations` (3) passes. A level that never settles
within `MaxTransitionMs` (300) or the segment edge is rejected (`info.nRejected.unbounded`).

Steps within `CommonModeWindowMs` (50) of each other form an event (`steps.eventId`,
`info.events`); an event on more than `CommonModeChannelFraction` (0.08) of channels (and at
least 2) is `commonMode`. Bad contacts are ranked worst-first by the binding one of three
criteria over the channel's own steps - neither common-mode nor, once the wrapper has flagged
them, eye-movement (the detector ranks over `~commonMode`; the wrapper re-ranks with
`RELAX_rank_pop_channels(steps, info)` after `RELAX_flag_eye_steps`, see the next section):
steps per minute (`StepsPerMinuteToRejectChannel`, now 6 -
each fill perturbs ~1 s of the channel after the highpass), fraction of samples inside a
transition (`MaxTransitionFractionToRejectChannel` 0.02), and fraction of samples where the local
statistic is over threshold (`MaxUnstableFractionToRejectChannel` 0.25 - fast telegraph noise no
single step describes). `MedianFilterOrder`, `SustainMs` and `MaxVarianceRatio` are gone; the
wrapper errors if a config still sets them.

Documented limits: pulses shorter than ~W/2 (75 ms) are not detected (left to extreme-period
marking); pulses >= ~75 ms become two opposite steps with a demeaned plateau; a 300 ms / 200 uV
ramp sits at the floor; the first and last ~W + g samples of every inter-boundary segment are a
dead zone for the *locator* (a bracket can still walk into it).

### Eye movements - `RELAX_pop_eye_mask.m`, `RELAX_flag_eye_steps.m` (new)

A horizontal saccade is a 20-60 ms same-sign level shift on a handful of peri-ocular channels -
population B to the detector. Repairing it would remove EOG structure the later ICA relies on.
The wrapper builds an eye-movement mask on a temporary copy with RELAX's own detectors
(`RELAX_blinks_IQR_method` after an average re-reference; `RELAX_horizontaleye` on the first
available `HEOGLeftpattern` / `HEOGRightpattern` pair band-passed 1-25 Hz), and
`RELAX_flag_eye_steps` marks every member of an event that spans >= 2 channels, includes an
eye electrode (`BlinkElectrodes` plus the HEOG pair the mask actually resolved) and overlaps the
mask. Flagged steps are neither filled nor excised - eye movement takes precedence over common
mode, because a rim-wide shift inside an eye period is EOG, not head movement. Any failure of
the mask degrades to "nothing flagged" with a warning.

**Tightened 2026-09-18** after an audit of the five validation files' saved step tables
(`docs/superpowers/reports/2026-09-18-eye-flag-audit.md`). The
flag rate of eligible events was 2-5x the mask fraction overall, so the mask does discriminate -
but the enrichment collapsed to 1.2-1.7x chance for 2-5-channel events on the files whose blink
mask covers 25-30 % (m497, m583), and on m155 every eye-electrode step (4 of 4, all 1.5-2.6 mV
in 3-channel clusters - contact pops by any definition) was flagged. Flagged eye-electrode
steps sit at median |A| 130-180 uV; only 0-1.4 % exceed 500 uV. The rule is now:
(2) at least *two* eye electrodes in the event, or one when the event is common-mode (an eye
electrode is incidental in a 30+-channel shift, and the point of flagging those is to keep the
EOG in place rather than excise it); (3) the mid-level crossing `steps.sample` of an eye
electrode - not any part of any member's `[onset, offset]` - inside the mask; (4) no member
over `PopEye_MaxAmplitudeMicrovolts` (500). And the HEOG detector's contribution to the mask
is opt-in (`PopEye_IncludeHEOG`, default 0): it compares two channels for opposite-polarity
excursions, and after the zero-phase 1-25 Hz bandpass a DC step is a bipolar excursion, so two
same-sign pops a few tens of ms apart on the pair satisfy the test and mark
`+-HorizontalEyeMovementFocus` around themselves. The pair is still resolved and still counts
as eye electrodes for rule (2) and for `RELAX_protect_eye_channels`. Eligible events rejected
by rule (4) alone are counted in `PopNEyeEventsOverAmplitude`.

**Order (changed 2026-09-17 after the first real-data run):** the mask and the flag run
*before* bad-contact rejection, on the still-complete montage, and the bad-contact ranking is
then recomputed without the eye-movement steps. Ranking first rejected E10, E46, E18, E37, E2,
E11, E38, E50 and E59 on m159 - a saccade is a genuine DC step on the periocular ring, so at a
few saccades per minute those electrodes exceed 6 steps/min - which would have removed the
blink and HEOG electrodes every later RELAX eye stage depends on, and would have removed them
before the eye mask itself could be built. The sample-level instability criterion is unaffected
by the re-ranking.

**Found by the first production batch (2026-09-17/18, array 19721732; report §6) and fixed:**

- *Eye electrodes exempt from rejection* (`RELAX_protect_eye_channels.m`). Where the HEOG
  detector marks almost nothing (m407: 1.1 %) saccade clusters outside blink periods stay
  unflagged and the rate criterion removed all 14 periocular electrodes; RELAX's blink stage
  then cannot run. An eye electrode (`BlinkElectrodes` ∪ HEOG pair) over a limit is now kept,
  its steps left alone (flagged like eye movement, so the interpolation ceiling cannot be
  exceeded on it), and counted in `PopNEyeElectrodesProtected`. Supersedes "rejected via
  pop_select" for those channels.
- *Per-segment degenerate detection, rejected as `flat`.* A railed channel (m544 day 2, E154:
  0 / +333333 µV, constant within every segment, 0 steps) is not globally constant, so the old
  degenerate test missed it; the per-segment robust demean then made it exactly flat and
  `pop_jointprob` crashed. `RELAX_detect_steps` now calls a channel degenerate when it is
  constant within every usable segment; `RELAX_rank_pop_channels` lists such channels first
  (reason `flat`), and `detrendSubSegment` returns a constant sub-segment unchanged
  (`treatmentCounts.constant`) so no channel can become globally flat.
- *Excisions never leave a sub-second remnant* (`MinSegmentSeconds`, 1.0). One cut on m292
  day 3 landed a few samples from a TMS boundary; after resampling the remnant was ≤ 1 sample
  and zapline's `pwelch_segmented`, which sized its window from the shortest segment, errored.
  A cut is now extended to the segment edge, or merged with the next cut, when the remnant
  would be shorter than `MinSegmentSeconds`; and `pwelch_segmented` sizes its window from the
  segments of at least 1 s (see `vendor/zapline-plus/VENDORING.md`).

### Repair - `RELAX_repair_steps.m` (rewrite; supersedes decision 7, and decision 1 for the
transition samples only)

The `nt_detrend` ladder (orders, `C = [1.69 1.24 1.09]`, overlap-add gating, robust demean,
fallbacks) is unchanged and still does the level and drift removal; what changed is the
partition and what happens to the transition:

1. Each channel is partitioned at {segment starts} U {its own `[onset, offset]` intervals} U
   {global excision intervals}; a legacy `steps` struct with only `.sample` reproduces the old
   single-sample partition bit-for-bit.
2. Every sub-segment is detrended in isolation, exactly as before.
3. **Each transition is filled** (decision 7 reversed, on evidence): a linear chord between the
   means of the `AnchorMs` (3) samples on either side, plus - `FillMethod 'spatial'`, the default -
   the spherical-spline estimate of the neighbours' anchor-relative detail
   (`RELAX_spline_weights`, numerically identical to EEGLAB `eeg_interp 'spherical'`). The basis
   excludes the channel itself, channels with their own interval within the anchor windows,
   channels without coordinates, and `BasisExcludeLabels` (the wrapper passes degenerate channels
   and contacts over `MaxImpedanceKOhm`); fewer than `MinBasisChannels` (16) falls back to the
   plain chord. Downstream analysis (active pipeline: 1-200 Hz order-4 Butterworth per segment,
   500 Hz, ZaplinePlus, `_boundary` epoching, average reference over ~203 channels, PICARD with
   PCA ~93): the highpass residual of a chord fill is 2.22 x the chord-deviation area (<= 0.02 uV
   at 10 ms, 2-5 uV at 100 ms, ~60 uV only for a 300 ms fill on a 150 uV blink, which the eye
   skip removes); the spatial detail cuts it 2-4 x for long fills; the ICA weight perturbation is
   ~f x kappa (<= 1 % at a 2 % per-channel fill fraction), below the sampling floor.
4. **Common-mode events are excised** (decision 10 resolved): the union of the event's
   intervals, clipped to segments, is deleted from all channels with `eeg_eegrej`, which inserts
   a `boundary` of that duration; events longer than `MaxExciseMs` (500) are left to
   extreme-period marking. Excising *every* pop was rejected because each new boundary shortens
   ZaplinePlus's Welch window file-wide (`pwelch_segmented`) and orphans up to 1 s of data on
   both sides in the `_boundary` epoching (2-10 % of a poppy file); for the handful of rim-wide
   events per file the cost is acceptable and no per-channel fill could be trusted.

`steps.sample` / `.onset` / `.offset` are on the pre-excision timeline; the wrapper stores them
with the excised list on `EEG.RELAX.PopScreen` / `EEG.RELAX.PopRepair`.

### Covariates added

`PopNEvents PopMedianRiseMs PopNChannelsRejectedByRate/ByFraction/ByUnstable
PopNCommonModeEvents PopNEyeMovementStepsSkipped PopEyeMaskFraction PopNIntervalsFilled
PopFractionSamplesInterpolated PopMaxChannelFractionInterpolated PopMedianTransitionMs
PopNFillsSpatial PopNFillsLinear PopNCommonModeExcised PopExcisedSeconds PopNCommonModeTooLong`,
alongside every covariate the first version wrote.

### Validation

`validate_pop_repair.m` (sbatch `run_validate_pop_repair.sh`) reproduces the stage on the
production 204-channel montage for the five surveyed files and reports the detection load, the
post-1 Hz-highpass residual at every filled transition for raw / legacy (commit f6649a8 at the
new sites) / linear / spatial versions against step-free controls, a leave-one-out test of the
spatial estimate on step-free windows, and an ICA agreement check (PICARD, matched ICs, ICLabel
classes, and a detrend-only arm plus a raw-vs-raw floor) raw vs repaired on one file. The stage
was switched **on** in the SLURM entry point on 2026-09-17 on the strength of that report
(`docs/superpowers/reports/2026-09-17-pop-repair-validation.md`).

## Addendum 2026-09-18 - screen order, exclusion, per-file isolation, blink fallback

Applied before the second production batch (all 26 participants, `OverwriteExistingCleanedFiles=1`),
in response to the first batch's failures (`docs/superpowers/reports/2026-09-17-pop-repair-validation.md` §6).

### Stage order (supersedes "Wrapper wiring" above)

The impedance screen and the bridge screen now run **before** the pop stage in
`RELAX_Wrapper_slurm.m`, directly after `EEG.allchan = EEG.chanlocs` and the budget init, on the
raw native-rate data. Budget order is impedance → bridge (does not spend) → pop. Rationale:

- The pop detector's channel-wide MAD gate, the common-mode channel count, the eye-mask
  blink median and the spherical-spline fill basis are all whole-montage statistics; a
  contact already known to be bad (impedance over `MaxImpedanceKOhm`, or bridged to its
  neighbour) should not vote in them. `PopBasisExclude` still covers the over-limit contacts
  the impedance screen had to leave in place once the budget was spent.
- Neither screen needs filtered data. The impedance screen reads `chanlocs.impedance`
  metadata. `eBridge`'s electrical distance is the variance of the pairwise difference per
  epoch, so a DC step or slow drift common to a bridged pair cancels rather than inflating
  it, and a step on one member only raises the ED (looks *less* bridged) - moving it ahead
  of the notch/bandpass cannot create false bridges.
- A recording that is excluded for excessive bridging (below) is now excluded before the
  most expensive stage rather than after it.

PREP `findNoisyChannels` (off in the shipped config) stays at its original position after the
filters.

### Excessive bridging excludes the recording

`RELAX_bridged_channels` already reports `excessiveBridging` when more than
`BridgeMaxProportionChannels` of the montage is bridged, but the wrapper used to warn and
strip the channels anyway (m396: 177 of 193 removed, then the blink method errored). It now
raises `error('RELAX:excludeRecording', ...)`. A net that bridged is an acquisition failure,
not something to clean around.

### Per-file `try/catch`

The body of the file loop is wrapped in `try ... catch fileErr`. On any error the wrapper
prints the extended report to stderr, warns `RELAX:fileNotProcessed`, appends a row
(`File, Outcome, Identifier, Message, When`; `Outcome` = `excluded` for
`RELAX:excludeRecording`, else `error`) to
`RELAXProcessed/PerParticipantStats/RELAX_unprocessed_files_<participant>.mat/.csv`, and
`continue`s to the participant's next file. A file that later succeeds has its row removed
before processing starts. The post-loop aggregation tolerates a participant whose every file
failed (empty aggregate tables, `RELAX_issues_to_check` column guards). One failing file no
longer ends the SLURM array task, and every file missing from `Cleaned_Data/` has a row saying
why.

### Blink method fallback (`RELAX_blinks_IQR_method.m`)

If fewer than two of `RELAX_cfg.BlinkElectrodes` survive the channel screens (m160 day 2:
the native extreme-value screen removed 6 of 8), the method adds the row one higher on the
GSN-256 net (`E12 E20 E27 E34 E21`) and retries; if that still leaves fewer than two it
takes the existing no-blink path (all-zero mask, `IQRmethodDetectedBlinks = 0`) instead of
erroring. Covariates on `RELAXProcessingExtremeRejections`: `BlinkElectrodesUsed` and
`BlinkElectrodeFallback` (0 primary set, 1 backup row added, 2 no blink marking).
Test: `tests/test_RELAX_blinks_IQR_fallback.m`.

## Addendum 2026-09-22 - median-slope gate

Supersedes this addendum's own 2026-09-17 "Detection" section (`RELAX_detect_steps.m`'s gating
statistic and the refinement default; repair, eye flagging and event grouping are unchanged). Full
design: `docs/superpowers/specs/2026-09-22-median-slope-gate-design.md`. Thresholds chosen from a
blind-scored survey on the four Sep26 files:
`docs/superpowers/reports/2026-09-23-slope-gate-survey.md`.

### What changed

- **Level-difference gate -> median-slope gate.** The 2026-09-17 detector compared medians of two
  `WindowMs` windows separated by a `TransitionMs` allowance and gated on how *far* the level moved.
  The new detector reads the same centred median filter across `TransitionMs` alone (an exact lag in
  samples, `Dslope(t) = mc(t + gPost) - mc(t - gPre)`) and gates on how *fast* it moved: excess over
  a `SlopeBaselineMs`-long local baseline slope, in local robust spread. A level-difference floor
  (`MinAmplitudeMicrovolts`) and a channel-wide gate on the same slope statistic (`ChannelZThreshold`)
  remain. This is what stops slow 1/f drift and slow waves, which reach the same *level* change as a
  pop given enough time, from being gated as if they moved as *fast* as one.
- **Sample-level instability channel criterion retired.** Its statistic
  (`info.unstableFraction`, samples passing the old local gate) has no counterpart under the new
  gate; `MaxUnstableFractionToRejectChannel` and the `PopNChannelsRejectedByUnstable` covariate are
  gone. `RELAX_rank_pop_channels` now scores channels on two criteria only: steps/min and fraction of
  samples in transition. A run against an existing `RELAXProcessed/` needs a fresh output directory,
  since the aggregate tables' columns changed.
- **Refinement is a single bracketing pass** (user decision 2026-09-25, after the Task 5 repair
  comparison in the survey report's section 11b): from the run's strongest sample, the step placed
  at the mid-level crossing and `[onset, offset]` where the trace settles within a band of the
  pre/post levels. A walk that never settles falls back to the run's own end rather than
  rejecting; only the amplitude-floor and segment-edge rejections remain. `RefineTransition` is a logical: `true`
  (default) is this pass, `false` reports the gate's own run.
- **Drift screens and iterative bracketing removed** (Task 4D, 2026-09-25, user). The plateau- and
  hold-drift screens (`MaxPlateauDriftFraction`, `HoldMs`, `MaxHoldDriftFraction`), the iterative
  bracketing with its `MaxTransitionMs` / `RefineIterations` rejections, and the `plateauDrift`,
  `holdDrift`, `rawZ` and `varianceRatio` report columns are gone. `info.nRejected` counts only
  `amplitude` and `edge`. The earlier sections of this spec describing them are historical.
- **Common-mode excision stays off by default** (`PopRepair_ExciseCommonMode=0`, unchanged from
  2026-09-17). Looking at the chosen point's common-mode events, the user sorted them into eye
  artifacts and possible reference-channel pops; the latter belong to the reference-pop pass
  (`Perform_reference_pop_repair`), not to excision.

### Shipped thresholds (`RELAX_SET_PARAMETERS_AND_RUN_slurm.m`)

| Parameter | Value |
| --- | --- |
| `PopDetect_WindowMs` | 151 |
| `PopDetect_TransitionMs` | 4 |
| `PopDetect_ZThreshold` | 20 |
| `PopDetect_ChannelZThreshold` | 35 |
| `PopDetect_MinAmplitudeMicrovolts` | 15 |
| `PopDetect_SlopeBaselineMs` | 1000 |
| `PopDetect_SlopeNoiseFloorMicrovolts` | 0.1 |
| `PopDetect_SlopeSpreadQuantile` | 0.85 |
| `PopDetect_RefineTransition` | `true` (single pass) |
| `PopDetect_StepsPerMinuteToRejectChannel` | 6 |
| `PopDetect_MaxTransitionFractionToRejectChannel` | 0.02 |
| `PopRepair_ExciseCommonMode` | `false` |
