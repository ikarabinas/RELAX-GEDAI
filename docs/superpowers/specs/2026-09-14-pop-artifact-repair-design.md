# Electrode-pop detection and repair for the RELAX-GEDAI SLURM pipeline

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
9. **A shared running channel-rejection budget** across the pop, impedance and PREP screens, so the
   three together cannot exceed `MaxProportionOfElectrodesThatCanBeDeleted` of the original montage.
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

Participating screens: **pop, impedance, PREP**. Two deliberate exclusions:

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
