# A median-slope gate for electrode-pop detection

Date: 2026-09-22
Status: implemented 2026-09-23/25, with changes. What shipped differs from this design: the
spread is a 0.85-quantile scale rather than the MAD, the baseline is 1000 ms rather than 2000 ms,
the local z is |Dslope - b| / s, and the plateau/hold drift screens were removed (Task 4D). For
the shipped behaviour read `docs/superpowers/plans/2026-09-22-median-slope-gate.md` (2026-09-23
revision onward) and `docs/superpowers/reports/2026-09-23-slope-gate-survey.md`. The safety
argument under "A wide centred span" assumed the MAD's 50 % breakdown. The 0.85 quantile breaks
down at 15 %, so a rise filling more than ~150 ms of the 1000 ms baseline masks itself (final
review, 2026-09-25; the user accepted the limit).
Scope: `RELAX_detect_steps.m`, `RELAX_rank_pop_channels.m`, and the plumbing that configures them
(`RELAX_Wrapper_slurm.m`, `RELAX_SET_PARAMETERS_AND_RUN_slurm.m`, `validate_pop_repair.m`); a new
`slope_gate_survey/`. Supersedes the detector's gating statistic as described in the 2026-09-17
addendum to `2026-09-14-pop-artifact-repair-design.md`; the bracketing, repair, ranking and
eye-flagging stages are unchanged.

## Problem

`RELAX_detect_steps` gates candidates on the **size** of a level change. For every sample `t` it
forms `D(t) = median(post window) - median(pre window)` and requires

    |D| / hypot(robust sd pre, robust sd post) > ZThreshold        (local)
    |D| / (1.4826 * MAD of D over quiet samples) > ChannelZThreshold  (channel)
    |D| >= MinAmplitudeMicrovolts                                  (floor)

Nothing in that gate asks how *fast* the level moved. A slow excursion - sweat drift, an electrode
drying out, a large slow wave - that happens to be big enough clears all three tests, and is then
handed to the bracketing stage, which sometimes succeeds in bracketing it. The downstream defences
(`unbounded` rejection when the smoothed trace never settles inside `MaxTransitionMs`, and
`MaxPlateauDriftFraction = 0.3`) catch some of these but not all, because a drift that flattens
briefly on both sides of a steep stretch looks like a settled step at every test the detector applies.

The cost of a false positive is not neutral. `RELAX_repair_steps` detrends the sub-segments either
side of the interval and fills the interval with a spline-plus-chord estimate, so a bracketed drift
has genuine low-frequency signal removed and replaced with an interpolation.

**Goal: reject drift and slow waves at the gate, without losing the pops the current gate finds.**
This is a false-positive problem. Recall is not to be traded for it.

## The statistic

The current statistic is already a difference of a moving median, at a long lag. For odd `W` with
`h = floor(W/2)`, `movmedian(x, [W-1 0])` evaluated at index `i` is exactly the centred median-filter
output at `i - h`. Writing `mc` for that centred output:

    D(t) = mPost(t+g) - mPre(t-g) = mc(t + g + h) - mc(t - g - h)

a lag of `2(g + h) = TransitionMs + WindowMs`, i.e. 250 ms at the shipped defaults. The design knob
is therefore the **lag**, not the operator. A lag-`L` difference of a median filter reports
`A * min(1, L/T)` for a step of amplitude `A` and rise `T`, and `s * L` for drift of slope `s`:
shrinking `L` buys drift rejection and charges for it in slow-pop amplitude.

This design sets the lag to `TransitionMs`, so that the abruptness timescale and the transition
allowance are one parameter:

    mc(t)     = mPre(t + h)                 % an index shift of a value already computed
    Dslope(t) = mc(t + g) - mc(t - g)       % lag 2g ~ TransitionMs
    D(t)      = mPost(t + g) - mPre(t - g)  % unchanged; feeds the amplitude floor and the sign

Two properties make this cheap and safe:

- **Support.** `Dslope` reads `x` over `t-g-h .. t+g+h`, strictly narrower than `D`'s
  `t-g-2h .. t+g+2h`. The existing `validMask` margins (`W + g` from each segment edge) remain
  correct and conservative; no new edge handling, no change to `minSegLen`.
- **Window overlap is harmless.** At the defaults the two median windows overlap (100 ms lag,
  150 ms windows). A median filter is edge-preserving, so `mc` has already switched to the
  post-level once its centre passes the transition; the overlap costs no amplitude and reduces the
  noise of the difference, because overlapping supports share samples.

### The gate

Local baseline and variability of `Dslope`, from a wide span centred on `t`:

    b(t) = median(|Dslope|)      over +/- SlopeBaselineMs/2   (omitnan)
    s(t) = 1.4826 * MAD(Dslope)  over the same span           (omitnan)
    zslope(t) = (|Dslope(t)| - b(t)) / max(s(t), SlopeNoiseFloorMicrovolts)

`SlopeBaselineMs` starts development at 2000 (confirmed or moved by the sweep).
`SlopeNoiseFloorMicrovolts` is 0.1, the same guard and the same value the current local z uses
against a window with no spread; it is not swept. `Dslope` is NaN outside `valid`, so both
moving-window estimates run with `omitnan` and a span that is entirely invalid yields NaN, which
fails the comparison and produces no candidate - the behaviour the current code gets from the
`valid` mask.

A candidate sample requires all three of:

    zslope(t) > ZThreshold                                   (local, excess over baseline slope)
    |Dslope(t)| / channelScaleSlope > ChannelZThreshold       (channel-wide)
    |D(t)| >= MinAmplitudeMicrovolts                          (absolute floor, unchanged)

`channelScaleSlope` is `1.4826 * MAD(Dslope)` over the channel's locally-quiet samples, by the same
two-pass construction `gateCandidates` uses today: the first pass forms the provisional mask
`loc = valid & zslope > ZThreshold & |D| >= MinAmplitudeMicrovolts`, the scale is taken over
`valid & ~loc`, and the second pass adds the channel test. `loc` therefore survives as an internal
variable with a new definition - what disappears is its *export* as `info.unstableFraction`. As
today, a channel whose quiet samples have no spread gets `eps` so the channel gate says nothing
rather than dividing by zero.

Rationale for each piece:

- **Excess over a local baseline, not a ratio to local variability alone.** Steady drift produces an
  almost-constant non-zero `Dslope`, whose MAD is near zero. Normalising only by variability would
  make the local z explode on exactly the case this design exists to reject. Subtracting `b` handles
  steady drift; dividing by `s` handles noisy drift.
- **A wide centred span, not the adjacent pre/post windows.** `Dslope` evaluated at `u = t - g`
  equals `mc(t) - mc(t - 2g)`, which already contains half the step: windows flanking `t` are
  contaminated by the event being scored, inflating `b` and `s` and desensitising the gate where it
  matters. Over a span of ~2 s the contaminated samples are a ~15-20 % minority, and both the median
  and the MAD have 50 % breakdown, so the contamination is ignored without any offset arithmetic.
- **The absolute floor stays on `D`, not on `Dslope`.** `D`'s longer lag measures the settled level
  change, which is what "50 microvolts" is meant to mean. Keeping it upstream also stops the new
  gate promoting small abrupt wiggles into repairable intervals.

### What does not change

Candidate runs keep their shape. `Dslope` sits at the full amplitude `A` for every `t` whose two
read-points straddle the transition, so a run spans about `2g - T` samples (~100 ms at defaults).
`signedRuns` (with the sign taken from `Dslope`), `RefractoryMs` gap-closing, `mid0` as the run
midpoint, `refineCandidate`'s bracketing, the level re-measurement, `plateauDrift`, `holdDrift`,
`varianceRatio`, event grouping and the common-mode flag are all untouched. The `steps` field names
and the rejection reasons in `info.nRejected` are unchanged; `steps.z` and `steps.channelZ` keep
their names and acquire their new meaning.

### Known limit introduced

A step whose rise `T` exceeds `2g` never reaches full amplitude; it reports `A * 2g/T`. At
`TransitionMs = 100` that attenuates the slow tail of this dataset's pops (rises run out to about
100 ms). Because `TransitionMs` now carries two jobs, the survey sweeps it (100-200 ms) rather than
assuming the current default survives.

## Cost

Per channel, four moving-window passes before and after. The two `movmad(x, ...)` passes that fed the
old local z retire; one `movmedian(|Dslope|)` and one `movmad(Dslope)` over `SlopeBaselineMs` take
their place. `summarize_profile.m` confirms the wall-clock claim on a real file.

## Retiring the instability criterion

`info.unstableFraction` was the fraction of samples passing the old **local** gate, and
`MaxUnstableFractionToRejectChannel = 0.25` rejected channels on it - the criterion that catches
telegraph noise whose plateaus are too short to bracket as steps. The statistic it is built on
disappears with the old local gate, and the fraction is not meaningful under the new one. The
criterion is removed rather than re-derived:

| File | Change |
| --- | --- |
| `RELAX_detect_steps.m` | drop the `MaxUnstableFractionToRejectChannel` parameter and `info.unstableFraction` (`loc` stays, redefined, to select quiet samples) |
| `RELAX_rank_pop_channels.m` | `scores` drops to two columns (rate, transition fraction); the `case 3` reason string goes; doc comment updated |
| `RELAX_Wrapper_slurm.m` | config default, pass-through argument, and the `PopNChannelsRejectedByUnstable` covariate (both the count and the field assignment) |
| `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` | config line |
| `validate_pop_repair.m` | stage config and pass-through |
| `tests/` | `test_RELAX_detect_steps.m` (telegraph case), `test_RELAX_rank_pop_channels.m`, `test_RELAX_protect_eye_channels.m` fixtures asserting the `locally unstable` reason |

**Hypothesis to be checked, not assumed:** telegraph-noise channels should now fall out by **step
rate** instead, because each telegraph transition is abrupt and clears the amplitude floor, and the
slope statistic suits them better than the level statistic did. The survey checks this on the
channels that previously rejected by instability.

**Migration.** Removing `PopNChannelsRejectedByUnstable` changes the schema of
`RELAXProcessingExtremeRejections`, and the wrapper reloads aggregate `.mat` checkpoints and appends
rows to them. The next production run therefore writes to a **fresh `RELAXProcessed/` directory**.
This is the right call independently: a detector change means every already-cleaned file was cleaned
by a different pipeline, so mixing old and new rows in one aggregate table would be misleading even
if the schema matched.

## Validation

The claim to test is specific: *the steps this gate stops finding are drift, and the pops it used to
find it still finds.* On the five surveyed recordings (`m159, m497, m208, m583, m155`, production
204-channel montage), reported at `docs/superpowers/reports/` depth.

The old algorithm is pinned as a renamed copy, `slope_gate_survey/legacy/RELAX_detect_steps_level.m`
(renamed so it cannot shadow the production function), and both detectors run in one process on one
load of each file - the pattern `validate_pop_repair`'s `legacyTools` arm already uses for the
superseded derivative detector. This gives the two-arm comparison with no branch in production code,
no risk of the arms seeing different data, and one load per file rather than two.

1. **Step-table diff.** Steps matched across arms by channel label and overlapping `[onset, offset]`,
   classified **retained / lost / gained**. Per step: amplitude, `riseMs`, `plateauDrift`,
   `holdDrift`, and a driftiness measure (excursion over the surrounding 2 s divided by the step
   amplitude).
2. **Blind adjudication.** No ground truth exists. Sample ~40 lost and ~40 gained steps, render 10 s
   windows around each, and score them pop / drift / unclear without knowing which arm produced them.
   *Acceptance:* the large majority of **lost** steps score as drift, and **gained** steps are neither
   numerous nor obviously spurious.
3. **Threshold sweep** over `ZThreshold` x `ChannelZThreshold` x `TransitionMs` x
   `SlopeBaselineMs`, reporting steps, channels rejected and interpolated fraction, so defaults are
   chosen at a knee rather than inherited from a statistic they were not tuned for.
4. **Channel rejection.** Which channels each arm deletes and why; specifically whether the channels
   that previously rejected by instability still reject by rate.
5. **Reference-pop path.** `RELAX_detect_reference_pops` inherits the new gate on the cross-channel
   median (`RefPop_MinAmplitudeMicrovolts = 10`); re-check the counts in the 2026-09-21 reference-pop
   survey.
6. **Stage-level rerun.** `sbatch run_validate_pop_repair.sh` with the chosen defaults, so the
   post-1 Hz-highpass residual at filled transitions, the leave-one-out spatial-fill check and the
   ICA agreement numbers are measured at the new detection sites.

All MATLAB runs go to compute nodes via `sbatch` / `srun`, never a login node.

## Implementation order

Tests first at each step (`tests/run_all_tests`, run under `sbatch`).

1. `RELAX_detect_steps.m`: `levelStatistic` returns `D` and `Dslope`; `gateCandidates` becomes the
   three-way gate; add `SlopeBaselineMs`; drop the instability parameter, field and `loc`;
   `info.channelScaleD` becomes `info.channelScaleSlope` (no reference outside the file today -
   confirm with a grep before renaming); rewrite the header's detector description
   and known-limits paragraph.
2. `RELAX_rank_pop_channels.m`.
3. Plumbing: `RELAX_Wrapper_slurm.m`, `RELAX_SET_PARAMETERS_AND_RUN_slurm.m`, `validate_pop_repair.m`.
4. `slope_gate_survey/`: pinned legacy detector, diff and sweep script, adjudication plot renderer,
   `sbatch` wrapper.
5. Defaults from the survey; addendum to `2026-09-14-pop-artifact-repair-design.md`; the stage-2
   paragraph in `CLAUDE.md`; a report under `docs/superpowers/reports/`.
6. Stage rerun, results into that report.

### New unit tests

- an ideal step is detected with the correct amplitude and bracket;
- **a 200 microvolt linear ramp over 2 s is rejected** (the motivating case; the current detector
  finds it);
- a 100 microvolt sigmoid over 50 ms is detected;
- **a genuine step riding on strong drift is still detected** - the failure a baseline-subtracting
  gate invites, and the test that keeps `SlopeBaselineMs` honest: too narrow and the step
  contaminates its own baseline, too wide and the gate stops being local;
- a random-walk-only channel yields no steps;
- a step just below `MinAmplitudeMicrovolts` is rejected;
- two steps 300 ms apart are both detected (mutual baseline contamination);
- a step whose rise exceeds `2g` is attenuated as documented;
- `info` no longer carries `unstableFraction`, and `info.params` echoes `SlopeBaselineMs`;
- degenerate / flat-channel handling is unchanged.

## Decisions taken, with the alternatives rejected

| Decision | Alternatives rejected |
| --- | --- |
| Slope replaces the per-sample gate | An acceptance test applied after bracketing (smaller blast radius, but leaves drift generating candidates) |
| `MinAmplitudeMicrovolts` stays upstream on `D` | Dropping it entirely; keeping it only in `refineCandidate` |
| Lag = `TransitionMs` | Lag 1 (demotes this dataset's 10-100 ms pops by `A/T`); a separate `SlopeSpanMs` parameter; an empirically swept lag with no tie to an existing parameter |
| Excess over local baseline, scaled by local variability | Ratio to local baseline alone (ignores how variable the local slope is); keeping `hypot(sdPre, sdPost)` (leaves all drift rejection to the channel gate) |
| Instability criterion removed | Keeping it on the old level-based gate; recomputing and re-tuning it on the new gate |
| Covariate removed, fresh output directory | Keeping the field as a constant 0; field-union logic in the wrapper's checkpoint path |
| Edit `RELAX_detect_steps.m` in place | A `GateStatistic` switch removed later; a `_slope` sibling file (the suffix convention is for dataset/execution variants, and CLAUDE.md warns siblings drift apart) |
