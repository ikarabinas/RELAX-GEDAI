# Boundary-aware preprocessing for the RELAX-GEDAI SLURM pipeline

Date: 2026-09-08
Status: approved design, not yet implemented
Scope: `RELAX_Wrapper_slurm.m` and the functions it calls

## Problem

The input files for this pipeline
(`/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Aug26/tms_artifact_removed/`)
are continuous recordings from which TMS pulse-train artifact periods have already been excised with
`eeg_eegrej`, leaving `boundary` events that mark genuine discontinuities in the signal. Most of the
RELAX pipeline treats the record as one unbroken time series, so filters, sliding windows, wavelet
transforms and mask-padding operations run straight across those discontinuities.

### Measured properties of the input data

Two representative files (`m121_dlpfc_day1`, `m155_dlpfc_day3`):

| Property | Value |
| --- | --- |
| Sample rate / channels | 1000 Hz / 257 |
| Duration | ~510-535 s |
| Boundary events | 47-60 per file, `duration` ~= 1869 samples |
| Segment lengths | almost all **7.93 s**, plus one long pre-TMS segment (36-135 s) |
| Other events | none - the files contain only `boundary` events |

The 7.93 s inter-train segment length is the dominant design constraint. A 1 Hz acausal Butterworth
has edge transients on the order of 1-2 s per end, so roughly a quarter of every segment is edge.

## Audit: which stages are already safe

Verified by reading the installed EEGLAB 2026.0.0 source. These need no change:

| Stage | Mechanism |
| --- | --- |
| `pop_resample` | splits at boundaries and resamples each segment (`pop_resample.m:135-137`) |
| `clean_data_with_zapline_plus_eeglab_wrapper` | passes `zaplineConfig.boundaries` into zapline |
| `pop_eegfiltnew` / `firfilt` | splits at discontinuities (`firfilt.m:64-71`), pads each segment with a repeated edge sample |
| `pop_epoch` | **drops** any epoch containing a boundary event (`pop_epoch.m:414-433`) |

`pop_epoch`'s behaviour is what protects `RELAX_epoching`,
`RELAX_excluding_channels_and_epoching`, `RELAX_excluding_extreme_values`, `RELAX_muscle` and
`RELAX_drift`: epochs that would straddle a gap never reach them. Step 6 below verifies this rather
than assuming it.

`RELAX_ICA_subtract` is also safe: it is ICA decomposition plus `pop_subcomp`, a purely spatial
projection with no time-domain filtering and no wavelets. ICA treats samples as exchangeable, so
discontinuities do not smear.

## Audit: which stages are broken

On the currently configured path (`FilterType='Butterworth'`, `NotchFilterType='ZaplinePlus'`,
`Do_MWF_*=0`, `Run_GEDAI=0`, `Perform_targeted_wICA=1`):

1. **`RELAX_filtbutter.m`** - a plain per-channel `filtfilt` over the whole record. Used for the
   1-200 Hz bandpass in the wrapper and for a 1-25 Hz bandpass inside `RELAX_blinks_IQR_method`.
2. **PREP `findNoisyChannels`** - see below; the most damaging of the set.
3. **`RELAX_blinks_IQR_method.m`** - inherits (1), plus `+/-400 ms` peak windows and `BlinkMaskFocus`
   mask spans that can cross boundaries.
4. **`RELAX_targeted_wICA.m`** - four boundary-sensitive time-domain operations: a level-5 `coif5`
   `swt`/`iswt` over the concatenated record, a 0.5-25 Hz `filtfilt` for blink localisation, a 15 Hz
   `filtfilt` for muscle components, and blink-mask padding loops plus a `movmean` that walk across
   boundaries. It also calls `RELAX_blinks_IQR_method` internally.
5. **The SER/ARR marking block**, `RELAX_Wrapper_slurm.m:366-382` - runs `RELAX_muscle`,
   `RELAX_horizontaleye`, `RELAX_drift` and two `RELAX_pad_brief_mask_periods` calls, all
   boundary-blind. See "dead code" below.

Latent (same class of bug, currently disabled): MWF delay-stacking in
`RELAX_perform_MWF_cleaning.m`, the Butterworth notch branch, `RELAX_wICA_on_ICLabel_artifacts.m`,
and GEDAI.

### The SER/ARR marking block is dead code

`RELAX_Wrapper_slurm.m:366-382` exists solely to build `NoiseMaskFullLengthR1`. The only consumer of
that field is `RELAX_metrics_final_SER_and_ARR`, whose call is **commented out** at
`RELAX_Wrapper_slurm.m:756`. The block additionally has a latent bug: `Marking_all_artifacts_for_SER_ARR`
is grown from nothing by logical indexing, so the NaN structure of the mask never survives.

Guarding this block off removes four of the broken functions from the active path at no scientific
cost, and saves compute.

### PREP is quantifiably unsafe on this data

`findNoisyChannels.m:78` states its own precondition: *"No segments of the EEG data have been
removed."* Specifically:

- It slices correlation windows as a flat tiling from sample 1 (`findNoisyChannels.m:197`), so
  windows straddle boundaries.
- `badTimeThreshold` defaults to **0.01** - a channel is rejected if more than 1% of its 1 s windows
  are bad.
- With 47-60 boundaries across ~535 one-second windows, **~9-11% of windows straddle a
  discontinuity** - nearly an order of magnitude over the threshold.
- It also runs a 45-50 Hz `filtfilt_fast` over the whole record (`findNoisyChannels.m:165-168`)
  before both the correlation and RANSAC stages, so every boundary step produces a transient.

The consequence is that channel rejection on this data is close to arbitrary. This is a fix, not an
accepted residual.

## Decisions taken

1. **Reflect-pad per segment.** Each inter-boundary segment is processed independently, with its own
   data mirrored onto each end before filtering and trimmed afterwards. Segments never exchange
   information, so a bad segment cannot contaminate its neighbours.
2. **Shared helpers plus in-place edits behind a flag.** New helper files, and modifications to the
   existing shared functions gated on `RELAX_cfg.RespectBoundaries` (default `0`). Sibling wrappers
   (`RELAX_Wrapper.m`, `RELAX_Wrapper_bart.m`) keep byte-identical behaviour until they opt in. This
   avoids the `_iTBS`/`_bart` duplication problem described in `CLAUDE.md`.
3. **Synthetic fixture with invariant assertions, then a real-file spot check** per step.
4. **Fix the active path; defer, document and guard the rest.**

## Design

### Shared primitives (new files)

**`RELAX_get_segments(EEG)`** -> `[nSeg x 2]` matrix of `[startSample endSample]`.

Splits at `boundary` events. Uses the integer-safe form

```matlab
starts = [1, floor(boundaries) + 1];
ends   = [floor(boundaries), EEG.pnts];
```

which handles both EEGLAB's `x.5` latency convention and tools that write integer boundary latencies,
without duplicating a sample across two segments. (Note this differs from the `ceil`/`floor` form in
`~/TMS_response_variability/iTBS_pipeline.m`, which duplicates a sample for integer latencies.)
Segments with `start > end` are dropped. Data with no event field, or with no boundary events,
yields a single whole-record segment - so every caller degrades to exactly the current behaviour on
ordinary continuous data.

**`RELAX_reflect_pad(x, p)`** -> odd reflection about the edge sample, matching `filtfilt`'s own
internal convention:

```matlab
left  = 2*x(:,1)   - x(:, p+1:-1:2);
right = 2*x(:,end) - x(:, end-1:-1:end-p);
```

Odd rather than even reflection matters for the 1 Hz highpass: even reflection creates a cusp at the
edge that the highpass then rings on, whereas odd reflection preserves the local trend. `p` is
clamped to `segLen-1`.

**`RELAX_apply_by_segments(data, segs, fn, p)`** - for each segment: reflect-pad by `p`, apply the
function handle `fn`, trim back to the original extent, write into the output.

Pad width is derived from the governing (lowest) frequency of the operation:

```matlab
p = min(round(3 * srate / f_low), segLen - 1);
```

For the 1 Hz highpass at 1000 Hz that is 3000 samples against a 7930-sample segment. The padding is
trimmed, so the only cost is roughly 3x working memory per segment.

**`RELAX_check_boundary_support(EEG, RELAX_cfg)`** - startup guard, see "Deferred stages" below.

**Design rule:** every caller re-derives its segment list from the *current* `EEG.event` and never
caches one across stages. `eeg_eegrej` at `RELAX_Wrapper_slurm.m:657` adds boundaries mid-run, so the
re-reference and wICA stages see more discontinuities than the filter stage did.

### Configuration flag

`RELAX_cfg.RespectBoundaries`, default `0`, set to `1` only in
`RELAX_SET_PARAMETERS_AND_RUN_slurm.m`. Every modified function checks the flag and takes the
existing code path when it is `0` or absent.

### Implementation steps

Each step is implemented and tested before the next begins.

#### Step 0 - Infrastructure and fixture

New files: `RELAX_get_segments.m`, `RELAX_reflect_pad.m`, `RELAX_apply_by_segments.m`,
`RELAX_check_boundary_support.m`, `tests/make_boundary_fixture.m`, `tests/run_all_tests.m`.

No behaviour change to the pipeline. Tests: segment splitting for integer and half-integer boundary
latencies; data with no events yields one segment; reflect-pad output shape and values; pad width
clamping on short segments.

#### Step 1 - Guard the dead SER/ARR block

Gate `RELAX_Wrapper_slurm.m:366-382` behind a new `RELAX_cfg.ComputeSER_ARR` flag defaulting to `0`,
with a comment recording that `RELAX_muscle`, `RELAX_horizontaleye`, `RELAX_drift` and
`RELAX_pad_brief_mask_periods` are boundary-blind and must be fixed before it is re-enabled.

Done early because it simplifies the path that every later step runs through.

Test: run the wrapper on the fixture and assert `NoiseMaskFullLengthR1` is absent and all other
outputs are unchanged.

#### Step 2 - `RELAX_filtbutter`

Add an optional 7th argument `respect_boundaries` (default `false`). When set, derive segments from
`EEG.event` and route the filter through `RELAX_apply_by_segments`. All existing call signatures are
preserved.

Optional, same lines: replace the per-channel loop with a single matrix `filtfilt(z1,p1,data')'`.
This is a meaningful speed-up at 257 channels x 535 s. Keep it only if the backward-compatibility
assertion still passes bitwise; otherwise drop it and leave the loop alone.

Callers to opt in: the bandpass/highpass/lowpass calls in `RELAX_Wrapper_slurm.m`.

#### Step 3 - PREP `findNoisyChannels`

New file `RELAX_findNoisyChannels_by_segments.m`. `findNoisyChannels` returns *channel* indices, so a
rearranged copy of the data can be fed to it and the result mapped straight back.

1. Take the already segment-wise-highpassed data.
2. Split into segments; trim each to a whole multiple of **5 s** - the LCM of
   `correlationWindowSeconds` (1) and `ransacWindowSeconds` (5) - and drop segments shorter than 5 s.
3. Demean each trimmed segment per channel, so the joins carry no DC step for PREP's 45-50 Hz FIR to
   ring on.
4. Concatenate into a temporary EEG struct and call `findNoisyChannels` **once**.
5. Apply the returned channel list to the real EEG with `pop_select`.

Every PREP window then lies wholly inside one original segment. Data budget: ~47 x 5 s plus the long
segment gives ~270-370 s of the original 535 s, ample for channel detection. The ~37% discarded from
each short segment is trim waste, not data removed from the recording.

Called from `RELAX_Wrapper_slurm.m:324` when the flag is set.

#### Step 4 - `RELAX_blinks_IQR_method`

Pass `RELAX_cfg.RespectBoundaries` through to its internal `RELAX_filtbutter` call. Clamp the
`+/-400 ms` peak-search windows and the `BlinkMaskFocus` mask spans so they cannot extend past the
end of the segment containing the blink.

This function stays in the active path because `RELAX_targeted_wICA.m:297` consumes
`EEG.RELAX.eyeblinkmask`.

#### Step 5 - `RELAX_targeted_wICA`

- Both `filtfilt` calls (0.5-25 Hz blink localisation, 15 Hz muscle) route through
  `RELAX_apply_by_segments`.
- `swt`/`iswt` applied per segment. `swt` at level 5 requires the length be divisible by 32; the
  existing code zero-pads globally. Per segment, reflect-pad by at least 512 samples (coif5's level-5
  effective support is about 30 * 2^4 = 480 samples, ~0.5 s), round up to a multiple of 32, transform,
  then trim. 7930-sample segments support this comfortably.
- Clamp the blink-padding loops and the `movmean` to segment extents.

#### Step 6 - Verify the epoching stages

Assertions only, no code change expected. On the fixture, confirm that `pop_epoch` drops every
boundary-crossing epoch; that `originallatency` bookkeeping stays aligned after
`pop_select('notrial', ...)`; and that the resulting `NaNsForNonEvents` mask is NaN wherever an epoch
was dropped. If any assertion fails, this becomes an implementation step.

#### Step 7 - Guards and documentation

`RELAX_check_boundary_support(EEG, RELAX_cfg)` **errors** at startup when boundary events are present
and any of the following is enabled:

- `Do_MWF_Once` / `Do_MWF_Twice` / `Do_MWF_Thrice` - delay-stacking in `RELAX_perform_MWF_cleaning`
  builds one lagged covariance across the whole record, so lags span boundaries
- `Run_GEDAI`
- `NotchFilterType = 'Butterworth'`
- `Perform_wICA_on_ICLabel` - wavelet-based, same problem as targeted wICA
- `ComputeSER_ARR` (the block guarded in Step 1)

`Perform_ICA_subtract` is safe and is **not** guarded.

`CLAUDE.md` gains a "Boundary/discontinuity handling" section recording which stages are
boundary-aware and by what mechanism (including the four EEGLAB functions that handle it natively),
and which are deferred.

## Testing

`tests/make_boundary_fixture.m` builds a synthetic `.set`: 8 channels, 1000 Hz, ~60 s, six boundaries
at realistic 7.93 s spacing plus one 35 s segment. Content is a 10 Hz sinusoid plus pink noise plus
known per-segment DC offsets and slow drifts, so cross-boundary smearing is directly measurable.

`tests/run_all_tests.m` runs the per-step `test_*.m` scripts and prints PASS/FAIL per assertion.

Three invariants recur at every step:

- **Backward compatibility** - with no boundary events present, output is bitwise-identical to the
  pre-change function.
- **Isolation** - scale segment *k+1* by 100x and assert segment *k*'s output is unchanged to machine
  precision. This is the assertion that actually proves nothing leaks across a gap.
- **Bookkeeping** - boundary count, boundary latencies, `pnts`, and mask/NaN alignment all preserved.

Step-specific checks:

| Step | Check |
| --- | --- |
| 0 | helper unit tests (above) |
| 1 | `NoiseMaskFullLengthR1` absent; all other wrapper outputs unchanged |
| 2 | a 10 Hz sine inside a segment survives the 1-200 Hz bandpass within 1% of standalone-segment filtering |
| 3 | rejected-channel list before vs. after, on a real cropped file |
| 4 | no blink mask span crosses a boundary; blink count per segment stable under the isolation test |
| 5 | `iswt(swt(x))` reconstructs `x` per segment within tolerance; `artifact_comp` length preserved |
| 6 | epoch-drop and `originallatency` assertions (above) |
| 7 | guard errors when each deferred flag is set with boundaries present, and does not when they are not |

After each step passes on the fixture, one spot check on a cropped real file before moving on.

## Files touched

New:

- `RELAX_get_segments.m`
- `RELAX_reflect_pad.m`
- `RELAX_apply_by_segments.m`
- `RELAX_check_boundary_support.m`
- `RELAX_findNoisyChannels_by_segments.m`
- `tests/make_boundary_fixture.m`, `tests/run_all_tests.m`, `tests/test_*.m`

Modified:

- `RELAX_filtbutter.m`
- `RELAX_blinks_IQR_method.m`
- `RELAX_targeted_wICA.m`
- `RELAX_Wrapper_slurm.m`
- `RELAX_SET_PARAMETERS_AND_RUN_slurm.m`
- `CLAUDE.md`

Deferred, guarded, unmodified:

- `RELAX_muscle.m`, `RELAX_drift.m`, `RELAX_horizontaleye.m`, `RELAX_pad_brief_mask_periods.m`
- `RELAX_perform_MWF_cleaning.m`
- `RELAX_wICA_on_ICLabel_artifacts.m`
- GEDAI integration
