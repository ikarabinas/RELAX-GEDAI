# Reference-electrode pop survey — 2026-09-21

Go/no-go evidence for the reference-pop repair stage (plan:
`docs/superpowers/plans/2026-09-21-reference-pop-repair.md`, Task 4). A reference pop is a
level change on the reference contact (E257 / Cz on these files), which by `x_i = s_i − r`
appears with the same sign, amplitude and waveform on every channel. The question the survey
answers: do such events exist in this dataset, how often, and do they look like that model
says — or is the common-mode step class something else (ground pops, partial movements)?

## Job

| | |
|---|---|
| Job | 19834319, partition `scu-cpu`, node c15, wall time 00:10:04 (64 GB requested) |
| Script | `survey_reference_pops.m` via `run_survey_reference_pops.sh` (commit 18b3d69) |
| Log | `logs/refpop_survey_19834319.log` |
| Data | `~/TMS_response_variability/store/mdd_dlpfc/preprocessed_Sep26/rest_treat_joined/` — `m159`, `m497`, `m208`, `m583`, `m155` `_dlpfc_day1_stitched.set`, production 204-channel montage (203 usable after E257 is dropped), 1000 Hz, 17.6–21.3 min each |
| Outputs | `reference_pop_survey.mat` (uncommitted, 38 MB), figures in `docs/superpowers/reports/reference_pop_survey/` |
| Extra | duplicate / overlap check on the saved candidates, job 19836367 (`docs/superpowers/reports/reference_pop_survey/check_dups.m`, output in `logs/check_dups_19836367.log`) |

Per-file cost of the new machinery: virtual channel + detector + uniformity ≈ 10–30 s at
either floor; the existing per-channel detector ≈ 50–80 s. The plan's 1–3 h estimate for the
survey was off by an order of magnitude.

## Method

`RELAX_reference_channel` builds the virtual reference as the cross-channel median of the
per-segment-demeaned usable channels. `RELAX_detect_reference_pops` runs `RELAX_detect_steps`
on it (shipped detector config, `MaxPlateauDriftFraction = Inf`) with the uniformity gate
**off**, so every candidate is kept together with its `RELAX_step_uniformity` diagnostics
(`signAgreement`, `amplitudeSpread` = 1.4826·MAD/|median|, `fractionNear` = fraction of
channels within 25 % of the median amplitude); the gate (`signAgreement ≥ 0.9`,
`amplitudeSpread ≤ 0.5`) is then applied after the fact. Two floors were run (30 µV, the
shipped `PopDetect_MinAmplitudeMicrovolts`, and 10 µV). The existing per-channel detector ran
at the shipped config (`CommonModeChannelFraction 0.15`) and its common-mode events were
matched to the reference candidates with `RELAX_match_intervals` (50 ms slack); the unmatched
common-mode events got the same uniformity diagnostics.

## Prevalence (floor 30 µV, gate 0.9 / 0.5)

| file | min | candidates | gated | CM events | CM = ref | median agree (gated) | median spread (gated) |
|---|---|---|---|---|---|---|---|
| m159 | 17.6 | 30 | 1 | 67 | 25 | 0.99 | 0.21 |
| m497 | 21.1 | 14 | 8 | 223 | 6 | 0.99 | 0.16 |
| m208 | 20.4 | 31 | 24 | 59 | 27 | 0.94 | 0.23 |
| m583 | 21.3 | 55 | 2 | 188 | 51 | 1.00 | 0.14 |
| m155 | 17.7 | 0 | 0 | 0 | 0 | – | – |
| **total** | 98.1 | 130 | **35** | 537 | 109 | | |

35 gated reference pops in four of five files (0.36/min overall, but concentrated: m208 alone
has 24 in 20 min, in 11 episodes; m497's 8 are 2 episodes). |amplitude| of the gated set
ranges 31–301 µV, median rise 70 ms (2–172 ms).

## Assumption check (gated candidates, floor 30)

| metric | bar | all 35 | m159 (1) | m497 (8) | m208 (24) | m583 (2) |
|---|---|---|---|---|---|---|
| median `signAgreement` | ≥ 0.95 | **0.98** ✅ | 0.99 | 0.99 | 0.945 | 1.00 |
| median `amplitudeSpread` | ≤ 0.25 | **0.20** ✅ | 0.21 | 0.16 | 0.225 | 0.14 |
| median `fractionNear` | ≥ 0.75 | **0.74** ❌ (17/35 ≥ 0.75) | 0.74 | 0.88 | 0.70 | 0.84 |

Visually (`summary_uniformity.png`): the candidates fall in two well-separated populations.
Top-left (agree > 0.9, spread < 0.5) is the reference-pop class; the other cluster sits at
agree ≈ 0.73, spread 0.8–1.3 and coincides with the unmatched common-mode events. The 0.9/0.5
gate sits in the gap between them; only a handful of points lie within 0.05 of either line.

What the two classes are, from the 24 candidate figures:

- **Gated = reference pops.** `m583_cand04` and `m583_cand24` are textbook: a 26–86 ms step
  of −66 µV on all 203 channels, black median and grey channels the same shape, the amplitude
  bar chart flat within ±15 µV of the red median line. `m497_cand05/06` likewise (spread 0.08,
  near 0.93).
- **Rejected = blinks seen through the Cz reference.** In every rejected figure (`m159_cand08/14/17/21`,
  `m583_cand06/46/47`, `m497_cand13/14`) a few frontal channels rise by +200…+800 µV while
  the other ~150 channels move −60…−90 µV *together* — the blink potential picked up at Cz and
  subtracted from every channel. The amplitude bars are a smooth gradient from −100 to +400,
  never a plateau, and the events come in opposite-sign pairs 250–700 ms apart (a pulse, not a
  step). Sign agreement ≈ 0.73 is the fraction of channels on the far side of the blink
  dipole. None of these should be subtracted from the montage — average re-referencing and the
  blink stage handle them later — and the gate rejects all of them (0 unmatched common-mode
  events pass the gate on m159/m497/m208, 2 on m583).

Why `fractionNear` misses its bar: m208. That recording is grossly unstable — individual
channels railing by thousands of µV (`m208_cand02/04`), common-mode ramps of ~1 mV/s
(`m208_cand26/27`) — so 150 ms window medians on the individual channels scatter around the
common step even when the step itself is uniform (bars in `m208_cand21`: −250…−350 on 190
channels, spread 0.25, near 0.71). On the two quieter files the gated pops have near 0.84–0.94.
The failure is a property of m208's noise floor, not evidence of a ground-pop / partial class:
the amplitude distributions of the m208 gated candidates are single plateaus with outliers,
not gradients.

Two further things the figures show that the metrics do not:

1. **Not every gated candidate is a step.** `m497_cand03/04/05/06` (663.5–665 s) are the rising
   and falling flanks of one ~2 s common-mode slow wave (±100 µV). They pass the gate because
   the movement is uniform across channels, and the detector accepts them because the shipped
   `MaxPlateauDriftFraction` is `Inf` (the plateau-flatness test is diagnostic only). The same
   holds for `m208_cand26/27`, steps riding on a 1 mV/s ramp. Subtracting only the *step*
   content of such a candidate (Task 5's design: transition waveform, then the constant
   amplitude to the segment end) leaves the rest of the wave and, where the falling flank is
   not also detected, introduces an offset.
2. **Gated candidates overlap and duplicate.** Job 19836367: m497 and m208 each have 3 gated
   pairs whose `[onset, offset]` intervals overlap; m583 has 2 exact duplicates (identical
   onset/offset, different z — both ungated here, but nothing prevents a gated one). Task 5's
   stacking rule (each step adds its own transition and amplitude) double-counts both cases.

## Common-mode overlap (the gain the stage offers)

109 of 537 common-mode events (20 %) in the per-channel detector's output are reference pops;
the other 80 % are the blink class above. Per file: m159 25/67, m497 6/223, m208 27/59,
m583 51/188. The reference detector saw 126 of its 130 floor-30 candidates as common-mode
events too (the 4 misses are on m159), so the per-channel stage already sees every reference
pop — today it excises those shorter than `PopRepair_MaxExciseMs` and leaves the longer ones;
the new stage would subtract them instead and keep the data.

## Floor 10 µV

Floor 10 gates 48 candidates vs 35 at floor 30: +1 on m159 (21 µV), +2 on m497 (25, 30 µV),
+10 on m583 (15–24 µV), 0 on m208. So sub-30 µV reference pops exist and are the majority on
m583. They are invisible to the per-channel stage (its floor is 30 µV) and, unrepaired, a
20 µV common step still becomes a multi-second bipolar excursion after the 1 Hz zero-phase
highpass on every channel.

## Virtual-channel noise

`robustSd` of the median vs the median channel's own `robustSd`: 78/115, 24/116, 210/259,
101/172, 151/167 µV — a ratio of 0.2–0.9, not the ≈ 1/√203 ≈ 0.07 the hop-noise concern
anticipated. The within-segment variability of these raw Cz-referenced files is dominated by
drift that is common to all channels (i.e. reference drift), which the median inherits by
construction. Independent per-channel noise is averaged away as intended; the detector's
channel-wide z (MAD over quiet samples) is what keeps that drift from flooding it, and the
candidate counts (0–55 per file) show it does.

## Decision: DISCUSS (prevalence says GO)

- Prevalence: 35 gated in 4/5 files clears "≥ 5 total and ≥ 2 files" by a wide margin → GO.
- Assumptions: 2 of 3 medians pass; `fractionNear` 0.74 vs 0.75, attributable to m208's
  noise. The figures agree with the model for the gated class and show the rejected class is
  blinks, not a ground-pop / variable-amplitude class. Per the plan a failed assumption means
  stop and discuss rather than proceed silently, and the two figure findings above bear
  directly on Task 5's subtraction rule, so this report stops here.

Recommendation: proceed to Tasks 5–7 with two changes to the design, which are the actual
outcome of the survey:

1. **Merge overlapping candidates and drop exact duplicates before subtracting** (in
   `RELAX_detect_reference_pops` or at the top of `RELAX_subtract_reference_pops`): union the
   `[onset, offset]` of overlapping candidates into one step whose amplitude is the post-level
   of the last minus the pre-level of the first. Without this the m497/m208 clusters are
   corrected twice.
2. **Gate on plateau flatness for the reference detector** (a finite
   `MaxPlateauDriftFraction`, e.g. 0.5, passed through `Detect`), or surface
   `steps.plateauDrift` in `refInfo.candidates` and gate on it, so slow-wave flanks are not
   treated as steps. Task 7's residual check will show whether this is enough.

Also worth a look, outside this plan: `RELAX_detect_steps` emitting the same `[onset, offset]`
twice on one channel with different z (m583 candidates 32/33 and 34/35).

## Parameters to carry into Task 6

| parameter | value | basis |
|---|---|---|
| `RefPop_MinChannels` | 8 | 203 channels used on every file; the floor only matters for degenerate montages |
| `RefPop_MinAmplitudeMicrovolts` | 30 (= `PopDetect_MinAmplitudeMicrovolts`) for the first pass | keeps Task 7's before/after comparison against the per-channel stage like-for-like; floor 10 finds 13 more real pops (mostly m583, 15–24 µV) and is the obvious follow-up once the subtraction is validated |
| `RefPop_MinSignAgreement` | 0.9 | sits in the gap between the two populations; `m159_cand05` (0.88 / 0.28, a genuine −62 µV step on ~170 channels) is the one visible casualty, so 0.85 is defensible too |
| `RefPop_MaxAmplitudeSpread` | 0.5 | same; no rejected-class point lies below it with agree ≥ 0.9 except 2 on m583 |
| `MaxPlateauDriftFraction` for the reference detector | finite (0.5 proposed) | recommendation 2 above; not in the original plan |

## After subtraction (Task 7, 2026-09-22)

**Run:** `run_survey_reference_pops.sh` → job **19844528** (19 min,
`logs/refpop_survey_19844528.log`), output `reference_pop_survey_task7.mat` (untracked),
figures in `reference_pop_survey/task7/`. The subtraction arm now does exactly what
`RELAX_Wrapper_slurm.m` does: `RELAX_detect_reference_pops` with the uniformity gate on
(0.9 / 0.5), floor 10 µV, `MaxPlateauDriftFraction 0.3`, `HoldMs 500`,
`MaxHoldDriftFraction Inf` (the level-hold gate was surveyed and not switched on —
`2026-09-21-hold-drift-survey.md`), overlapping candidates merged, then
`RELAX_subtract_reference_pops` with the detector's own median trace and `WindowMs 150`;
afterwards both detectors run again on the corrected data. Per-site diagnostics:
`reference_pop_survey/task7/residuals.m` (job 19845189) and `still_cm.m` (job 19845212).

| file | subtracted (merged) | residual gated | common-mode events before → after | at subtracted sites still common-mode | per-channel steps before → after |
|---|---|---|---|---|---|
| m159 | 2 (0) | 0 ✅ | 67 → 67 ✅ | 0 ✅ | 8790 → 8786 ✅ |
| m497 | 7 (2) | **1** ❌ | 186 → 183 ✅ | **1** ❌ | 22 334 → 22 039 ✅ |
| m208 | 18 (3) | **1** ❌ | 51 → 40 ✅ | **6** ❌ | 8252 → 6355 ✅ (−23 %) |
| m583 | 11 (1) | 0 ✅ | 177 → 173 ✅ | 0 ✅ | 21 755 → 21 275 ✅ |
| m155 | 0 | 0 ✅ | 0 → 0 ✅ | 0 ✅ | 1199 → 1199 ✅ |

Nothing got worse anywhere: no file gained common-mode events or per-channel steps, and
m208 — the file with the most reference pops — lost a quarter of its per-channel steps and
11 of 51 common-mode events. Two of the four criteria fail on two files, for two reasons
that the per-site diagnostic makes specific:

**The residual gated candidates are the slow-wave flanks, re-bracketed.** m497's is at
664.74 s (−28 µV), between the subtracted 664.55 and 664.91 s treads of the flank the
first survey flagged; m208's is at 260.96 s (−85 µV), where the plateau gate had correctly
rejected candidate 17 (plateau 0.43) and subtracting the 259.74 / 260.63 s treads left a
remainder that brackets with plateau 0.23. Subtracting "steps" off a wave leaves more wave.
Neither is a pop the stage missed; both are the known residual of gating on 150 ms windows.

**Sites that stay common-mode are the ones whose per-channel amplitudes were never
uniform.** From `still_cm.m` (per-channel common-mode event overlapping each subtracted
step, before → after):

| m208 site | amplitude | spread | before | after |
|---|---|---|---|---|
| 222.40 / 222.73 s | −231 / +215 µV | 0.20 / 0.46 | 1 ev, 188 ch, median \|amp\| 248 | 1 ev, 93 ch, median 245 |
| 227.51 s | −233 | 0.40 | 176 ch, 269 | 67 ch, 162 |
| 281.40 s | −95 | 0.38 | 99 ch, 77 | 50 ch, 122 |
| 288.52 s | −301 | 0.25 | 180 ch, 310 | 48 ch, 244 |
| 720.73 s | −193 | 0.32 | 68 ch, 153 | 73 ch, 155 (slow-wave flank) |
| the other 12 sites | | 0.02–0.44 | 46–195 ch | **none** |
| m497 868.41 s | −44 | 0.48 | 155 ch, 75 (one channel at 94 mV) | 217 ch, 71 |

Every site with spread ≤ 0.19 is gone after subtraction, and so are several with spread
0.3–0.44 (719.92, 720.42, 855.31 s). What remains is where a large fraction of channels
carried a step of a very different size from the median — the pair at 222 s keeps a
245 µV median step on 93 channels after a 231 µV subtraction, i.e. those channels had
roughly twice the reference amplitude. That is not a reference pop in the sense of the
design (identical on every channel); it passes the 0.5 spread gate because the gate is an
IQR-type measure and the distribution is wide rather than shifted. On those sites the
per-channel stage still sees a common-mode event and excises it, exactly as it did before
the subtraction — with fewer channels involved.

**Verdict on the plan's criteria:** not met on m497 and m208
(`nResidualGated`, `nSubtractedStillCommonMode`). Per Task 7 step 3,
`Perform_reference_pop_repair` stays **0** and the decision is the user's. The evidence
for switching on anyway: the subtraction removed 11 of 18 m208 sites and all 7 uniform
m497 sites completely, reduced common-mode events and per-channel steps on every file that
had reference pops, and left the remaining sites to the per-channel excision that already
handled them. The evidence for holding: the residuals are a real class (wide-spread
"reference" events and slow-wave flanks) that a one-pass subtraction turns from one clean
excision into a subtraction plus a smaller excision, and nobody has yet looked at what the
corrected data look like at those sites at the ICA stage. Two knobs would shrink the
residual class before a second look: `RefPop_MaxAmplitudeSpread` 0.5 → ~0.25 (loses
719.92 / 720.42 / 855.31 s on m208, which subtracted cleanly, so it is a trade), and a
merged-event test for slow waves (see the hold-drift survey's closing paragraph).
