# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

RELAX-GEDAI is a fork of the RELAX EEG artifact-cleaning pipeline (MATLAB/EEGLAB) that adds GEDAI
(Generalized Eigenvalue De-Artifacting Instrument) lead-field filtering as an additional/alternative
denoising stage. It is a script-driven MATLAB toolbox, not a packaged library: there is no build step,
package manager, or test suite — correctness is verified by running the pipeline on real EEG data and
inspecting the cleaning-quality metrics it outputs.

This is a lab research codebase used by multiple people (originally Neil Bailey's RELAX, extended by
Isabella Karabinas for GEDAI, further adapted here for TMS-EEG / SLURM use). Expect dataset-specific,
hardcoded paths inside scripts (`RELAX_cfg.myPath`, `addpath(...)` calls to EEGLAB/plugins) — these are
meant to be edited per-run/per-user, not treated as bugs.

`~/TMS_response_variability` is a separate repo holding most of the code for a scientific experiment
that consumes this pipeline's output — RELAX-GEDAI cleans the raw TMS-EEG recordings, and that
experiment's code is the downstream consumer/analysis layer. When working across both repos, treat
this repo as the preprocessing stage and `~/TMS_response_variability` as where the resulting cleaned
data gets used for the actual study.

## Running the pipeline

There is no CLI entry point; you run these as MATLAB scripts (typically via `run(...)` or by hitting
"Run" in the MATLAB editor):

- `RELAX_SET_PARAMETERS_AND_RUN.m` — the canonical entry point. Sets `addpath`s for EEGLAB/PREP/
  MWF/FieldTrip/FastICA, configures `RELAX_cfg` (paths + all pipeline parameters), builds the file
  list, then calls `RELAX_Wrapper(RELAX_cfg)`.
- `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` / `RELAX_Wrapper_slurm.m` — SLURM job-array variant. Instead
  of processing every file in one MATLAB process, it groups files by participant ID and uses
  `SLURM_ARRAY_TASK_ID` (via `getenv`) to select a single participant per array task, and gives each
  array task its own `parpool` `JobStorageLocation` under `$TMPDIR` to avoid clobbering the shared
  default parpool job-storage dir when many tasks start concurrently. Submit with
  `sbatch run_relax-gedai_slurm.sh` (adjust `--array=1-N%M` to the participant count).
- `run_relax-gedai.sh` — non-array SLURM submission wrapping the plain (non-slurm) entry point.
- `RELAX_SET_PARAMETERS_AND_RUN_bart.m` / `RELAX_Wrapper_bart.m` — a parallel variant configured for
  a different dataset/experiment ("bart"). Same architecture as the main pipeline; parameter/path
  differences only.
- `AARATEP_RELAX_tmseeg_pipeline.m` — combines steps from the AARATEP TMS-EEG pipeline (in the
  `AARATEPPipeline` git submodule — run `git submodule update --init` if that directory is empty)
  with RELAX steps, for preprocessing continuous TMS-EEG recordings.
- `convert_files_to_eeglab_format/` — standalone conversion utilities (`save_mff_to_set*.m`,
  `save_fif_to_set.m`, `run_set_conversion.py`) to turn raw EGI `.mff` / `.fif` recordings into the
  EEGLAB `.set` format the main pipeline expects. Run the Python driver with
  `python run_set_conversion.py` after loading a MATLAB module (it shells out to MATLAB).
- `find_missing_files.m` — diff input `.set` files against `RELAXProcessed/Cleaned_Data` output to
  find files a run failed to process.
- `recover_metrics.m` — rebuilds aggregate cleaning-metric `.mat`/`.csv` files from individual
  processed files, for when a run was interrupted before it could write the aggregate files itself.

Because everything is a script with a hardcoded `RELAX_cfg.myPath`/`addpath` section, when adapting a
script for a new user/dataset, update the paths at the top rather than adding config indirection.

### EEGLAB plugin mode

The toolbox also installs as an EEGLAB plugin (`eegplugin_RELAX.m` registers menu items;
`pop_RELAX.m` is the GUI-driven entry point equivalent to `RELAX_SET_PARAMETERS_AND_RUN.m`;
`pop_RELAX_epoch_the_clean_data.m` / `RELAX_epoch_the_clean_data_Wrapper.m` handle the post-cleaning
epoching step from the GUI; `pop_RELAX_help.m` / `pop_RELAX_citation.m` are informational dialogs).

### External dependencies (not vendored, expected on the MATLAB path)

- EEGLAB (with the PREP plugin, ICLabel, PICARD, and FastICA plugins)
- MATLAB Signal Processing Toolbox and Statistics and Machine Learning Toolbox
- MWF artifact removal toolbox (github.com/exporl/mwf-artifact-removal)
- FieldTrip
- GEDAI (github.com/neurotuning/GEDAI-master or ikarabinas/GEDAI-master) when `RELAX_cfg.Run_GEDAI=1`
- ERPLAB, if `RELAX_cfg.NotchFilterType='PMnotch'`

zapline-plus is **no longer** an external dependency when `RELAX_cfg.NotchFilterType='ZaplinePlus'`: the
complete upstream tree (`f4628a3`) is vendored under `vendor/zapline-plus/`, with 4 files patched (3
carrying local boundary customizations, 1 an unrelated carried-forward fix) plus one new file
(`pwelch_segmented.m`) with no upstream counterpart — see
`vendor/zapline-plus/VENDORING.md`. Entry points must call `RELAX_add_vendored_zapline()` after
`eeglab()` so this directory is prepended to the path ahead of the EEGLAB plugin's copy; without that
call, MATLAB resolves the plugin's copy instead, which is boundary-aware but still has the nfft-grid
bug (see "Boundary/discontinuity handling" below). Leaving the zapline-plus
EEGLAB plugin installed is harmless and is what keeps zapline-plus available on EEGLAB's own menus.

## Pipeline architecture (`RELAX_Wrapper.m` / `RELAX_Wrapper_slurm.m` / `RELAX_Wrapper_bart.m`)

All three wrapper variants follow the same stage sequence over a per-file loop
(`for FileNumber = RELAX_cfg.FilesToProcess`), driven entirely by fields on the `RELAX_cfg` struct
built in the corresponding `*_SET_PARAMETERS_AND_RUN*.m` script. Each stage stamps its own metrics
onto `EEG.RELAXProcessing*` substructs, which get folded into per-participant rows of aggregate
tables (`RELAXProcessingMWFStepOneAllParticipants`, `CleanedMetrics`, `RawMetrics`, etc.) and
periodically checkpointed to `.mat` files in `RELAX_cfg.OutputPath` (`<myPath>/RELAXProcessed/`) —
these are reloaded at the start of a run so a run can resume/accumulate across invocations.

Stage order:

1. **Load & prep** — `pop_loadset`, optional channel-location lookup, delete irrelevant channels,
   split off "keep but don't clean" auxiliary electrodes (`RELAX_cfg.electrodes_2_keep_but_not_clean`).
2. **Notch filter** — Butterworth, ZaplinePlus, or ERPLAB PMnotch, selected by
   `RELAX_cfg.NotchFilterType`.
3. **Band-pass filter & optional downsample** — `RELAX_filtbutter.m` or `pop_eegfiltnew`; order
   relative to MWF controlled by `RELAX_cfg.LowPassFilterBeforeMWF` (low-pass-before-MWF increases
   rank-deficiency risk and is not the recommended default).
4. **Bad-channel rejection** — PREP's `findNoisyChannels` (off by default in the SLURM path via
   `RELAX_cfg.Perform_PREP_bad_channel_rejection=0`; see "Boundary/discontinuity handling" below).
5. **Extreme-period marking** — `RELAX_excluding_channels_and_epoching.m` (or the `_iTBS` variant for
   TMS data with pulse artifacts) then `RELAX_excluding_extreme_values.m`; blink detection via
   `RELAX_blinks_IQR_method.m`.
6. **Artifact-type marking** — `RELAX_muscle.m`, `RELAX_horizontaleye.m`, `RELAX_drift.m` build masks
   of artifact vs. clean data (padded via `RELAX_pad_brief_mask_periods.m`) used both for the MWF
   cleaning template and for post-hoc SER/ARR (signal-to-error/artifact-to-residue ratio) metrics.
7. **MWF cleaning** — `RELAX_perform_MWF_cleaning.m`, run up to three times
   (`RELAX_cfg.Do_MWF_Once/Twice/Thrice`), each round re-marking artifacts on the residual.
8. **Average re-reference** — `RELAX_average_rereference.m`, after extreme-outlier periods are
   rejected.
9. **Residual-artifact ICA cleaning** — mutually-exclusive options selected by `RELAX_cfg` flags:
   `Perform_targeted_wICA` (recommended; `RELAX_targeted_wICA.m`), `Perform_wICA_on_ICLabel`
   (`RELAX_wICA_on_ICLabel_artifacts.m`), or `Perform_ICA_subtract` (`RELAX_ICA_subtract.m`) — all
   use ICLabel component classification, thresholded by `RELAX_cfg.ICLabel_thresholds`.
10. **GEDAI lead-field filtering** (this fork's addition) — gated by `RELAX_cfg.Run_GEDAI`. Rejected
    channels are interpolated first, then `GEDAI(EEG, ...)` is called (external GEDAI toolbox, not in
    this repo); optionally with a per-participant reference covariance matrix loaded by
    `loadGEDAI_ppt_refCOV.m` (matched to the EEG filename by its `<id>_<site>` prefix). GEDAI's
    per-file metrics (SENSAI score, ENOVA, epochs rejected) are folded into `CleanedMetrics`.
11. **Compute cleaned metrics, interpolate remaining missing channels, save** — cleaned data goes to
    `RELAXProcessed/Cleaned_Data/`, intermediate rounds to `RELAXProcessed/<N>xMWF/` (if
    `RELAX_cfg.saveround*`-style flags are set), GEDAI-rejected-artifact copies to
    `RELAXProcessed/GEDAI_artifacts/`, and aggregate stats/metrics `.mat` files directly under
    `RELAXProcessed/`.

After cleaning, `RELAX_epoch_the_clean_data_Wrapper.m` / `pop_RELAX_epoch_the_clean_data.m` /
`RELAX_EPOCH_CLEAN_DATA_FOR_ANALYSIS.m` handle epoching the cleaned continuous data for downstream
analysis (interpolating channels, re-coding triggers, rejecting still-artifactual epochs).

### Naming convention for pipeline variants

Suffixed files (`_slurm`, `_bart`, `_iTBS`) are dataset/execution-mode-specific copies of a base
file with the same stage structure — not newer/older versions of each other. When fixing a bug that
affects the shared algorithm, check whether the same fix is needed in the sibling variants
(`RELAX_Wrapper.m` / `RELAX_Wrapper_slurm.m` / `RELAX_Wrapper_bart.m`, and
`RELAX_excluding_channels_and_epoching.m` / `_iTBS.m`) rather than assuming one is canonical.

`vendor/` follows a different convention: `vendor/zapline-plus/` mirrors its upstream repository
byte-for-byte except where `vendor/zapline-plus/VENDORING.md` documents a local customization — so an
upstream re-sync stays diffable against the vendored commit rather than against a hand-edited copy.

### Skip-already-processed behavior

Both the plain and SLURM entry-point scripts scan `RELAXProcessed/Cleaned_Data/` for filenames
already containing `RELAX` matching an input file's basename, and drop already-processed files from
`RELAX_cfg.files`/`RELAX_cfg.dirList` before the run — reruns of an entry-point script are therefore
incremental by default, scoped per-participant in the SLURM variant.

### Boundary/discontinuity handling

The TMS-EEG inputs to `RELAX_Wrapper_slurm.m` contain `boundary` events marking periods excised by
TMS artifact removal, leaving only ~8 s of continuous data between them. `RELAX_cfg.RespectBoundaries=1`
makes the pipeline process each inter-boundary segment independently, with odd reflection padding at
each segment edge. Design and analysis:
`docs/superpowers/specs/2026-09-08-boundary-aware-preprocessing-design.md`.

Helpers: `RELAX_get_segments.m` (segment list from `EEG.event`), `RELAX_reflect_pad.m` (odd reflection),
`RELAX_apply_by_segments.m` (apply a function per segment), `RELAX_check_boundary_support.m` (startup
guard — see caveat below). Callers must always re-derive the segment list from the *current* `EEG` —
`eeg_eegrej` adds boundaries part-way through the wrapper.

`RELAX_check_boundary_support` runs once, at startup, on the freshly loaded file — it only sees
discontinuities already present in the INPUT file. It is not re-checked later in the run, and
`eeg_eegrej` adds more boundaries mid-pipeline (rejecting extreme periods) on essentially every file,
so a clean pass at startup does not mean every later stage actually saw and handled every boundary a
given run ends up with — it means the input file didn't already contain a combination the guard knows
is unsafe.

Boundary-aware without any change from us (verified against EEGLAB 2026.0.0 source): `pop_resample`
(splits at boundaries), `clean_data_with_zapline_plus_eeglab_wrapper` (passes `boundaries` to zapline),
`pop_eegfiltnew`/`firfilt` (splits at discontinuities), and `pop_epoch` (drops any epoch containing a
boundary — which is what protects `RELAX_epoching`, `RELAX_excluding_channels_and_epoching`,
`RELAX_excluding_extreme_values`, `RELAX_muscle` and `RELAX_drift`). `RELAX_ICA_subtract` is also safe:
it is a purely spatial projection.

The zapline-plus boundary-awareness above holds **only for the vendored copy** under
`vendor/zapline-plus/`, reached via `RELAX_add_vendored_zapline()` (see "External dependencies" above).
Its boundary customizations (C1-C6, `vendor/zapline-plus/VENDORING.md`) were originally authored
directly against the installed EEGLAB plugin copy, so that copy also passes `boundaries` to zapline —
but it still has the frequency-grid-collapse crash on short inter-boundary segments that
`docs/superpowers/plans/2026-09-09-zapline-nfft-grid-fix.md` fixed only in the vendored tree. If
`RELAX_add_vendored_zapline()` is not called before the notch-filter stage, MATLAB can resolve
`clean_data_with_zapline_plus_eeglab_wrapper` from the plugin path instead, which is boundary-aware but
not crash-fixed.

The vendored copy's **adaptive chunker** (C8, `vendor/zapline-plus/VENDORING.md`) is also boundary-aligned:
covariance epochs never straddle a boundary, an inter-boundary segment shorter than `minChunkLength/2`
is never subdivided, chunk edges either fall on a boundary or leave at least `minChunkLength/2` on both
sides, and the narrowband detection filter runs per segment. This too holds only for the vendored copy
reached via `RELAX_add_vendored_zapline()`; the EEGLAB plugin copy remains sample-1-anchored.

Made boundary-aware here: `RELAX_filtbutter.m` (optional 7th argument; this includes the Butterworth
notch — its bandstop branch pads by the stop-band width rather than the centre frequency, since a
notch rings at its bandwidth), `RELAX_blinks_IQR_method.m`, `RELAX_targeted_wICA.m`.

**Still boundary-blind — `RELAX_check_boundary_support` errors if any of these is enabled on data with
boundaries already present in the input file:** PREP `findNoisyChannels`
(`Perform_PREP_bad_channel_rejection`, off in the shipped entry-point config — the wrapper itself
defaults an absent flag to on, for backward compatibility — it tiles correlation windows flatly from
sample 1 against a 1% `badTimeThreshold`, so ~9-11% of windows straddle a discontinuity on this data),
MWF delay-stacking (`Do_MWF_*`), `Run_GEDAI`, `Perform_wICA_on_ICLabel`, and the SER/ARR marking block
(`ComputeSER_ARR`, off by default — its only consumer `RELAX_metrics_final_SER_and_ARR` is commented
out, and `RELAX_muscle`, `RELAX_horizontaleye`, `RELAX_drift` and `RELAX_pad_brief_mask_periods` would
all need fixing first).

Tests live in `tests/`. Run them with:
`matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"`
