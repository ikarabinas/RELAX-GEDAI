# CLAUDE.md

## What this is

RELAX-GEDAI is a fork of the RELAX EEG artifact-cleaning pipeline (MATLAB/EEGLAB). It is a
script-driven toolbox, not a packaged library: no build step, no package manager, no automated
correctness gate beyond `tests/` — real correctness is verified by running the pipeline on real EEG
data and inspecting the cleaning-quality metrics it outputs.

This is a lab research codebase used by multiple people (originally Neil Bailey's RELAX, extended by
Isabella Karabinas for GEDAI, further adapted here for TMS-EEG / SLURM use). Dataset-specific
hardcoded paths inside scripts (`RELAX_cfg.myPath`, `addpath(...)` to EEGLAB/plugins) are meant to be
edited per-run/per-user, not treated as bugs. When adapting a script for a new user/dataset, edit the
paths at the top rather than adding config indirection.

`~/TMS_response_variability` is a separate repo holding the scientific experiment that consumes this
pipeline's output: this repo is the preprocessing stage, that one is the analysis layer.

## Running the pipeline

No CLI entry point; these are MATLAB scripts (`run(...)` or the editor's Run button):

- `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` / `RELAX_Wrapper_slurm.m` — **the current path.** SLURM
  job-array variant: groups files by participant, selects one participant per array task from
  `SLURM_ARRAY_TASK_ID`, and gives each task its own `parpool` `JobStorageLocation` under `$TMPDIR`.
  Submit with `sbatch run_relax-gedai_slurm.sh` (adjust `--array=1-N%M` to the participant count).
  The index → participant mapping is pinned in `<myPath>/participant_list.mat`
  (`RELAX_participant_list.m`): existing participants never change index, later arrivals are
  appended; the entry script prints the mapping, so widen `--array` when the count grows.
- `RELAX_SET_PARAMETERS_AND_RUN.m` — deprecated single-process entry point; sets addpaths, builds
  `RELAX_cfg` and the file list, calls `RELAX_Wrapper(RELAX_cfg)`. `run_relax-gedai.sh` submits it.
- `convert_files_to_eeglab_format/` — standalone conversion of raw EGI `.mff` / `.fif` to EEGLAB
  `.set` (`save_mff_to_set*.m`, `save_fif_to_set.m`; drive with `python run_set_conversion.py` after
  `module load matlab`).
- `find_missing_files.m` — diff input `.set` files against `RELAXProcessed/Cleaned_Data` output.
- `recover_metrics.m` — rebuild aggregate metric `.mat`/`.csv` from individual processed files after
  an interrupted run.
- `RELAX_QC_figures.m` (`sbatch run_RELAX_QC_figures.sh [extract|draw|all] [file,...]`) — lab-facing
  artifact-removal QC figure set from a finished `RELAX_Wrapper_slurm.m` run: example traces, spectra,
  paired blink/muscle metrics, data taken per file, blink-locked topography, into
  `RELAXProcessed/QC_figures/` with `captions.md`. Stage 1 caches per-file summaries (`qc_cache.mat`);
  stage 2 redraws from it in seconds. Its "before" is the wICA input rebuilt from `Extremes_Rejected`
  (`RELAX_QC_baseline.m`), so it refuses configs where that is wrong (MWF/GEDAI on, low-pass after the
  save) — `RELAX_QC_check_config.m`. Design: `docs/superpowers/specs/2026-10-15-qc-figures-design.md`.
- `RELAX_residual_screen.m` (`sbatch run_RELAX_residual_screen.sh [score|sheets|all] [file,...]`) — decides
  which cleaned files are clean enough for source localization. `score` measures muscle proportion and
  residual blink ratio per file and flags robust-z outliers over the cleaned files (muscle strength is kept as
  a covariate only), and records why each interpolated channel was removed (`RELAX_screen_channel_reasons.m`); it also
  fails a file whose interpolation hole is deeper than `holeDepthLimit` (1.4× the median electrode spacing
  from the nearest real electrode, `RELAX_screen_hole_depth.m`) or whose interpolated-channel count reached
  the pipeline's `MaxProportionOfElectrodesThatCanBeDeleted` limit. It writes `RELAXProcessed/Residual_screen/residual_screen.mat/.csv` with a `README.md`;
  `sheets` draws a review sheet into `review/` for each failing, unscreened or sampled-passing file (or just
  the named files); its muscle windows are the epochs with the highest pipeline muscle slope
  (`RELAX_screen_muscle_slopes.m`, via `RELAX_metrics_muscle`). Settings are at the top of the script.
- `validate_pop_repair.m` (`sbatch run_validate_pop_repair.sh`) — real-data validation of the
  electrode-pop stage on the four Sep26 surveyed files (m155 was dropped as too noisy). Superseded
  derivative-detector tooling is parked in
  `variance_ratio_survey/legacy_tooling/`.

Tests: `matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"`

### External dependencies (not vendored, expected on the MATLAB path)

- EEGLAB (with PREP, ICLabel, PICARD, FastICA plugins)
- MATLAB Signal Processing and Statistics and Machine Learning Toolboxes
- MWF artifact removal toolbox (github.com/exporl/mwf-artifact-removal)
- FieldTrip
- GEDAI (github.com/neurotuning/GEDAI-master or ikarabinas/GEDAI-master) when `RELAX_cfg.Run_GEDAI=1`
- ERPLAB, if `RELAX_cfg.NotchFilterType='PMnotch'`
- EEGLAB `ebridge` plugin, if `RELAX_cfg.Perform_bridged_channel_rejection=1`

Vendored, so **not** dependencies:

- `vendor/zapline-plus/` — full upstream tree (`f4628a3`), 4 files patched + `pwelch_segmented.m`
  added; see `VENDORING.md` there. Entry points must call `RELAX_add_vendored_zapline()` after
  `eeglab()` to prepend it ahead of the EEGLAB plugin's copy, which is boundary-aware but still has
  the nfft-grid bug. Leaving the plugin installed is harmless.
- `vendor/noisetools/` — deliberately partial: just the 9-file runtime closure of `nt_detrend`, which
  `RELAX_repair_steps.m` calls. Nothing shadows it, but entry points must still call
  `RELAX_add_vendored_noisetools()` after `eeglab()`.

## Pipeline architecture (`RELAX_Wrapper_slurm.m`, plus deprecated `RELAX_Wrapper.m` / `RELAX_Wrapper_bart.m`)

All three wrappers share the stage sequence below over a per-file loop
(`for FileNumber = RELAX_cfg.FilesToProcess`), driven by `RELAX_cfg` fields set in the matching
`*_SET_PARAMETERS_AND_RUN*.m`. Each stage stamps metrics onto `EEG.RELAXProcessing*` substructs,
folded into per-participant rows of aggregate tables (`RELAXProcessingMWFStepOneAllParticipants`,
`CleanedMetrics`, `RawMetrics`, …) and periodically checkpointed to `.mat` in
`RELAX_cfg.OutputPath` (`<myPath>/RELAXProcessed/`); these reload at the start of a run so runs
resume/accumulate.

1. **Load & prep** — `pop_loadset`, optional channel-location lookup, delete irrelevant channels,
   split off "keep but don't clean" auxiliaries (`RELAX_cfg.electrodes_2_keep_but_not_clean`).
2. **Impedance screen, electrode-pop detection & repair, bridge screen** (`RELAX_Wrapper_slurm.m`
   only) — see "Stage 2 in the SLURM wrapper" below.
3. **Notch filter** — Butterworth, ZaplinePlus, or ERPLAB PMnotch per `RELAX_cfg.NotchFilterType`.
   Butterworth/PMnotch run before band-pass & downsampling; ZaplinePlus runs after downsampling.
4. **Band-pass filter & optional downsample** — `RELAX_filtbutter.m` or `pop_eegfiltnew`; order
   relative to MWF set by `RELAX_cfg.LowPassFilterBeforeMWF` (low-pass-before-MWF raises
   rank-deficiency risk; not the recommended default).
5. **PREP bad-channel rejection** — `findNoisyChannels`, `RELAX_cfg.Perform_PREP_bad_channel_rejection`
   (off by default in the SLURM path; boundary-blind). In `RELAX_Wrapper_slurm.m` the impedance and
   bridge screens have moved to stage 2, so only PREP remains here. **PREP does not spend from the
   shared channel-rejection budget**, so enabling it lets total deletions across all three screens
   exceed `RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted`.
6. **Extreme-period marking** — `RELAX_excluding_channels_and_epoching.m` (or `_iTBS` for TMS pulse
   artifacts, or `_boundary`, used automatically when `RELAX_cfg.RespectBoundaries=1`) then
   `RELAX_excluding_extreme_values.m`; blink detection via `RELAX_blinks_IQR_method.m`. If fewer than
   two `RELAX_cfg.BlinkElectrodes` survive the channel screens, the blink method falls back to the
   row one higher on the net (`E12 E20 E27 E34 E21`), and failing that reports no blinks through the
   existing no-blink path rather than erroring. The choice is a covariate (`BlinkElectrodesUsed`,
   `BlinkElectrodeFallback` 0/1/2 on `RELAXProcessingExtremeRejections`) because a file cleaned
   without blink marking is cleaned differently.
   In `RELAX_Wrapper_slurm.m`, `RELAX_cfg.RejectMuscleChannelsAfterCleaning=1` switches the epoching
   stage's muscle-channel screen off and reruns it after GEDAI (or after ICA when GEDAI is off) on the
   cleaned data, with the same thresholds and the shared budget (`RELAX_muscle_channels_to_reject.m`),
   so muscle-dominated electrodes are available to ICA and GEDAI. Covariates:
   `RejectMuscleChannelsAfterCleaning` and `NumberOfMuscleContaminatedChannelsBeforeCleaning` on
   `RELAXProcessingExtremeRejections`, where `NumberOfMuscleContaminatedChannelsRecomendedToDelete` then
   counts the cleaned data. Comparison tooling: `muscle_deferral_comparison/`.
   `RELAX_cfg.ExtremeMinBadChannelsToRejectEpoch=k` (absent/1 = RELAX's rule) rejects an epoch in
   `RELAX_excluding_extreme_values.m` only when at least k channels trip a per-channel criterion
   (voltage shift, absolute voltage, joint probability, drift); whole-montage kurtosis still rejects
   alone. On the ~200-channel net the any-one-channel rule mostly condemned epochs on one channel
   (`extreme_rejection_survey/summary.txt`). Single-channel kurtosis is computed but has never rejected
   anything, upstream included, and still doesn't. Covariates on `RELAXProcessingExtremeRejections`:
   `ExtremeMinBadChannelsToRejectEpoch`, `ProportionExcludedIfAnyChannelRejected`,
   `EpochsKeptWithFlaggedChannels`, `FlaggedChannelEpochsKept` (flagged channel-epochs passed on to ICA/GEDAI).
   `RELAX_cfg.InterpolateExtremeChannelEpochsAfterCleaning=1` (`RELAX_Wrapper_slurm.m` only, absent/0 =
   off) repairs the flagged channel-epochs those kept epochs carry into the cleaned data. Right after
   the first pass, `RELAX_lock_extreme_thresholds.m` fixes each channel's statistics from all
   first-pass epochs, average-referenced: the shift bound without the blink-aware bound, the
   joint-probability histogram, and the kurtosis mean and SD. After the deferred muscle screen and
   before the cleaned metrics, `RELAX_second_pass_channel_epochs.m` re-epochs the cleaned data the same
   way and flags channel-epochs by shift, absolute voltage, joint probability (improbable side only),
   drift and single-channel kurtosis (z > 10). It replaces each flagged channel over its whole 1 s
   epoch with a spherical spline of the epoch's other channels, blended over 0.25 s inside each edge.
   An epoch with more than `SecondPassMaxProportionFlaggedChannels` (0.20) of its channels flagged
   raises `RELAX:secondPassTooManyFlaggedChannels`. Covariates: `InterpolateExtremeChannelEpochsAfterCleaning`,
   `SecondPassEpochsEvaluated`, `SecondPassEpochsWithInterpolatedChannels`,
   `SecondPassChannelEpochsInterpolated`, `SecondPassProportionChannelTimeInterpolated`,
   `SecondPassFlaggedBy<ShiftHigh|ShiftLow|AbsVolt|JointProb|Drift|Kurtosis>` (NaN while off). Design:
   `docs/superpowers/specs/2026-11-06-post-cleaning-channel-epoch-interpolation-design.md`.
7. **Artifact-type marking** — `RELAX_muscle.m`, `RELAX_horizontaleye.m`, `RELAX_drift.m` build
   artifact/clean masks (padded by `RELAX_pad_brief_mask_periods.m`) used both for the MWF template
   and for post-hoc SER/ARR metrics.
8. **MWF cleaning** — `RELAX_perform_MWF_cleaning.m`, up to three rounds
   (`RELAX_cfg.Do_MWF_Once/Twice/Thrice`), re-marking artifacts on the residual each round.
9. **Average re-reference** — `RELAX_average_rereference.m`, after extreme-outlier period rejection.
10. **Residual-artifact ICA cleaning** — mutually exclusive: `Perform_targeted_wICA` (recommended,
    `RELAX_targeted_wICA.m`), `Perform_wICA_on_ICLabel`, or `Perform_ICA_subtract`; all use ICLabel
    classification thresholded by `RELAX_cfg.ICLabel_thresholds`.
11. **GEDAI lead-field filtering** (fork addition) — gated by `RELAX_cfg.Run_GEDAI`, parameters in
    `RELAX_cfg.GEDAI_*` (defaults match `~/TMS_response_variability/iTBS_pipeline.m`; the GEDAI fork at
    `~/gitRepos/GEDAI-master` is added to the path by the entry script). In `RELAX_Wrapper_slurm.m`
    GEDAI runs on the surviving channels only; refCOV is built from the participant's raw leadfield
    (`RELAX_load_ppt_leadfield.m`, `<ppt>_leadfield.mat` in `RELAX_cfg.GEDAI_LeadfieldPath`), rows
    selected by channel label and average-referenced over that subset. The deprecated wrappers still
    interpolate first and use `loadGEDAI_ppt_refCOV.m`. GEDAI metrics (SENSAI, ENOVA, epochs
    rejected) fold into `CleanedMetrics`.
12. **Cleaned metrics, interpolate remaining missing channels, save** — `RELAXProcessed/Cleaned_Data/`,
    intermediate rounds to `RELAXProcessed/<N>xMWF/`, GEDAI-rejected artifacts to
    `RELAXProcessed/GEDAI_artifacts/`, aggregates directly under `RELAXProcessed/`.

After cleaning, `RELAX_epoch_the_clean_data_Wrapper.m` / `pop_RELAX_epoch_the_clean_data.m` /
`RELAX_EPOCH_CLEAN_DATA_FOR_ANALYSIS.m` epoch the cleaned continuous data for analysis (interpolating
channels, re-coding triggers, rejecting still-artifactual epochs).

### Stage 2 in the SLURM wrapper

Runs at the native sample rate, before any filtering — a zero-phase highpass turns a DC step into a
multi-second bipolar excursion in both directions. Order within the stage: impedance screen → pop
stage → bridge screen. Impedance runs first so neither the pop detector's whole-montage statistics
nor PREP's correlation statistics ever see a channel already known bad
(`RELAX_cfg.Perform_impedance_bad_channel_rejection`, `RELAX_impedance_bad_channels.m`; spends from
the shared channel budget via its optional 4th argument `maxRejectableCount`). The bridge screen
(`RELAX_cfg.Perform_bridged_channel_rejection`, `RELAX_bridged_channels.m`, needs `ebridge`) runs
last and does not spend from the budget. A recording whose bridge screen reports `excessiveBridging`
(more than `RELAX_cfg.BridgeMaxProportionChannels`) is **excluded**, not stripped and cleaned: the
wrapper raises `RELAX:excludeRecording`.

Per-file `try/catch` around the file loop records both exclusions and errors in
`RELAXProcessed/PerParticipantStats/RELAX_unprocessed_files_<participant>.mat/.csv` (columns
`File, Outcome, Identifier, Message, When`; `Outcome` = `excluded` or `error`) and moves to the next
file, so one bad file no longer kills the array task and every file missing from `Cleaned_Data/` has
a row saying why. A file that succeeds on rerun has its row removed.

The **pop stage** itself is gated by `RELAX_cfg.Perform_pop_artifact_repair` (absent/0 = off; new
behaviour must never switch itself on for a dataset never checked for pops). Its sub-steps
(reference pops, detection, eye masking, channel ranking, repair), parameters, covariates and design
documents are in the `pop-artifact-stage` skill (`.claude/skills/pop-artifact-stage/SKILL.md`).

### Naming convention for pipeline variants

Suffixed files (`_slurm`, `_bart`, `_iTBS`, `_boundary`) are dataset/execution-mode-specific copies of
a base file with the same stage structure — not newer/older versions of each other. When fixing a bug
in shared algorithm code, check whether the sibling variants need the same fix
(`RELAX_Wrapper.m` / `_slurm.m` / `_bart.m`, and `RELAX_excluding_channels_and_epoching.m` /
`_iTBS.m` / `_boundary.m`) rather than assuming one is canonical. `_boundary.m`'s channel-rejection
logic (everything after epoch placement) was copied from `_iTBS.m` rather than factored into a shared
helper, so a fix there likely needs porting to both.

All of stage 2's fork-specific machinery — pop detection/repair and the shared channel-rejection
budget — exists **only** in `RELAX_Wrapper_slurm.m`. In `RELAX_Wrapper.m` / `RELAX_Wrapper_bart.m` the
impedance/bridge/PREP screens still cap themselves independently at
`MaxProportionOfElectrodesThatCanBeDeleted` each. Porting means adding the budget init, the stage
block, and the impedance screen's 4th argument to those wrappers.

`vendor/` follows a different convention: `vendor/zapline-plus/` mirrors upstream byte-for-byte except
where its `VENDORING.md` documents a local customization, so a re-sync stays diffable against the
vendored commit. `vendor/noisetools/` is a deliberately partial mirror.

### Boundary/discontinuity handling

The TMS-EEG inputs contain `boundary` events marking periods excised by TMS artifact removal, leaving
~8 s of continuous data between them. `RELAX_cfg.RespectBoundaries=1` processes each inter-boundary
segment independently, with odd reflection padding at each edge. Design:
`docs/superpowers/specs/2026-09-08-boundary-aware-preprocessing-design.md`.

Helpers: `RELAX_get_segments.m` (segment list from `EEG.event`), `RELAX_reflect_pad.m`,
`RELAX_apply_by_segments.m`, `RELAX_check_boundary_support.m` (startup guard on the input file only —
`eeg_eegrej` adds boundaries mid-pipeline that it never re-checks). Callers must always re-derive the
segment list from the *current* `EEG`, never a cached one.

Boundary-aware already: `pop_resample`, `pop_eegfiltnew`/`firfilt`, `pop_epoch` (and everything it
gates: `RELAX_excluding_channels_and_epoching`, `RELAX_excluding_extreme_values`, `RELAX_muscle`,
`RELAX_drift`), `RELAX_ICA_subtract`, `RELAX_bridged_channels.m` (by construction — its epoching goes
through `RELAX_get_segments`), and zapline-plus **only via the vendored copy**. Note `pop_epoch`'s
protection costs something: it epochs on a fixed global 0.5 s grid and silently drops
boundary-straddling epochs, leaving a sub-1 s sliver un-epoched on both sides of every boundary
(stays `NaN`, later deleted by `eeg_eegrej` whether or not it was extreme) — hence `_boundary.m`.

Made boundary-aware in this fork: `RELAX_filtbutter.m` (optional 7th argument; covers the Butterworth
notch), `RELAX_blinks_IQR_method.m`, `RELAX_targeted_wICA.m`,
`RELAX_excluding_channels_and_epoching_boundary.m` (epochs each segment independently, tiling epoch
placement evenly across it instead of relying on a global grid; used automatically when
`RespectBoundaries=1`).

**Still boundary-blind — `RELAX_check_boundary_support` errors if any is enabled on input data with
boundaries:** PREP `findNoisyChannels` (`Perform_PREP_bad_channel_rejection`), MWF delay-stacking
(`Do_MWF_*`), `Run_GEDAI`, `Perform_wICA_on_ICLabel`, and the SER/ARR marking block (`ComputeSER_ARR`).

### Skip-already-processed behavior

Both entry-point scripts scan `RELAXProcessed/Cleaned_Data/` for filenames containing `RELAX` that
match an input basename and drop those files from `RELAX_cfg.files`/`RELAX_cfg.dirList` before the
run, so reruns are incremental by default (per-participant in the SLURM variant). Set
`RELAX_cfg.OverwriteExistingCleanedFiles=1` in `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` to reprocess and
silently overwrite; defaults to `0`. The plain entry-point script has no such flag.

### Saving covariates for downstream ML classification tasks

Nearly every RELAX threshold is computed from within-recording statistics, so the same `RELAX_cfg`
cleans a noisy session harder than a quiet one, in a way that tracks that session's noise level. A
model predicting treatment response from session-to-session EEG change can then learn the cleaning
rather than the neural signal. Any new stage must therefore record per-file covariates describing how
hard it cleaned (see `RELAXProcessing_wICA`, `RELAXProcessingExtremeRejections`).
