# Vendored zapline-plus

Upstream: https://github.com/MariusKlug/zapline-plus
Vendored at: f4628a3 (2023-04-14, HEAD at time of vendoring)
Customizations originally authored against: fb0dcbb (2022-05-19) == the installed
EEGLAB plugin zapline-plus1.2.1 (22 of its 26 .m files are byte-identical to it;
the 4 listed below are not).

This is a COMPLETE mirror of the upstream tree, not a partial one. Every file is
upstream's, byte for byte, except:

  - the 4 customized files listed below,
  - this file,
  - `pwelch_segmented.m`, a file with no upstream counterpart at all (see "Changes
    made by the nfft-grid-fix plan" below), and `zapline_boundary_segments.m`,
    `zapline_covariance_epochs.m`, `zapline_chunk_indices.m` and
    `zapline_filter_by_segments.m`, likewise with no upstream counterpart (see
    "Changes made by the boundary-aligned-chunking plan" below).

(As of this commit, `pwelch_segmented` has been moved out of the bottom of
`clean_data_with_zapline_plus.m`, where it started as a local subfunction, into
that file of its own - Task 2 of this plan.)

`zapline_filter_by_segments.m` calls `RELAX_reflect_pad` from the repo root, so the
vendored tree is not standalone: it depends on the RELAX-GEDAI repo layout, not just
on MATLAB/EEGLAB toolboxes. It is only ever reached through `RELAX_add_vendored_zapline()`,
which is called from entry points that already have the repo root on the path.

`RELAX_add_vendored_zapline.m` prepends this directory to the MATLAB path so these
copies shadow the EEGLAB plugin's. Because the mirror is complete, RELAX no longer
depends on the plugin being installed - though leaving it installed is harmless and
is what keeps zapline-plus available from EEGLAB's own menus. `eegplugin_zapline_plus.m`
is mirrored for completeness but is inert here: EEGLAB only runs eegplugin_* files it
finds by scanning its own plugins folder, not ones sitting on the path.

Upstream is GPL-3, the same license as RELAX; `LICENSE` is mirrored as that requires.

To re-sync with a newer upstream, re-run the procedure in
`docs/superpowers/plans/2026-09-09-zapline-nfft-grid-fix.md`, Task 1 Step 1.

## Local customizations

| ID | File | Change |
|----|------|--------|
| C1 | `..._eeglab_wrapper.m` | Read `boundary` event latencies from `EEG.event` into `zaplineConfig.boundaries`; harden the `zaplineConfig` existence test with `|| isempty(...)`. |
| C2 | `clean_data_with_zapline_plus.m` | Add the `boundaries` parameter (parser, local, and `zaplineConfig` record). |
| C3 | `clean_data_with_zapline_plus.m` | Comment out the "data set is short" `winSizeCompleteSpectrum` auto-adjust. `winSizeCompleteSpectrum` is unused once C4 replaces every `pwelch` call. |
| C4 | `clean_data_with_zapline_plus.m` | Replace 5 `pwelch` calls with `pwelch_segmented`, which never lets an FFT window straddle a discontinuity. |
| C5 | `clean_data_with_zapline_plus.m` | Slice the boundary list per chunk and rebase it to chunk-local indices. |
| C6 | `nt_zapline_plus.m` | Apply `nt_smooth` independently within each inter-boundary segment. |
| C7 | `nt_smooth.m` | `filter` -> `filtfilt` (zero-phase smoothing). Not boundary-related; carried forward unchanged. |
| C8 | `clean_data_with_zapline_plus.m` | Adaptive-chunk block: replaced the sample-1-anchored narrowband filter, covariance tiling and chunk-edge arithmetic with boundary-aligned equivalents built on `zapline_boundary_segments`, `zapline_covariance_epochs`, `zapline_chunk_indices` and `zapline_filter_by_segments`. Chunk edges are now the first sample of the following epoch (upstream used the last sample of the current one). `'MinPeakDistance', minChunkLength` was removed from `findpeaks` because it counts epoch indices; `minChunkLength` is now enforced in samples inside `zapline_chunk_indices`. Also removed upstream's commented-out `%% plot` block (dead code plotting the narrowband covariance distances). Also added `numel(distances) >= 3` and non-empty-`proms` guards around the `findpeaks` calls - a genuine behaviour change: upstream errored outright when fewer than 3 covariance epochs were available, whereas this now degrades to a single whole-record chunk instead. The narrowband detection filter now also receives odd reflection padding at the *record* edges (not just at internal boundaries), so `data_narrowfilt` differs from upstream in the first and last ~`filtord` samples even on data with no boundaries at all - that perturbs the first and last epoch covariances, hence `distances`, hence potentially the selected peak set. |
| C9 | `clean_data_with_zapline_plus.m` | Analytics: added `analyticsResults.chunkIndices`, a cell of the chunk layout per noise frequency. |

## Changes made by the nfft-grid-fix plan

(`docs/superpowers/plans/2026-09-09-zapline-nfft-grid-fix.md`)

- `pwelch_segmented` moved out of `clean_data_with_zapline_plus.m` into its own
  file so it can be unit-tested (`tests/test_pwelch_segmented.m`).
- `pwelch_segmented` no longer aliases pwelch's `nfft` argument with its own
  frequency-vector output.
- All 5 `pwelch_segmented` call sites in `clean_data_with_zapline_plus.m` pass an
  explicit `specNfft`, computed once from `detailedFreqBoundsUpper` (pinning the
  frequency grid to a bin spacing fine enough for the fine noise-peak search to
  always find bins, independent of how short the data's inter-boundary segments
  are), instead of `[]`.

## Changes made by the boundary-aligned-chunking plan

(`docs/superpowers/plans/2026-09-10-boundary-aligned-zapline-chunking.md`, C8-C9 above)

- Four new files with no upstream counterpart: `zapline_boundary_segments.m`,
  `zapline_covariance_epochs.m`, `zapline_chunk_indices.m` and
  `zapline_filter_by_segments.m`. See each file's own header comment for its
  contract.
- `clean_data_with_zapline_plus.m`'s adaptive-chunk block rewritten to use them
  (C8), and `analyticsResults.chunkIndices` added to record the resulting chunk
  layout per noise frequency (C9).

## Inherited upstream changes

Rebasing the 4 customizations from their `fb0dcbb` baseline (what the EEGLAB plugin
runs) onto `f4628a3` pulls in upstream's own changes to the customized files across
that range, beyond the customizations themselves. Both are noted here so that a
difference between this pipeline's output and a pre-vendoring run is not mistakenly
attributed to the nfft-grid fix in Tasks 2-4 when it is really one of these.

**`cleanData = data;` moved to after flat-channel removal (upstream commit
1038c89) - a real behaviour change, kept deliberately.** In `fb0dcbb`,
`cleanData = data;` sits at line 278, before the flat-channel strip at lines
285-290. With that ordering, the per-chunk write-back
(`cleanData(chunkIndices(iChunk):...,:) = cleanData_chunk;`, line 555) assigns a
narrowed (flat-channels-removed) chunk into `cleanData`, which is still full width,
while the flat-channel reinsertion later (line 949) assumes `cleanData` already
excludes them - a latent shape mismatch on any file with flat channels. Upstream
`1038c89` moved `cleanData = data;` to after the flat-channel strip, which is now
where this vendored copy has it too (see U1 in Task 1's verification steps). This is
a genuine upstream bugfix, not a cosmetic reordering, and we are deliberately
keeping it rather than reverting to match `fb0dcbb`: RELAX runs ZaplinePlus before
any bad-channel rejection, so flat channels reaching zapline is a realistic case on
this data, not a hypothetical one. **If cleaned output on flat-channel data differs
between a pre-vendoring run and this pipeline, attribute that difference to this
hunk - not to the nfft-grid fix Tasks 2-4 make. Task 5's before/after comparison
should account for it.**

**`noisefreqs='line'` search widened (upstream commit 18369dc) - dead code here,
left as-is.** Upstream commit 18369dc widened the `noisefreqs='line'` search but its
code reads `idx = (f > 49 & f < 50) | (f > 59 & f < 60);` while its comment states
the intent was 49-51 and 59-61 Hz. The range is one-sided and excludes 60.0 Hz
itself. RELAX never sets `noisefreqs='line'` (it leaves it empty for automatic
detection), so this branch is dead here. Left as upstream wrote it.
