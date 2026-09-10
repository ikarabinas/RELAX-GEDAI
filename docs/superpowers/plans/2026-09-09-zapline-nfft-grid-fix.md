# Zapline-plus Frequency-Grid Fix (vendored, rebased on upstream HEAD)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop `clean_data_with_zapline_plus` from dying at the fine noise-peak search on boundary-heavy TMS-EEG data, by pinning the pwelch frequency grid instead of inheriting it from the shortest inter-boundary segment — delivered as a vendored copy of zapline-plus rebased onto current upstream.

**Architecture:** Vendor the complete zapline-plus tree (32 files: 26 `.m`, `LICENSE`, upstream `README.md`, 4 PNGs) into `vendor/zapline-plus/` at upstream HEAD (2023-04-14), with the four locally-customized files rebased from the 2022-05-19 baseline they were authored against. The vendored directory is prepended to the MATLAB path by the RELAX entry points so it shadows the EEGLAB plugin copy. `pwelch_segmented` is promoted from a file-local subfunction to its own file so it becomes unit-testable, then given two fixes: stop aliasing pwelch's `nfft` argument with its frequency-vector output, and accept an explicit `nfft` derived from `detailedFreqBoundsUpper`.

**Tech Stack:** MATLAB (Signal Processing Toolbox), EEGLAB 2026.0.0, zapline-plus 1.2.1, SLURM.

**Spec:** No separate spec document. The root-cause analysis is inlined in "Background" below. Related design doc: `docs/superpowers/specs/2026-09-08-boundary-aware-preprocessing-design.md`.

## Global Constraints

- MATLAB only; no build step or package manager. Run the suite with:
  `matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"`
- **Never run MATLAB on a login node** (`scu-login*`). Wrap every MATLAB invocation in
  `srun --partition=scu-cpu --cpus-per-task=2 --mem=8G --time=00:20:00 --chdir=$PWD bash -c 'module load matlab; ...'`.
- Vendored upstream baseline: `https://github.com/MariusKlug/zapline-plus` commit `f4628a3` (2023-04-14, current HEAD).
- Local customizations were authored against commit `fb0dcbb` (2022-05-19). Exactly 4 of the installed plugin's 26 `.m` files differ from `fb0dcbb`; the other 22 are byte-identical, so that commit is the confirmed baseline.
- Vendor the **whole** upstream tree, not just the customized files: all 26 `.m` files plus `LICENSE`, upstream's `README.md`, and the 4 example PNGs. zapline-plus is GPL-3 (same as RELAX), so `LICENSE` must travel with the redistributed copy. The PNGs account for ~3 MB of the 3.4 MB total; they are upstream README illustrations and are included for a faithful mirror.
- Upstream ships its own `README.md`. Keep it untouched and record the vendoring provenance in `vendor/zapline-plus/VENDORING.md` instead, so the mirror stays diffable against upstream.
- Do **not** change the defaults of `detailedFreqBoundsUpper` (`[-0.05 0.05]`) or `detailedFreqBoundsLower` (`[-0.4 0.1]`). The fix adapts the frequency grid to them, not the reverse.
- Do **not** modify anything under `~/tools/eeglab2026.0.0/`. The installed plugin stays exactly as it is; the vendored copy shadows it.
- **Out of scope** (deliberately deferred): decoupling the pwelch *window length* from the shortest inter-boundary segment. That is a separate change; see "Deferred" at the end.

---

## Background: the root cause

`RELAX_Wrapper_slurm.m:327` calls zapline after downsampling to 500 Hz. On `m208_dlpfc_day1` the input carries 60 `boundary` events leaving 61 segments (min 7.930 s, median 7.931 s, max 147.0 s). A 7.93 s segment at 500 Hz is 3965 samples.

`pwelch_segmented` sizes its window from the *shortest* segment in the chunk, so `window = hanning(3965)`. MATLAB's `pwelch` then picks `nfft = max(256, 2^nextpow2(3965)) = 4096`, giving a bin spacing of `500/4096 = 0.12207 Hz`. The bins bracketing 60 Hz land at 59.9365 and 60.0586 — both outside the `±0.05 Hz` fine-peak window. So:

```matlab
this_freq_idx_detailed = f>noisefreq-0.05 & f<noisefreq+0.05;   % all false
this_freqs_detailed    = f(this_freq_idx_detailed);             % 0x1
maxFinePower           = max(mean(pxx_chunk(this_freq_idx_detailed,:),2));  % []
noisePeaks(iChunk)     = this_freqs_detailed(... == maxFinePower);          % 0 elements -> ERROR
```

producing *"Unable to perform assignment because the left and right sides have a different number of elements."*

Upstream is unaffected because it uses `pwelch(chunk, hanning(length(chunk)), ...)` — a window of the whole chunk (≥ `minChunkLength` = 30 s), giving `nfft = 16384` and 3 bins inside the window.

A second, independent defect lives in the same function: `f` is used as both the 4th input argument (pwelch's `nfft` slot) and the output frequency vector. After iteration 1 it holds a 2049-element vector, so iterations 2..N pass a *vector* where `nfft` belongs, which switches pwelch to its Goertzel path. Measured ratio of Goertzel to FFT output: **0.5 at every interior frequency**, 1.0 at DC and Nyquist — the one-sided power doubling is not applied. Segment 1 therefore contributes twice the weight of every other segment to `pxx_total`.

**Chosen grid:** target a bin spacing of `min(abs(detailedFreqBoundsUpper))/1.5` = 0.0333 Hz, i.e. `nfft = 2^nextpow2(srate/0.0333)` = **16384** at 500 Hz. This reproduces upstream's effective grid exactly (upstream: `2^nextpow2(30*500)` = 16384), puts 3 bins in the search window, and lands within 0.003 dB of the peak-height asymptote. Measured alternatives: 8192 gives 1 bin and is 0.16 dB low; 32768 and 65536 give no measurable improvement at 2× and 4× the memory.

Zero-padding sets the frequency *grid*, not the resolution bandwidth — that stays governed by the 3965-sample window (~0.18 Hz Hann −3 dB). That is acceptable: measured peak-location error from an 8 s window zero-padded to 16384 was 0.003–0.010 Hz for lines at 60.031 and 59.977 Hz, well inside the search window.

---

## File Structure

**Created:**
- `vendor/zapline-plus/` — full mirror of upstream `f4628a3` (32 files). Of these:
  - 22 `.m` files copied verbatim (`nt_pca.m`, `nt_dss0.m`, `find_next_noisefreq.m`, `pop_zapline_plus.m`, `eegplugin_zapline_plus.m`, …), plus `LICENSE`, upstream `README.md`, and 4 PNGs.
  - `clean_data_with_zapline_plus.m` — upstream + customizations C2–C5 + this plan's `specNfft` threading.
  - `clean_data_with_zapline_plus_eeglab_wrapper.m` — upstream + customization C1.
  - `nt_zapline_plus.m` — upstream + customization C6.
  - `nt_smooth.m` — upstream + customization C7.
- `vendor/zapline-plus/pwelch_segmented.m` — new; extracted from the subfunction, then fixed. The only `.m` file with no upstream counterpart.
- `vendor/zapline-plus/VENDORING.md` — records the upstream baseline, the rebase target, and each local customization.
- `RELAX_add_vendored_zapline.m` — prepends the vendor dir to the path and asserts it wins precedence.
- `tests/test_pwelch_segmented.m`
- `tests/test_zapline_boundary_crash.m`
- `tests/test_RELAX_add_vendored_zapline.m`

**Modified:**
- `RELAX_SET_PARAMETERS_AND_RUN_slurm.m` — call `RELAX_add_vendored_zapline` after `eeglab;`.
- `RELAX_SET_PARAMETERS_AND_RUN.m`, `RELAX_SET_PARAMETERS_AND_RUN_bart.m`, `pop_RELAX.m`, `AARATEP_RELAX_tmseeg_pipeline.m` — same one-line call (sibling-variant rule from `CLAUDE.md`).
- `CLAUDE.md` — document the vendored plugin and the path requirement.

---

### Task 1: Vendor the full zapline-plus tree and make it win path precedence

The rebase has been verified to apply with zero conflicts (only line offsets). Upstream changed one file in two hunks — moving `cleanData = data;` after flat-channel removal, and a 2D rewrite of the `noisefreqs='line'` candidate search — and neither touches a customized region.

**Files:**
- Create: `vendor/zapline-plus/` — the complete upstream tree at `f4628a3` (32 files), with 4 of them rebased-and-patched
- Create: `vendor/zapline-plus/VENDORING.md`
- Create: `RELAX_add_vendored_zapline.m`
- Create: `tests/test_RELAX_add_vendored_zapline.m`
- Modify: `RELAX_SET_PARAMETERS_AND_RUN_slurm.m:116` (after `eeglab;`)

**Interfaces:**
- Produces: `RELAX_add_vendored_zapline()` — no arguments, no return value; errors if the vendored copy does not resolve first.

- [ ] **Step 1: Produce the rebased files**

```bash
cd ~/.claude_scratch
rm -rf zapline-upstream && git clone https://github.com/MariusKlug/zapline-plus.git zapline-upstream
cd zapline-upstream
L=/home/nga4004/tools/eeglab2026.0.0/plugins/zapline-plus1.2.1
V=~/gitRepos/RELAX-GEDAI/vendor/zapline-plus
PATCHED="clean_data_with_zapline_plus.m clean_data_with_zapline_plus_eeglab_wrapper.m nt_zapline_plus.m nt_smooth.m"

# Confirm the baseline: exactly these 4 .m files may differ from fb0dcbb
git checkout --quiet fb0dcbb
for f in *.m; do cmp -s "$f" "$L/$f" || echo "DIFFERS: $f"; done
# Expected output: exactly the 4 files in $PATCHED, nothing else.

# Capture the customizations as patches against the baseline they were authored on
mkdir -p /tmp/zp_patches
for f in $PATCHED; do diff -u "$f" "$L/$f" > /tmp/zp_patches/$f.patch; done

# Mirror the WHOLE upstream tree at HEAD, excluding git metadata
git checkout --quiet f4628a3
mkdir -p "$V"
git archive f4628a3 | tar -x -C "$V"

# Rebase the 4 customized files onto it
cd "$V"
for f in $PATCHED; do patch --no-backup-if-mismatch -p0 "$f" < /tmp/zp_patches/$f.patch; done
```

Expected: every file patches successfully; `clean_data_with_zapline_plus.m` reports hunks succeeding at offsets (1 and 42 lines). **Any "FAILED" or "Hunk ... with fuzz" output means upstream moved since this plan was written — stop and re-derive, do not force.**

- [ ] **Step 2: Verify the rebase kept both sides**

```bash
cd ~/gitRepos/RELAX-GEDAI/vendor/zapline-plus

# The mirror must be complete and, apart from the 4 patched files, byte-identical to upstream
ls | wc -l                                                             # expect 32
diff <( (cd ~/.claude_scratch/zapline-upstream && git ls-tree -r --name-only f4628a3) | sort ) \
     <(ls | sort)                                                      # expect no output
for f in *.m; do
  cmp -s "$f" ~/.claude_scratch/zapline-upstream/"$f" || echo "patched: $f"
done                                                                   # expect exactly the 4 customized files

grep -n -B2 '^cleanData = data;' clean_data_with_zapline_plus.m       # upstream U1: must sit AFTER flat-channel removal
grep -n 'spectraChunk_allChans' clean_data_with_zapline_plus.m         # upstream U2: 2D candidate search present
grep -c 'pwelch_segmented(' clean_data_with_zapline_plus.m             # customization C4: expect 5 calls + 1 definition = 6
grep -n 'chunkBoundaries =' clean_data_with_zapline_plus.m             # customization C5
grep -n "addOptional(p, 'boundaries'" clean_data_with_zapline_plus.m   # customization C2
grep -n 'filtfilt' nt_smooth.m                                          # customization C7
```

- [ ] **Step 3: Write `vendor/zapline-plus/VENDORING.md`**

Do not overwrite upstream's own `README.md` — this is a separate file alongside it.

```markdown
# Vendored zapline-plus

Upstream: https://github.com/MariusKlug/zapline-plus
Vendored at: f4628a3 (2023-04-14, HEAD at time of vendoring)
Customizations originally authored against: fb0dcbb (2022-05-19) == the installed
EEGLAB plugin zapline-plus1.2.1 (22 of its 26 .m files are byte-identical to it;
the 4 listed below are not).

This is a COMPLETE mirror of the upstream tree, not a partial one. Every file is
upstream's, byte for byte, except:

  - the 4 customized files listed below,
  - `pwelch_segmented.m`, which has no upstream counterpart,
  - this file.

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
| C1 | `..._eeglab_wrapper.m` | Read `boundary` event latencies from `EEG.event` into `zaplineConfig.boundaries`; harden the `zaplineConfig` existence test with `\|\| isempty(...)`. |
| C2 | `clean_data_with_zapline_plus.m` | Add the `boundaries` parameter (parser, local, and `zaplineConfig` record). |
| C3 | `clean_data_with_zapline_plus.m` | Comment out the "data set is short" `winSizeCompleteSpectrum` auto-adjust. `winSizeCompleteSpectrum` is unused once C4 replaces every `pwelch` call. |
| C4 | `clean_data_with_zapline_plus.m` | Replace 5 `pwelch` calls with `pwelch_segmented`, which never lets an FFT window straddle a discontinuity. |
| C5 | `clean_data_with_zapline_plus.m` | Slice the boundary list per chunk and rebase it to chunk-local indices. |
| C6 | `nt_zapline_plus.m` | Apply `nt_smooth` independently within each inter-boundary segment. |
| C7 | `nt_smooth.m` | `filter` -> `filtfilt` (zero-phase smoothing). Not boundary-related; carried forward unchanged. |

## Changes made by this plan

- `pwelch_segmented` moved out of `clean_data_with_zapline_plus.m` into its own
  file so it can be unit-tested (`tests/test_pwelch_segmented.m`).
- `pwelch_segmented` no longer aliases pwelch's `nfft` argument with its own
  frequency-vector output.
- All 5 call sites pass an explicit `specNfft` derived from
  `detailedFreqBoundsUpper`.

## Inherited upstream quirk

Upstream commit 18369dc widened the `noisefreqs='line'` search but its code reads
`idx = (f > 49 & f < 50) | (f > 59 & f < 60);` while its comment states the intent
was 49-51 and 59-61 Hz. The range is one-sided and excludes 60.0 Hz itself. RELAX
never sets `noisefreqs='line'` (it leaves it empty for automatic detection), so
this branch is dead here. Left as upstream wrote it.
```

- [ ] **Step 4: Write the failing path-precedence test**

```matlab
function test_RELAX_add_vendored_zapline()
% The vendored zapline-plus copies must shadow the EEGLAB plugin's.

    RELAX_add_vendored_zapline();

    vendorDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'vendor', 'zapline-plus');

    % The mirror is complete, so EVERY zapline-plus function must resolve from it -
    % not just the customized ones. A partial shadow is the failure mode that would
    % pair a patched clean_data_with_zapline_plus with the plugin's unpatched nt_smooth.
    mFiles = dir(fullfile(vendorDir, '*.m'));
    % 26 upstream .m files. Task 2 adds pwelch_segmented.m, making it 27; the
    % floor stays at 26 so this test is correct both before and after that.
    assert(numel(mFiles) >= 26, ...
        'expected the full upstream tree (26 .m files), found %d', numel(mFiles));
    for k = 1:numel(mFiles)
        [~, fname] = fileparts(mFiles(k).name);
        resolved = which(fname);
        assert(strncmp(resolved, vendorDir, numel(vendorDir)), ...
            '%s resolved to %s, expected the vendored copy in %s', fname, resolved, vendorDir);
    end

    % Calling it twice must be harmless
    RELAX_add_vendored_zapline();
    assert(strncmp(which('nt_smooth'), vendorDir, numel(vendorDir)), 'must be idempotent');

    fprintf('test_RELAX_add_vendored_zapline: PASS\n');
end
```

- [ ] **Step 5: Run it and watch it fail**

```bash
srun --partition=scu-cpu --cpus-per-task=2 --mem=8G --time=00:20:00 --chdir=$PWD \
  bash -c 'module load matlab; matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,\"tests\")); test_RELAX_add_vendored_zapline"'
```

Expected: FAIL with `Unrecognized function or variable 'RELAX_add_vendored_zapline'`.

- [ ] **Step 6: Implement `RELAX_add_vendored_zapline.m`**

```matlab
function RELAX_add_vendored_zapline()
% RELAX_ADD_VENDORED_ZAPLINE - put the vendored zapline-plus ahead of the EEGLAB plugin.
%
% The EEGLAB plugin copy of zapline-plus is boundary-blind and crashes on data whose
% inter-boundary segments are short (see
% docs/superpowers/plans/2026-09-09-zapline-nfft-grid-fix.md). vendor/zapline-plus
% holds rebased, fixed copies of the four files that needed changing. Call this AFTER
% eeglab(), because eeglab() prepends its own plugin directories to the path.

    vendorDir = fullfile(fileparts(mfilename('fullpath')), 'vendor', 'zapline-plus');
    if ~isfolder(vendorDir)
        error('RELAX_add_vendored_zapline:missing', 'Vendored zapline-plus not found at %s', vendorDir);
    end

    addpath(vendorDir);  % addpath prepends, so this wins over the plugin directory

    resolved = which('clean_data_with_zapline_plus');
    if ~strncmp(resolved, vendorDir, numel(vendorDir))
        error('RELAX_add_vendored_zapline:shadowed', ...
            ['clean_data_with_zapline_plus resolves to\n  %s\ninstead of the vendored copy in\n  %s\n' ...
             'Something re-prepended the EEGLAB plugin path after this call - move the ' ...
             'RELAX_add_vendored_zapline() call to after the last eeglab() invocation.'], resolved, vendorDir);
    end
end
```

- [ ] **Step 7: Run the test and confirm it passes**

Same command as Step 5. Expected: `test_RELAX_add_vendored_zapline: PASS`.

- [ ] **Step 8: Wire it into the entry points**

In `RELAX_SET_PARAMETERS_AND_RUN_slurm.m`, immediately after line 144's `addpath('/home/nga4004/gitRepos/RELAX-GEDAI/');`:

```matlab
% Put the vendored, boundary-fixed zapline-plus ahead of the EEGLAB plugin copy.
% Must come after eeglab(), which prepends its own plugin directories.
RELAX_add_vendored_zapline();
```

Apply the identical two lines to `RELAX_SET_PARAMETERS_AND_RUN.m`, `RELAX_SET_PARAMETERS_AND_RUN_bart.m`, `pop_RELAX.m`, and `AARATEP_RELAX_tmseeg_pipeline.m`, each after that script's RELAX-GEDAI `addpath`. (Sibling-variant rule, `CLAUDE.md`.)

- [ ] **Step 9: Commit**

```bash
git add vendor/zapline-plus RELAX_add_vendored_zapline.m tests/test_RELAX_add_vendored_zapline.m \
        RELAX_SET_PARAMETERS_AND_RUN*.m pop_RELAX.m AARATEP_RELAX_tmseeg_pipeline.m
git commit -m "Vendor the full zapline-plus tree at upstream f4628a3, shadowing the EEGLAB plugin"
```

---

### Task 2: Extract `pwelch_segmented` to its own file and pin both bugs with failing tests

**Files:**
- Create: `vendor/zapline-plus/pwelch_segmented.m`
- Modify: `vendor/zapline-plus/clean_data_with_zapline_plus.m` (delete the subfunction and its now-superfluous `end`)
- Create: `tests/test_pwelch_segmented.m`

**Interfaces:**
- Consumes: `RELAX_add_vendored_zapline()` from Task 1.
- Produces: `[pxx, f] = pwelch_segmented(data, window, noverlap, nfft, fs, boundaries)` — `data` is samples×channels, `nfft` is a **scalar or `[]`** (never a frequency vector), `boundaries` is a vector of sample indices or `[]`.

- [ ] **Step 1: Move the subfunction into its own file, byte-for-byte unchanged**

Cut the `function [pxx, f] = pwelch_segmented(...) ... end` block from the end of `vendor/zapline-plus/clean_data_with_zapline_plus.m` into `vendor/zapline-plus/pwelch_segmented.m`. Then delete the bare `end` that C4 added to terminate the main function (a file whose only function is the main one must not be `end`-terminated if the original wasn't — check with `mlint`).

- [ ] **Step 2: Confirm the extraction did not break the file**

```bash
srun --partition=scu-cpu --cpus-per-task=2 --mem=8G --time=00:20:00 --chdir=$PWD \
  bash -c 'module load matlab; matlab -nodisplay -nosplash -batch "addpath(pwd); RELAX_add_vendored_zapline; checkcode(which(\"clean_data_with_zapline_plus\")); checkcode(which(\"pwelch_segmented\"))"'
```

Expected: no syntax errors (style warnings are fine).

- [ ] **Step 3: Write the failing tests**

```matlab
function test_pwelch_segmented()
% Boundary-aware Welch PSD used by the vendored zapline-plus.

    RELAX_add_vendored_zapline();

    fs     = 250;
    segLen = 1983;   % 7.93 s - the m208_dlpfc_day1 inter-train segment length
    nChan  = 4;

    % ---------------------------------------------------------------------
    % 1. Every segment must contribute with equal weight.
    % Three identical segments must average to the PSD of any one of them.
    % The pre-fix code assigned pwelch's frequency-vector output back into the
    % variable it passed as pwelch's nfft argument, so segments 2..N ran through
    % pwelch's Goertzel path, which omits the one-sided power doubling and
    % returns exactly half the PSD. That made the average 2/3 of the truth.
    % ---------------------------------------------------------------------
    rng(0);
    one    = randn(segLen, nChan);
    data   = [one; one; one];
    bounds = [segLen + 1; 2*segLen + 1];

    [pxxN, ~] = pwelch_segmented(data, [], [], 8192, fs, bounds);
    [pxx1, ~] = pwelch(one, hanning(segLen), [], 8192, fs);
    relErr = max(abs(pxxN(:) - pxx1(:))) / max(pxx1(:));
    assert(relErr < 1e-10, ...
        'identical segments must average to one segment''s PSD (rel err %.4g)', relErr);

    % ---------------------------------------------------------------------
    % 2. An explicit nfft must give the +/-0.05 Hz fine-peak window real bins.
    % ---------------------------------------------------------------------
    nfft = 2^nextpow2(round(fs / (0.05/1.5)));   % 8192 at fs = 250
    [~, f] = pwelch_segmented(data, [], [], nfft, fs, bounds);
    assert(sum(f > 59.95 & f < 60.05) >= 1, ...
        'the +/-0.05 Hz fine-peak window must contain at least one bin');

    % ---------------------------------------------------------------------
    % 3. Regression witness: omitting nfft reintroduces the crash condition.
    % The window-derived grid on 7.93 s segments has bins at 59.937 and 60.059
    % and none in between. Kept so nobody "simplifies" the call sites back.
    % ---------------------------------------------------------------------
    [~, fDefault] = pwelch_segmented(data, [], [], [], fs, bounds);
    assert(sum(fDefault > 59.95 & fDefault < 60.05) == 0, ...
        'expected the window-derived grid to be too coarse - this is what passing nfft fixes');

    % ---------------------------------------------------------------------
    % 4. nfft below the window length must be raised, never wrap the data.
    % ---------------------------------------------------------------------
    [~, fSmall] = pwelch_segmented(data, [], [], 256, fs, bounds);
    dfSmall = fSmall(2) - fSmall(1);
    assert(dfSmall <= fs / 2^nextpow2(segLen) + eps, ...
        'nfft smaller than the window must be raised so pwelch zero-pads instead of wrapping');

    % ---------------------------------------------------------------------
    % 5. Continuous data must behave exactly like a plain full-window pwelch.
    % ---------------------------------------------------------------------
    cont     = randn(4000, nChan);
    [pc, ~]  = pwelch_segmented(cont, [], [], [], fs, []);
    [pr, ~]  = pwelch(cont, hanning(4000), [], 2^nextpow2(4000), fs);
    assert(isequal(size(pc), size(pr)) && max(abs(pc(:) - pr(:))) < 1e-12, ...
        'boundary-free data must match a plain full-window pwelch');

    % ---------------------------------------------------------------------
    % 6. Known limitation, characterised deliberately: the window still
    % collapses to the SHORTEST segment, so a short fragment at a chunk edge
    % degrades the whole chunk's spectrum. Deferred - see the plan's "Deferred"
    % section. Update this test when that is addressed.
    % ---------------------------------------------------------------------
    ragged  = [randn(50, nChan); one; one];
    rbounds = [51; 51 + segLen];
    [~, fRagged] = pwelch_segmented(ragged, [], [], [], fs, rbounds);
    assert((fRagged(2) - fRagged(1)) > 0.5, ...
        'documenting that a 50-sample fragment still collapses the window');

    fprintf('test_pwelch_segmented: PASS\n');
end
```

- [ ] **Step 4: Run and confirm the right tests fail**

```bash
srun --partition=scu-cpu --cpus-per-task=2 --mem=8G --time=00:20:00 --chdir=$PWD \
  bash -c 'module load matlab; matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,\"tests\")); test_pwelch_segmented"'
```

Expected: FAIL on assertion 1 with a relative error near 0.33. (Assertions 2–6 may pass already; only 1 is guaranteed red at this point, because the pre-fix code does honour a scalar `nfft` on its first iteration.)

- [ ] **Step 5: Commit the red state**

```bash
git add vendor/zapline-plus/pwelch_segmented.m vendor/zapline-plus/clean_data_with_zapline_plus.m tests/test_pwelch_segmented.m
git commit -m "Extract pwelch_segmented to its own file; add failing tests for the nfft aliasing and grid bugs"
```

---

### Task 3: Fix the `nfft` aliasing and accept an explicit `nfft`

**Files:**
- Modify: `vendor/zapline-plus/pwelch_segmented.m`
- Test: `tests/test_pwelch_segmented.m` (unchanged from Task 2)

**Interfaces:**
- Produces: the signature from Task 2, now honoured — `nfft` is read-only input, `f` is output-only.

- [ ] **Step 1: Rewrite `vendor/zapline-plus/pwelch_segmented.m`**

```matlab
function [pxx, f] = pwelch_segmented(data, window, noverlap, nfft, fs, boundaries)
% PWELCH_SEGMENTED - Welch PSD that never lets an FFT window straddle a discontinuity.
%
% Averages pwelch estimates computed independently within each stretch of continuous
% data delimited by BOUNDARIES, so the broadband splatter of a discontinuity never
% enters the spectrum.
%
%   data       - samples x channels
%   window     - window vector, or [] to derive one from the segment layout
%   noverlap   - passed through to pwelch ([] = pwelch's default 50%)
%   nfft       - FFT length, scalar or []. MUST NOT be a frequency vector: pwelch
%                treats a vector here as a request for its Goertzel algorithm, which
%                omits the one-sided power doubling and so returns half the PSD.
%                nfft sets the frequency GRID (spacing fs/nfft), which is independent
%                of the WINDOW length that sets the true resolution bandwidth. Callers
%                searching a narrow frequency window must pass an nfft fine enough
%                that the window actually contains bins.
%   fs         - sampling rate in Hz
%   boundaries - sample indices of discontinuities, or [] for continuous data
%
% RELAX-GEDAI local addition to zapline-plus. See vendor/zapline-plus/VENDORING.md.

    if isempty(window)
        if isempty(boundaries)
            window = hanning(size(data, 1));
        else
            bounds_vec = boundaries(:)';
            segment_lengths = diff(unique(sort([1, bounds_vec, size(data, 1) + 1])));
            max_win_length = max(1, min(segment_lengths));
            window = hanning(max_win_length);
        end
    end

    % pwelch wraps the data modulo nfft when nfft is shorter than the window, which
    % silently coarsens the grid instead of zero-padding. Never allow that.
    if isempty(nfft)
        nfft = 2^nextpow2(length(window));
    else
        nfft = max(nfft, 2^nextpow2(length(window)));
    end

    if isempty(boundaries) || length(boundaries) < 2
        [pxx, f] = pwelch(data, window, noverlap, nfft, fs);
        return;
    end

    boundaries = unique(sort([1; boundaries(:); (size(data, 1) + 1)]));
    pxx_total = 0;
    total_segments = 0;
    f = [];

    for i = 1:length(boundaries) - 1
        start_idx = boundaries(i);
        end_idx   = boundaries(i + 1) - 1;
        seg_len   = end_idx - start_idx + 1;

        if seg_len >= length(window)
            % NOTE: nfft is passed in unchanged every iteration, and the frequency
            % vector comes back into a separate variable. Assigning f into the nfft
            % argument would switch pwelch to Goertzel from the second segment on.
            [pxx_seg, f] = pwelch(data(start_idx:end_idx, :), window, noverlap, nfft, fs);
            pxx_total = pxx_total + pxx_seg;
            total_segments = total_segments + 1;
        end
    end

    if total_segments > 0
        pxx = pxx_total / total_segments;
    else
        [pxx, f] = pwelch(data, window, noverlap, nfft, fs);
    end
end
```

- [ ] **Step 2: Run the tests and confirm they pass**

```bash
srun --partition=scu-cpu --cpus-per-task=2 --mem=8G --time=00:20:00 --chdir=$PWD \
  bash -c 'module load matlab; matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,\"tests\")); test_pwelch_segmented"'
```

Expected: `test_pwelch_segmented: PASS`.

- [ ] **Step 3: Run the whole suite for regressions**

```bash
srun --partition=scu-cpu --cpus-per-task=2 --mem=16G --time=00:40:00 --chdir=$PWD \
  bash -c 'module load matlab; matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,\"tests\")); run_all_tests"'
```

Expected: all test files pass.

- [ ] **Step 4: Commit**

```bash
git add vendor/zapline-plus/pwelch_segmented.m
git commit -m "Fix pwelch_segmented: stop aliasing nfft with the frequency output, honour explicit nfft"
```

---

### Task 4: Thread `specNfft` through the five call sites

**Files:**
- Modify: `vendor/zapline-plus/clean_data_with_zapline_plus.m` (5 call sites + 1 new computation)
- Create: `tests/test_zapline_boundary_crash.m`

**Interfaces:**
- Consumes: `pwelch_segmented(..., nfft, ...)` from Task 3.

- [ ] **Step 1: Write the end-to-end failing test**

```matlab
function test_zapline_boundary_crash()
% Regression for the reported m208_dlpfc_day1 failure: clean_data_with_zapline_plus
% died in its fine noise-peak search on data whose inter-boundary segments are ~8 s,
% because the frequency grid derived from those segments had no bin inside the
% +/-0.05 Hz search window.
%
%   Unable to perform assignment because the left and right sides have a
%   different number of elements.
%   Error in clean_data_with_zapline_plus (line 492)

    RELAX_add_vendored_zapline();

    fs     = 250;
    segLen = round(7.93 * fs);   % same 0.12207 Hz grid as the 500 Hz production case
    nSeg   = 20;                 % ~159 s, enough for several >=30 s adaptive chunks
    nChan  = 8;

    rng(7);
    t      = (0:segLen-1)' / fs;
    line   = 5 * sin(2*pi*60.02*t);
    data   = zeros(nSeg*segLen, nChan);
    bounds = zeros(nSeg-1, 1);
    for k = 1:nSeg
        rows = (k-1)*segLen + (1:segLen);
        data(rows, :) = repmat(line, 1, nChan) + randn(segLen, nChan);
        if k > 1
            bounds(k-1) = rows(1);
        end
    end

    cfg = struct('plotResults', 0, 'boundaries', bounds);
    [cleaned, cfgOut, analytics] = clean_data_with_zapline_plus(data, fs, cfg);

    assert(isequal(size(cleaned), size(data)), 'output size must match input');
    assert(all(isfinite(cleaned(:))), 'cleaned data must be finite');

    % The 60 Hz line must actually be attenuated
    [pRaw,   f] = pwelch(data(1:segLen,:),    hanning(segLen), [], 8192, fs);
    [pClean, ~] = pwelch(cleaned(1:segLen,:), hanning(segLen), [], 8192, fs);
    band = f > 59.5 & f < 60.5;
    attenDb = 10*log10(mean(mean(pRaw(band,:))) / mean(mean(pClean(band,:))));
    assert(attenDb > 6, 'expected >6 dB attenuation at 60 Hz, got %.2f dB', attenDb);

    fprintf('test_zapline_boundary_crash: PASS\n');
end
```

- [ ] **Step 2: Run it and confirm it reproduces the reported crash**

```bash
srun --partition=scu-cpu --cpus-per-task=2 --mem=8G --time=00:20:00 --chdir=$PWD \
  bash -c 'module load matlab; matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,\"tests\")); test_zapline_boundary_crash"'
```

Expected: FAIL with *"Unable to perform assignment because the left and right sides have a different number of elements"* inside `clean_data_with_zapline_plus`. **If it fails with any other message, stop — the fixture is not reproducing the real bug.**

- [ ] **Step 3: Compute `specNfft`**

In `vendor/zapline-plus/clean_data_with_zapline_plus.m`, immediately before the `disp('Computing initial spectrum...')` line:

```matlab
% Frequency grid for every spectrum that feeds the fine noise-peak search below.
% That search indexes f with detailedFreqBoundsUpper (default +/-0.05 Hz around the
% noise frequency). pwelch's default nfft follows the WINDOW length, and
% pwelch_segmented sizes its window from the shortest stretch of continuous data - on
% boundary-heavy recordings that can be short enough that the bin spacing exceeds the
% search window, leaving it with no bins at all and making the peak search return
% empty. Pin the grid to a target bin spacing instead, so the window always holds
% several bins however the data is segmented. This sets the GRID, not the resolution
% bandwidth, which remains governed by the window length.
specTargetDf = min(abs(detailedFreqBoundsUpper)) / 1.5;
specNfft     = 2^nextpow2(round(srate / specTargetDf));   % 16384 at srate 500
```

- [ ] **Step 4: Pass it at all five call sites**

Replace the 4th argument `[]` with `specNfft` in each:

```matlab
[pxx_raw_log,f]=pwelch_segmented(data,[],[],specNfft,srate, boundaries);
...
                [pxx_chunk,f]=pwelch_segmented(chunk,[],[],specNfft,srate, this_zaplineConfig_chunk.boundaries);
...
        [pxx_raw]=pwelch_segmented(data,[],[],specNfft,srate, boundaries);
        [pxx_clean,f]=pwelch_segmented(cleanData,[],[],specNfft,srate, boundaries);
        [pxx_removed]=pwelch_segmented(data-cleanData,[],[],specNfft,srate, boundaries);
```

Verify all five were changed: `grep -c 'pwelch_segmented(.*specNfft' vendor/zapline-plus/clean_data_with_zapline_plus.m` must print `5`.

- [ ] **Step 5: Run the test and confirm it passes**

Same command as Step 2. Expected: `test_zapline_boundary_crash: PASS`.

- [ ] **Step 6: Confirm the NaN metric is fixed too**

`this_freq_idx_noise` (the `proportionRemovedNoise` computation) uses the same `±0.05 Hz` window against the whole-data spectrum, so it was returning `NaN`. Add to the end of `test_zapline_boundary_crash`, before the `fprintf`:

```matlab
    assert(cfgOut.detailedFreqBoundsUpper(2) == 0.05, 'guard: default fine bounds assumed');
    assert(all(isfinite(analytics.proportionRemovedNoise)), ...
        'proportionRemovedNoise must be finite - it was NaN when the +/-0.05 Hz index was empty');
```

`cfgOut` and `analytics` are already captured by the call in Step 1 - do not invoke
`clean_data_with_zapline_plus` a second time. Re-run. Expected: PASS.

- [ ] **Step 7: Run the whole suite**

```bash
srun --partition=scu-cpu --cpus-per-task=2 --mem=16G --time=00:40:00 --chdir=$PWD \
  bash -c 'module load matlab; matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,\"tests\")); run_all_tests"'
```

- [ ] **Step 8: Commit**

```bash
git add vendor/zapline-plus/clean_data_with_zapline_plus.m tests/test_zapline_boundary_crash.m
git commit -m "Pin the zapline spectrum grid to detailedFreqBoundsUpper instead of the segment length"
```

---

### Task 5: Verify on the real failing file

**Files:** none modified — this is the acceptance gate.

- [ ] **Step 1: Find m208's array index**

```bash
srun --partition=scu-cpu --cpus-per-task=1 --mem=4G --time=00:10:00 --chdir=$PWD \
  bash -c 'module load matlab; matlab -nodisplay -nosplash -batch "L=load(\"/home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Aug26/tms_artifact_removed/participant_list.mat\"); idx=find(strcmp(L.uniqueParticipants,\"m208\")); fprintf(\"m208 is array task %d of %d\n\", idx, numel(L.uniqueParticipants))"'
```

Note the index — do **not** rely on the `task_id = 7` interactive fallback at `RELAX_SET_PARAMETERS_AND_RUN_slurm.m:192`.

- [ ] **Step 2: Confirm the day1 output is absent so the skip logic won't skip it**

```bash
ls /home/nga4004/TMS_response_variability/store/mdd_dlpfc/preprocessed_Aug26/tms_artifact_removed/RELAXProcessed/Cleaned_Data/ | grep m208_dlpfc_day1 || echo "not yet processed - good"
```

- [ ] **Step 3: Submit only that participant**

```bash
sbatch --array=<idx> run_relax-gedai_slurm.sh
```

- [ ] **Step 4: Confirm the zapline stage cleared**

Watch the job's stdout. Expected: the `Removing noise at 60Hz...`, `N chunks will be created.` and `Cleaning with ZapLine-plus done!` lines appear with **no** "different number of elements" error, and the run proceeds past `RELAX_Wrapper_slurm.m:327` into bad-channel rejection.

- [ ] **Step 5: Sanity-check the cleaning metrics**

Once the file completes, confirm `EEG.etc.zapline.analyticsResults.proportionRemovedNoise` is finite and `NremoveFinal` is neither all-zero nor pinned at the maximum for every chunk.

- [ ] **Step 6: Update `CLAUDE.md`**

Add to the "External dependencies" section: zapline-plus is no longer an external dependency — the complete upstream tree (`f4628a3`) is vendored under `vendor/zapline-plus/`, with 4 files carrying local boundary customizations and one new file (`pwelch_segmented.m`). Entry points must call `RELAX_add_vendored_zapline()` after `eeglab()`; without it MATLAB resolves the EEGLAB plugin's boundary-blind copies instead. Leaving the plugin installed is harmless and keeps zapline-plus on EEGLAB's menus. Add to the "Boundary/discontinuity handling" section that zapline-plus is boundary-aware **via the vendored copy only**. Note in the naming-convention section that `vendor/` mirrors upstream byte-for-byte except where `VENDORING.md` says otherwise, so upstream re-syncs stay diffable.

- [ ] **Step 7: Commit**

```bash
git add CLAUDE.md
git commit -m "Document the vendored zapline-plus and its path requirement"
```

---

## Flagged for your decision

Two things surfaced during the rebase that this plan carries forward **unchanged**, because they are yours and not part of the reported bug. Say so if you want either revisited:

1. **`nt_smooth.m`: `filter` -> `filtfilt` (C7).** This is not boundary-related. It makes zapline's bias smoother zero-phase, which changes what `x - xx` contains and therefore what the DSS step sees — a real behavioural divergence from published zapline-plus. `nt_zapline_plus` calls `nt_smooth` with three arguments, so `nodelayflag` is 0 and upstream's causal delay is *not* compensated; your change effectively removes that delay. That reads as a deliberate improvement, so it is carried forward as-is.
2. **`winSizeCompleteSpectrum` is dead (C3).** Once C4 replaced every `pwelch` call, nothing reads it, and the commented-out "data set is short" block that adjusted it is likewise inert. Kept commented rather than deleted, to keep the rebase diff minimal.

## Deferred

Not in this plan, by your scoping. The pwelch window is still `hanning(min segment length in the chunk)`, so a chunk edge landing a few samples before a boundary produces a very short window and the whole chunk's spectrum degrades. With `specNfft` in place this no longer crashes — measured behaviour with a 2-sample window is a spectrum flat to 1 part in 2500 across 55–65 Hz, whose argmax pins to the lowest bin, which drives `foundNoise` to 0 and silently falls that chunk back to fixed, non-adaptive removal at the global `noisefreq`. Assertion 6 of `test_pwelch_segmented` characterises this so it is visible rather than forgotten. The fix is to choose the window from a target length, keep only segments at least that long, and weight the average by the number of Welch windows per segment.
