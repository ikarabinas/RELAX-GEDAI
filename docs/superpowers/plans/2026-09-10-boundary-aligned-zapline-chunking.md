# Boundary-Aligned Zapline Chunking Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the vendored zapline-plus adaptive chunker boundary-aligned, so no stretch of continuous data between two `boundary` events is ever split into a fragment shorter than the covariance/spectral machinery can use.

**Architecture:** Three new self-contained helper functions in `vendor/zapline-plus/` replace the flat, sample-1-anchored tiling inside `clean_data_with_zapline_plus.m`'s adaptive-chunk block (lines 411-476). The record is first cut into inter-boundary *segments*; the narrowband detection filter is applied per segment with odd reflection padding instead of across the whole array; covariance *epochs* are tiled inside each segment (one whole-segment epoch when the segment is short, ~1 s epochs with the last absorbing the remainder when it is long); and chunk edges are chosen from epoch starts, constrained so an edge either **is** a segment start or leaves at least `minChunkLength/2` seconds of that segment on both sides. `minChunkLength` is enforced in samples rather than epoch indices.

**Tech Stack:** MATLAB (Signal Processing Toolbox: `bandpass`, `filtord`, `findpeaks`, `pdist`), EEGLAB, vendored zapline-plus @ upstream `f4628a3`.

**Spec:** `docs/superpowers/specs/2026-09-08-boundary-aware-preprocessing-design.md` (parent boundary-handling design) plus the Background and Design sections of this document, which are the binding spec for this change.

---

## Background

`RELAX_Wrapper_slurm.m:326-332` calls zapline-plus with `chunkLength` unset, so zapline takes the **adaptive** chunking path (`clean_data_with_zapline_plus.m:411-476`). That block does three things, all anchored at sample 1 and all blind to `boundary` events:

1. **Line 413** narrowband-filters the whole array in one call:
   `data_narrowfilt = bandpass(data,[noisefreq-detectionWinsize/2 noisefreq+detectionWinsize/2],srate)`.
   The filter runs straight across every discontinuity, so the samples on each side of a boundary are contaminated by data from the other side.
2. **Lines 415-431** tile fixed `segmentLength`-second covariance epochs flatly from sample 1
   (`1+segmentLength*srate*(iSegment-1) : segmentLength*srate*iSegment`). Epochs straddle boundaries.
3. **Line 463** converts peak locations to chunk edges by the arithmetic `chunkIndices(2:end-1) = locs*segmentLength*srate`, which is only meaningful because every epoch is exactly `segmentLength*srate` samples long.

On the TMS-EEG inputs (`RELAX_cfg.RespectBoundaries=1`), m208_dlpfc_day1 has 60 boundaries and 61 inter-boundary segments of 7.930 s each. A chunk edge landing part-way through one of those segments leaves a fragment that can be a handful of samples long. That fragment then becomes the shortest segment in its chunk, and `pwelch_segmented` sizes its Hann window from the shortest segment — so the chunk's spectrum collapses to a near-flat estimate, `foundNoise` goes to 0, and the chunk silently falls back to fixed, non-adaptive removal. Before `docs/superpowers/plans/2026-09-09-zapline-nfft-grid-fix.md` pinned the FFT grid, the same fragment produced the hard crash *"Unable to perform assignment because the left and right sides have a different number of elements"* at line 492.

This plan removes the cause rather than tolerating the symptom.

### Measurements taken while designing this plan (do not re-derive)

- `bandpass(x,[57.5 62.5],500)` returns an **order-202 linear-phase FIR**, and the order is **independent of signal length** (identical for 3965-sample and 50000-sample input).
- `bandpass` is **single-pass with delay compensation, NOT zero-phase**: `max|bandpass(x) - singlepass_delay_compensated(x)| = 0` exactly, while `max|bandpass(x) - filtfilt(d,x)| = 0.45` against an output rms of 0.24. **Substituting `filtfilt` would be a large, silent behaviour change.** The per-segment filter must call `bandpass` itself.
- `filtfilt` on that filter needs `3*202+1 = 607` samples (1.21 s); one m208 segment is 3965 samples (7.93 s), so segments are comfortably long enough.
- With RELAX's config: `chunkLength = 0`, `segmentLength = 1`, `minChunkLength = 30`, `detectionWinsize = 5`. So the short-segment threshold `minChunkLength/2` is **15 s**, and every m208 segment (7.93 s) falls under it.

## Design

Let `E = round(segmentLength*srate)` (nominal epoch, 500 samples at 500 Hz), `W = max(round(minChunkLength/2*srate), E)` (the whole-segment threshold and the chunk-edge margin, 7500 samples at 500 Hz), and `M = max(E, filtord+1)` (the minimum usable epoch, 500 samples at 500 Hz).

**Segments.** Each element of `boundaries` is the **first sample of a new segment** (the wrapper writes `round(latency)`, and EEGLAB places boundaries at `x.5`, so `round` yields the first post-gap sample). Segments tile the record exactly.

**Covariance epochs.** For a segment of length `L`:
- `L < M` → emit **no** epoch. The segment is too short for a usable covariance and would otherwise contribute a rank-starved matrix that reads as a huge distance from its neighbours, manufacturing a spurious chunk edge.
- `M <= L < W` → emit **exactly one** epoch spanning the whole segment. This is the case that matters on this data.
- `L >= W` → emit `n = floor(L/E)` epochs tiled from the **segment's own start**; epochs `1..n-1` are exactly `E` samples, epoch `n` runs to the segment end. Every epoch is therefore in `[E, 2E)` and the last one ends exactly on the boundary.

**Chunk edges.** A peak at distance index `loc` proposes the edge `epochs(loc+1, 1)` — the **first sample of the following epoch**, not the last sample of the current one. (Upstream's `locs*segmentLength*srate` is the epoch's last sample, which under boundary alignment would place the edge one sample *before* a boundary and create a 1-sample fragment — precisely the bug being removed.) An edge is **legal** iff it is the enclosing segment's first sample, or it leaves at least `W` samples of that segment on both sides.

**minChunkLength.** `'MinPeakDistance', minChunkLength` is dropped from `findpeaks`, because it counts epoch indices and epochs are no longer uniformly 1 s (on m208 it would silently become 240 s and cut ~490 s of data into ~3 chunks). Instead the same rule is applied in the sample domain: candidates are visited in **descending peak height** — MATLAB's own `MinPeakDistance` rule — and an edge is accepted only if it is at least `minChunkLength*srate` samples from every already-accepted edge.

**Invariant this establishes:** *no chunk edge ever cuts an inter-boundary segment into a piece shorter than `minChunkLength/2` seconds; on data whose segments are all shorter than that threshold, every chunk is a whole number of complete segments.*

## Global Constraints

- `vendor/zapline-plus/` mirrors upstream `MariusKlug/zapline-plus` @ `f4628a3` byte-for-byte **except** where `vendor/zapline-plus/VENDORING.md` documents a local customization. Every change in this plan must be recorded there.
- New vendored helpers live **inside** `vendor/zapline-plus/`, following the precedent of `pwelch_segmented.m`.
- Vendored code may call `RELAX_reflect_pad` (the vendored tree is only ever reached through `RELAX_add_vendored_zapline()`, a RELAX function, so RELAX is on the path by construction). This dependency must be stated in `VENDORING.md`. Do **not** re-implement odd reflection padding.
- Do **not** transpose the full `data` array. It is `[samples x channels]` and can be 257 x 245000 doubles (~504 MB); helpers operate in that orientation and transpose only individual segments.
- Do **not** replace `bandpass` with `filtfilt` or with `RELAX_filtbutter`. See the measurement above.
- All line numbers in this plan are relative to commit `8a87fd7` ("Updated zapline-plus configuration and impedance based channel screening"). Verify the anchor text before editing by line number.
- The working tree carries uncommitted user work: modified `run_relax-gedai_slurm.sh`, and untracked `tests/test_RELAX_append_table_row.m` and `tests/test_RELAX_impedance_bad_channels.m` (which `run_all_tests` picks up by glob). **Never** `git add -A`, stash, revert, or clean. Stage only the specific files each task names.
- Every task ends with the full suite green:
  `matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"`
- Load MATLAB with `module load matlab` first. Do not run heavy work on `scu-login*` nodes; the test suite is light enough to run directly.
- Tests go in `tests/test_*.m` as parameterless functions; `run_all_tests` discovers them by glob. Every test must call `RELAX_add_vendored_zapline();` first.

## File Structure

**Create:**
- `vendor/zapline-plus/zapline_boundary_segments.m` — boundary vector → `[nSegments x 2]` inclusive `[start end]` rows tiling the record.
- `vendor/zapline-plus/zapline_covariance_epochs.m` — segments → boundary-respecting covariance epoch list plus each epoch's owning segment index.
- `vendor/zapline-plus/zapline_chunk_indices.m` — peak locations → legal, minimum-spaced `chunkIndices`.
- `vendor/zapline-plus/zapline_filter_by_segments.m` — apply a filter function per segment with odd reflection padding, in `[samples x channels]` orientation.
- `tests/test_zapline_chunk_segmentation.m`, `tests/test_zapline_chunk_indices.m`, `tests/test_zapline_filter_by_segments.m`, `tests/test_zapline_boundary_aligned_chunks.m`.

**Modify:**
- `vendor/zapline-plus/clean_data_with_zapline_plus.m:411-476` (the adaptive-chunk block) and the `res*`/`analyticsResults` accumulators near lines 280, 969-980 and 1065.
- `vendor/zapline-plus/VENDORING.md`, `CLAUDE.md`.

---

### Task 1: Segment and covariance-epoch geometry

**Files:**
- Create: `vendor/zapline-plus/zapline_boundary_segments.m`
- Create: `vendor/zapline-plus/zapline_covariance_epochs.m`
- Test: `tests/test_zapline_chunk_segmentation.m`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `segments = zapline_boundary_segments(nSamples, boundaries)` → `[nSegments x 2]` double, inclusive 1-based `[start end]`, contiguous, `segments(1,1)==1`, `segments(end,2)==nSamples`.
  - `[epochs, epochSegment] = zapline_covariance_epochs(segments, epochSamples, wholeSegmentSamples, minEpochSamples)` → `epochs` `[nEpochs x 2]` inclusive `[start end]`; `epochSegment` `[nEpochs x 1]` row index into `segments`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_zapline_chunk_segmentation.m`:

```matlab
function test_zapline_chunk_segmentation()
% Boundary-aligned covariance epochs for the vendored zapline-plus adaptive chunker.
% No epoch may straddle a boundary, and no inter-boundary segment short enough to
% be a single epoch may ever be subdivided.

    RELAX_add_vendored_zapline();

    srate  = 500;
    segLen = 3965;          % 7.93 s - the m208_dlpfc_day1 inter-train segment length
    nSeg   = 5;
    nSamp  = nSeg*segLen;
    bounds = (1:nSeg-1)'*segLen + 1;   % each boundary is the FIRST sample of a new segment

    % --- 1. segments tile the record exactly, one row per continuous stretch ---
    segments = zapline_boundary_segments(nSamp, bounds);
    assert(isequal(size(segments), [nSeg 2]), 'expected %d segments, got %d', nSeg, size(segments,1));
    assert(segments(1,1) == 1 && segments(end,2) == nSamp, 'segments must span the record');
    assert(all(segments(2:end,1) == segments(1:end-1,2) + 1), 'segments must be contiguous');
    assert(all(ismember(bounds, segments(:,1))), 'every boundary must start a segment');

    % --- 2. no boundaries at all degrades to a single whole-record segment ---
    assert(isequal(zapline_boundary_segments(nSamp, []), [1 nSamp]), ...
        'continuous data must produce exactly one segment');

    % --- 3. a segment below the whole-segment threshold gets exactly ONE epoch
    %        spanning it, and is never subdivided. 7.93 s < 15 s. ---
    epochSamples        = 1*srate;
    wholeSegmentSamples = 15*srate;
    minEpochSamples     = epochSamples;
    [epochs, epochSegment] = zapline_covariance_epochs(segments, ...
        epochSamples, wholeSegmentSamples, minEpochSamples);
    assert(size(epochs,1) == nSeg, 'each short segment must yield exactly one epoch (got %d)', size(epochs,1));
    assert(isequal(epochs, segments), 'a short segment''s epoch must span the whole segment');
    assert(isequal(epochSegment(:), (1:nSeg)'), 'epochSegment must map each epoch to its segment');

    % --- 4. no epoch may escape its segment, for a ragged layout ---
    raggedSamp = 100123;
    ragged = zapline_boundary_segments(raggedSamp, [3001; 3801; 40001; 40501]);
    [rEpochs, rSegIdx] = zapline_covariance_epochs(ragged, ...
        epochSamples, wholeSegmentSamples, minEpochSamples);
    for k = 1:size(rEpochs,1)
        seg = ragged(rSegIdx(k), :);
        assert(rEpochs(k,1) >= seg(1) && rEpochs(k,2) <= seg(2), ...
            'epoch %d [%d %d] escapes its segment [%d %d]', k, rEpochs(k,1), rEpochs(k,2), seg(1), seg(2));
    end

    % --- 5. a long segment is tiled from ITS OWN start, every epoch is at least
    %        segmentLength, and the last one ends exactly on the boundary ---
    longSeg = ragged(5,:);                       % 40501 .. 100123, 59623 samples
    inLong  = find(rSegIdx == 5);
    assert(rEpochs(inLong(1),1) == longSeg(1), 'tiling must start at the segment start');
    assert(rEpochs(inLong(end),2) == longSeg(2), 'the last epoch must end exactly on the boundary');
    lens = rEpochs(inLong,2) - rEpochs(inLong,1) + 1;
    assert(all(lens >= epochSamples), 'every epoch must be at least segmentLength (min was %d)', min(lens));
    assert(all(lens < 2*epochSamples), 'only the last epoch absorbs a remainder, so all stay under 2x');
    assert(lens(end) > epochSamples, 'this fixture was chosen to leave a real remainder');
    assert(isequal(rEpochs(inLong(2:end),1), rEpochs(inLong(1:end-1),2)+1), 'epochs must tile the segment');

    % --- 6. segments too short for a usable covariance emit no epoch at all,
    %        rather than contributing a rank-starved one that fakes a distance peak ---
    tiny = zapline_boundary_segments(nSamp, [201; 401]);   % two 200-sample slivers
    [tEpochs, tSegIdx] = zapline_covariance_epochs(tiny, ...
        epochSamples, wholeSegmentSamples, minEpochSamples);
    for k = 1:size(tEpochs,1)
        assert(tEpochs(k,2) - tEpochs(k,1) + 1 >= epochSamples, 'no epoch may be shorter than segmentLength');
    end
    assert(~any(ismember(tSegIdx, [1 2])), 'the two sub-500-sample segments must be dropped entirely');
end
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
module load matlab
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_zapline_chunk_segmentation"
```

Expected: FAIL with `Unrecognized function or variable 'zapline_boundary_segments'`.

- [ ] **Step 3: Create `vendor/zapline-plus/zapline_boundary_segments.m`**

```matlab
function segments = zapline_boundary_segments(nSamples, boundaries)
% ZAPLINE_BOUNDARY_SEGMENTS - stretches of continuous data between discontinuities.
%
%   nSamples   - number of samples in the record
%   boundaries - sample indices of discontinuities. Each value is the FIRST sample
%                of a new segment: the eeglab wrapper writes round(EEG.event.latency),
%                and EEGLAB places boundaries at x.5, so round() lands on the first
%                post-gap sample. Same convention as pwelch_segmented.
%
%   segments   - [nSegments x 2] inclusive 1-based [start end] rows. Contiguous and
%                exhaustive: segments(1,1) == 1, segments(end,2) == nSamples, and each
%                row starts one sample after the previous row ends. Empty or unusable
%                boundaries degrade to a single whole-record segment, so continuous
%                data behaves exactly as it did before.
%
% RELAX-GEDAI local addition to zapline-plus. See vendor/zapline-plus/VENDORING.md.

    if isempty(boundaries)
        segments = [1 nSamples];
        return;
    end

    b = unique(round(boundaries(:)'));
    b = b(b > 1 & b <= nSamples);

    if isempty(b)
        segments = [1 nSamples];
        return;
    end

    starts = [1, b];
    ends   = [b - 1, nSamples];
    segments = [starts(:), ends(:)];
end
```

- [ ] **Step 4: Create `vendor/zapline-plus/zapline_covariance_epochs.m`**

```matlab
function [epochs, epochSegment] = zapline_covariance_epochs(segments, epochSamples, wholeSegmentSamples, minEpochSamples)
% ZAPLINE_COVARIANCE_EPOCHS - epochs for adaptive chunk detection that respect boundaries.
%
% Upstream zapline tiles fixed-length covariance epochs flatly from sample 1, so an
% epoch can straddle a discontinuity and a chunk edge derived from it can slice a
% stretch of continuous data into a useless fragment. This tiles inside each segment
% instead.
%
%   segments            - [nSegments x 2] from zapline_boundary_segments
%   epochSamples        - nominal epoch length in samples (segmentLength*srate)
%   wholeSegmentSamples - segments shorter than this become ONE epoch spanning the
%                         whole segment, so they can never be subdivided
%   minEpochSamples     - segments shorter than this emit NO epoch. A covariance from
%                         too few samples is rank-starved and reads as a large distance
%                         from its neighbours, manufacturing a spurious chunk edge.
%
%   epochs              - [nEpochs x 2] inclusive [start end]. Every epoch lies wholly
%                         inside one segment and is at least epochSamples long.
%   epochSegment        - [nEpochs x 1] row index into segments for each epoch.
%
% RELAX-GEDAI local addition to zapline-plus. See vendor/zapline-plus/VENDORING.md.

    % Both thresholds are floored at epochSamples so "every epoch is at least
    % segmentLength long" holds however the caller parameterises this.
    minEpochSamples     = max(minEpochSamples, epochSamples);
    wholeSegmentSamples = max(wholeSegmentSamples, epochSamples);

    epochs       = zeros(0, 2);
    epochSegment = zeros(0, 1);

    for k = 1:size(segments, 1)
        segStart = segments(k, 1);
        segEnd   = segments(k, 2);
        segLen   = segEnd - segStart + 1;

        if segLen < minEpochSamples
            continue;
        end

        if segLen < wholeSegmentSamples
            epochs(end+1, :)       = [segStart segEnd];       %#ok<AGROW>
            epochSegment(end+1, 1) = k;                       %#ok<AGROW>
            continue;
        end

        nEpochs = floor(segLen / epochSamples);
        for j = 1:nEpochs
            thisStart = segStart + (j-1)*epochSamples;
            if j < nEpochs
                thisEnd = segStart + j*epochSamples - 1;
            else
                % The last epoch absorbs the remainder so it ends exactly on the
                % boundary and is never shorter than epochSamples.
                thisEnd = segEnd;
            end
            epochs(end+1, :)       = [thisStart thisEnd];     %#ok<AGROW>
            epochSegment(end+1, 1) = k;                       %#ok<AGROW>
        end
    end
end
```

- [ ] **Step 5: Run the test to verify it passes**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_zapline_chunk_segmentation"
```

Expected: PASS, no warnings.

- [ ] **Step 6: Run the full suite**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: every test file passes. Note `tests/test_RELAX_add_vendored_zapline.m` enumerates `vendor/zapline-plus/*.m` and asserts each resolves from the vendor directory — the two new files must satisfy it.

- [ ] **Step 7: Commit**

```bash
git add vendor/zapline-plus/zapline_boundary_segments.m \
        vendor/zapline-plus/zapline_covariance_epochs.m \
        tests/test_zapline_chunk_segmentation.m
git commit -m "Add boundary-respecting segment and covariance-epoch geometry for zapline chunking"
```

---

### Task 2: Legal, minimum-spaced chunk indices

**Files:**
- Create: `vendor/zapline-plus/zapline_chunk_indices.m`
- Test: `tests/test_zapline_chunk_indices.m`

**Interfaces:**
- Consumes: `zapline_boundary_segments`, `zapline_covariance_epochs` from Task 1, with exactly the signatures given there.
- Produces: `chunkIndices = zapline_chunk_indices(nSamples, epochs, epochSegment, segments, peakLocs, peakHeights, minChunkSamples, edgeMarginSamples)` → column vector, `chunkIndices(1)==1`, `chunkIndices(end)==nSamples+1`, strictly increasing. Consumed by `clean_data_with_zapline_plus.m` as `data(chunkIndices(i):chunkIndices(i+1)-1,:)`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_zapline_chunk_indices.m`:

```matlab
function test_zapline_chunk_indices()
% Chunk edges for the vendored zapline-plus adaptive chunker must never cut a stretch
% of continuous data into a fragment shorter than the edge margin, and minChunkLength
% must be enforced in samples rather than epoch indices.

    RELAX_add_vendored_zapline();

    srate  = 500;
    segLen = 3965;          % 7.93 s
    nSeg   = 40;
    nSamp  = nSeg*segLen;
    bounds = (1:nSeg-1)'*segLen + 1;

    epochSamples        = 1*srate;
    wholeSegmentSamples = 15*srate;
    minChunkSamples     = 30*srate;

    segments = zapline_boundary_segments(nSamp, bounds);
    [epochs, epochSegment] = zapline_covariance_epochs(segments, ...
        epochSamples, wholeSegmentSamples, epochSamples);

    % Propose an edge after every epoch, so nothing but the legality and spacing
    % rules can be what constrains the result.
    locs = (1:size(epochs,1)-1)';
    pks  = ones(size(locs));
    chunkIndices = zapline_chunk_indices(nSamp, epochs, epochSegment, segments, ...
        locs, pks, minChunkSamples, wholeSegmentSamples);

    % --- 1. shape ---
    assert(chunkIndices(1) == 1 && chunkIndices(end) == nSamp+1, ...
        'chunk indices must span the record');
    assert(all(diff(chunkIndices) > 0), 'chunk indices must be strictly increasing');

    % --- 2. On 7.93 s segments every epoch IS a whole segment, so every chunk edge
    %        must land exactly on a boundary. This is the invariant. ---
    interior = chunkIndices(2:end-1);
    assert(all(ismember(interior, bounds)), ...
        'chunk edges cut inside continuous segments: %s', mat2str(setdiff(interior, bounds)'));

    % --- 3. minChunkLength is enforced in SAMPLES: every chunk is at least 30 s even
    %        though epochs are 7.93 s each. ---
    chunkLens = diff(chunkIndices);
    assert(all(chunkLens >= minChunkSamples), ...
        'shortest chunk was %d samples, below the %d-sample minimum', min(chunkLens), minChunkSamples);

    % --- 4. ...and it must not degenerate. Interpreting minChunkLength as 30 EPOCHS
    %        would be 238 s here and would permit at most 2 chunks. ---
    assert(numel(chunkIndices)-1 >= 5, ...
        'expected many chunks from sample-domain spacing, got %d', numel(chunkIndices)-1);

    % --- 5. Inside a long segment an edge is legal only away from both its edges. ---
    longSamp = 100000;
    longSegs = zapline_boundary_segments(longSamp, 5001);
    [lEpochs, lSegIdx] = zapline_covariance_epochs(longSegs, ...
        epochSamples, wholeSegmentSamples, epochSamples);
    lLocs   = (1:size(lEpochs,1)-1)';
    lChunks = zapline_chunk_indices(longSamp, lEpochs, lSegIdx, longSegs, ...
        lLocs, ones(size(lLocs)), minChunkSamples, wholeSegmentSamples);
    for e = lChunks(2:end-1)'
        seg = longSegs(find(longSegs(:,1) <= e & longSegs(:,2) >= e, 1), :);
        assert(e == seg(1) || (e - seg(1) >= wholeSegmentSamples && seg(2) - e + 1 >= wholeSegmentSamples), ...
            'edge %d cuts segment [%d %d] closer than %d samples to an edge', ...
            e, seg(1), seg(2), wholeSegmentSamples);
    end
    assert(any(lChunks(2:end-1) > 5001), 'a long segment must still be splittable in its interior');

    % --- 6. Taller peaks win ties the way findpeaks'' MinPeakDistance does. ---
    % Both candidates sit far enough into the record that the first-chunk guard
    % cannot delete the survivor, so this isolates the spacing rule.
    twoLocs = [5; 6];
    twoPks  = [1; 5];                       % the second candidate is much taller
    twoChunks = zapline_chunk_indices(nSamp, epochs, epochSegment, segments, ...
        twoLocs, twoPks, minChunkSamples, wholeSegmentSamples);
    assert(numel(twoChunks) == 3, 'the two candidates are 3965 samples apart, so only one survives');
    assert(twoChunks(2) == epochs(7,1), 'the taller peak''s edge must be the survivor');

    % --- 7. No proposed peaks at all gives a single whole-record chunk. ---
    none = zapline_chunk_indices(nSamp, epochs, epochSegment, segments, [], [], ...
        minChunkSamples, wholeSegmentSamples);
    assert(isequal(none, [1; nSamp+1]), 'with no peaks the record must be one chunk');
end
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_zapline_chunk_indices"
```

Expected: FAIL with `Unrecognized function or variable 'zapline_chunk_indices'`.

- [ ] **Step 3: Create `vendor/zapline-plus/zapline_chunk_indices.m`**

```matlab
function chunkIndices = zapline_chunk_indices(nSamples, epochs, epochSegment, segments, peakLocs, peakHeights, minChunkSamples, edgeMarginSamples)
% ZAPLINE_CHUNK_INDICES - chunk edges that never fragment a stretch of continuous data.
%
%   nSamples          - number of samples in the record
%   epochs, epochSegment - from zapline_covariance_epochs
%   segments          - from zapline_boundary_segments
%   peakLocs          - findpeaks locations in the DISTANCE vector: peakLocs(i)==k means
%                       the covariance jumped between epoch k and epoch k+1
%   peakHeights       - the corresponding peak values, used to resolve spacing conflicts
%   minChunkSamples   - minimum chunk length in samples (minChunkLength*srate)
%   edgeMarginSamples - an edge inside a segment must leave at least this many samples
%                       of that segment on both sides
%
%   chunkIndices      - column vector, [1; edges; nSamples+1], consumed as
%                       data(chunkIndices(i):chunkIndices(i+1)-1,:)
%
% RELAX-GEDAI local addition to zapline-plus. See vendor/zapline-plus/VENDORING.md.

    candidates = zeros(0, 1);
    heights    = zeros(0, 1);

    for i = 1:numel(peakLocs)
        nextEpoch = peakLocs(i) + 1;
        if nextEpoch > size(epochs, 1)
            continue;
        end

        % The edge is the FIRST sample of the following epoch. Upstream used the last
        % sample of the current one, which under boundary alignment lands one sample
        % before a boundary and leaves a 1-sample fragment behind.
        edge = epochs(nextEpoch, 1);
        seg  = segments(epochSegment(nextEpoch), :);

        isSegmentStart = (edge == seg(1));
        leavesRoom     = (edge - seg(1) >= edgeMarginSamples) && ...
                         (seg(2) - edge + 1 >= edgeMarginSamples);

        if isSegmentStart || leavesRoom
            candidates(end+1, 1) = edge;             %#ok<AGROW>
            heights(end+1, 1)    = peakHeights(i);   %#ok<AGROW>
        end
    end

    % Reproduce findpeaks' MinPeakDistance rule in the SAMPLE domain: keep the tallest
    % peak, then drop every candidate within minChunkSamples of one already kept.
    % findpeaks itself cannot do this any more, because its distance is measured in
    % epoch indices and boundary-aligned epochs are not all the same length.
    [~, order] = sort(heights, 'descend');
    accepted = zeros(0, 1);
    for i = reshape(order, 1, [])
        if isempty(accepted) || all(abs(accepted - candidates(i)) >= minChunkSamples)
            accepted(end+1, 1) = candidates(i);      %#ok<AGROW>
        end
    end
    accepted = unique(accepted);

    chunkIndices = [1; accepted(:); nSamples + 1];

    % Upstream's guards: the first and last chunks must also reach the minimum length.
    if numel(chunkIndices) >= 3 && chunkIndices(2) - chunkIndices(1) < minChunkSamples
        chunkIndices(2) = [];
    end
    if numel(chunkIndices) >= 3 && chunkIndices(end) - chunkIndices(end-1) < minChunkSamples
        chunkIndices(end-1) = [];
    end
end
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_zapline_chunk_indices"
```

Expected: PASS, no warnings.

- [ ] **Step 5: Run the full suite**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

- [ ] **Step 6: Commit**

```bash
git add vendor/zapline-plus/zapline_chunk_indices.m tests/test_zapline_chunk_indices.m
git commit -m "Constrain zapline chunk edges to boundaries or well inside long segments"
```

---

### Task 3: Per-segment narrowband filtering

**Files:**
- Create: `vendor/zapline-plus/zapline_filter_by_segments.m`
- Test: `tests/test_zapline_filter_by_segments.m`

**Interfaces:**
- Consumes: `zapline_boundary_segments` from Task 1; `RELAX_reflect_pad(x, p)` from the repo root, which takes and returns `[nChannels x nSamples]`.
- Produces: `out = zapline_filter_by_segments(data, segments, fn, padWidth, minSegmentSamples)` where `data` and `out` are `[nSamples x nChannels]` and `fn` takes and returns `[nSamples x nChannels]`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_zapline_filter_by_segments.m`:

```matlab
function test_zapline_filter_by_segments()
% The narrowband filter that drives zapline's adaptive chunk detection must not smear
% data across a discontinuity.

    RELAX_add_vendored_zapline();

    srate  = 500;
    segLen = 3965;
    nSamp  = 3*segLen;
    bounds = [segLen+1; 2*segLen+1];
    band   = [57.5 62.5];

    rng(0);
    data     = randn(nSamp, 4);
    segments = zapline_boundary_segments(nSamp, bounds);

    % bandpass() designs an order-202 linear-phase FIR here, independent of signal
    % length, and applies it single-pass with delay compensation. Reuse it as-is.
    [~, d] = bandpass(data(:,1), band, srate);
    p  = filtord(d);
    fn = @(x) bandpass(x, band, srate);

    out = zapline_filter_by_segments(data, segments, fn, p, p+1);
    assert(isequal(size(out), size(data)), 'output must keep the input size');

    % --- 1. Changing one segment must not change any other segment's output. ---
    perturbed = data;
    perturbed(1:segLen, :) = randn(segLen, 4);
    outP = zapline_filter_by_segments(perturbed, segments, fn, p, p+1);
    tail = segLen+1:nSamp;
    assert(max(abs(out(tail,:) - outP(tail,:)), [], 'all') < 1e-12, ...
        'segment 1 leaked across the boundary into segments 2-3');

    % --- 2. Each segment must equal filtering that segment on its own. ---
    for k = 1:size(segments,1)
        s = segments(k,1);
        e = segments(k,2);
        alone = zapline_filter_by_segments(data(s:e,:), [1 e-s+1], fn, p, p+1);
        assert(max(abs(out(s:e,:) - alone), [], 'all') < 1e-12, ...
            'segment %d differs from filtering it in isolation', k);
    end

    % --- 3. Regression witness: the naive whole-record call DOES leak, which is the
    %        behaviour this helper exists to replace. ---
    naive  = bandpass(data,      band, srate);
    naiveP = bandpass(perturbed, band, srate);
    leak = max(abs(naive(tail,:) - naiveP(tail,:)), [], 'all');
    assert(leak > 1e-6, 'expected the naive whole-record filter to leak (leak was %.3g)', leak);

    % --- 4. Continuous data must match a plain whole-record call away from the record
    %        edges, where the reflection padding has no influence (FIR support is p). ---
    cont = zapline_filter_by_segments(data, [1 nSamp], fn, p, p+1);
    core = 2*p+1 : nSamp-2*p;
    assert(max(abs(cont(core,:) - naive(core,:)), [], 'all') < 1e-10, ...
        'away from the record edges, one-segment filtering must match a plain call');

    % --- 5. A segment too short to filter is passed through unchanged, with a warning
    %        rather than a crash. Such segments emit no covariance epoch either, so the
    %        unfiltered data never reaches a covariance. ---
    shortSegs = zapline_boundary_segments(nSamp, [101; 201]);
    lastwarn('');
    % evalc swallows the expected warning text so the suite's output stays pristine,
    % while lastwarn still records that the warning fired.
    [~] = evalc('outShort = zapline_filter_by_segments(data, shortSegs, fn, p, p+1);');
    assert(isequal(outShort(1:200,:), data(1:200,:)), ...
        'segments below the filter minimum must be copied through unchanged');
    assert(~isempty(lastwarn), 'skipping a segment must warn');
end
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_zapline_filter_by_segments"
```

Expected: FAIL with `Unrecognized function or variable 'zapline_filter_by_segments'`.

- [ ] **Step 3: Create `vendor/zapline-plus/zapline_filter_by_segments.m`**

```matlab
function out = zapline_filter_by_segments(data, segments, fn, padWidth, minSegmentSamples)
% ZAPLINE_FILTER_BY_SEGMENTS - apply a filter within each segment, never across one.
%
%   data              - [nSamples x nChannels], zapline's own orientation. The full
%                       array is never transposed: it can be 257 x 245000 doubles.
%   segments          - [nSegments x 2] from zapline_boundary_segments
%   fn                - function handle taking and returning [nSamples x nChannels]
%   padWidth          - odd-reflection pad in samples applied to each end of a segment,
%                       clamped to segLen-1. One filter order is enough: an FIR of
%                       order N cannot see further than N samples.
%   minSegmentSamples - segments shorter than this are copied through UNCHANGED with a
%                       warning. Callers must ensure such segments are excluded
%                       downstream; zapline_covariance_epochs drops them.
%
%   out               - [nSamples x nChannels], same size as data.
%
% Odd reflection is delegated to RELAX_reflect_pad, which works in [nChannels x nSamples],
% so only one segment at a time is transposed. The vendored tree is only ever reached
% through RELAX_add_vendored_zapline(), so that function is always on the path.
%
% RELAX-GEDAI local addition to zapline-plus. See vendor/zapline-plus/VENDORING.md.

    out = data;

    for k = 1:size(segments, 1)
        segStart = segments(k, 1);
        segEnd   = segments(k, 2);
        segLen   = segEnd - segStart + 1;

        if segLen < minSegmentSamples
            warning('zapline_filter_by_segments:shortSegment', ...
                'Segment %d (%d samples) is shorter than the filter needs; left unfiltered.', ...
                k, segLen);
            continue;
        end

        p        = min(padWidth, segLen - 1);
        padded   = RELAX_reflect_pad(data(segStart:segEnd, :).', p).';
        filtered = fn(padded);
        out(segStart:segEnd, :) = filtered(p+1 : p+segLen, :);
    end
end
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_zapline_filter_by_segments"
```

Expected: PASS with no output. Assertion 5 provokes a `zapline_filter_by_segments:shortSegment` warning on purpose and captures it with `evalc`, so nothing is printed.

- [ ] **Step 5: Run the full suite**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

- [ ] **Step 6: Commit**

```bash
git add vendor/zapline-plus/zapline_filter_by_segments.m tests/test_zapline_filter_by_segments.m
git commit -m "Apply zapline's narrowband detection filter per inter-boundary segment"
```

---

### Task 4: Wire the boundary-aligned chunker into clean_data_with_zapline_plus

**Files:**
- Modify: `vendor/zapline-plus/clean_data_with_zapline_plus.m` — the adaptive-chunk block at lines 411-476, plus the `res*` accumulator declaration near line 280, the per-noisefreq accumulator block near lines 969-980, and the `analyticsResults` assembly near line 1065
- Modify: `vendor/zapline-plus/VENDORING.md`
- Modify: `CLAUDE.md`
- Test: `tests/test_zapline_boundary_aligned_chunks.m`

**Interfaces:**
- Consumes: `zapline_boundary_segments`, `zapline_covariance_epochs` (Task 1), `zapline_chunk_indices` (Task 2), `zapline_filter_by_segments` (Task 3), with exactly the signatures given in those tasks.
- Produces: `analyticsResults.chunkIndices` — a `1 x nNoisefreqs` cell array whose i-th entry is the `chunkIndices` column vector used for noise frequency i. This is what makes the invariant observable to tests and to anyone auditing a run.

- [ ] **Step 1: Write the failing end-to-end test**

Create `tests/test_zapline_boundary_aligned_chunks.m`:

```matlab
function test_zapline_boundary_aligned_chunks()
% End-to-end: on data whose continuous stretches are all shorter than the
% whole-segment threshold, every zapline chunk must be a whole number of complete
% inter-boundary segments. Fixture mirrors m208_dlpfc_day1's 7.93 s inter-train
% segments at half the sampling rate.

    RELAX_add_vendored_zapline();

    srate  = 250;
    segLen = 1983;                 % 7.93 s
    nSeg   = 30;
    nChan  = 8;
    nSamp  = nSeg*segLen;
    bounds = (1:nSeg-1)'*segLen + 1;

    rng(0);
    t    = (0:nSamp-1)'/srate;
    data = randn(nSamp, nChan) + 8*sin(2*pi*60*t)*ones(1, nChan);

    cfg = struct('plotResults', 0, 'noisefreqs', 60, 'boundaries', bounds);
    [cleaned, ~, analytics] = clean_data_with_zapline_plus(data, srate, cfg);

    % --- 1. The chunk layout is now recorded, per noise frequency. ---
    assert(isfield(analytics, 'chunkIndices'), 'analyticsResults must record chunkIndices');
    assert(iscell(analytics.chunkIndices) && ~isempty(analytics.chunkIndices), ...
        'chunkIndices must be one cell per noise frequency');
    for iFreq = 1:numel(analytics.chunkIndices)
        ci = analytics.chunkIndices{iFreq};
        assert(ci(1) == 1 && ci(end) == nSamp+1, ...
            'noise frequency %d: chunk indices must span the record', iFreq);
        assert(all(diff(ci) > 0), 'noise frequency %d: chunk indices must increase', iFreq);
    end

    % --- 2. THE INVARIANT: no chunk edge falls inside a continuous segment. ---
    for iFreq = 1:numel(analytics.chunkIndices)
        ci = analytics.chunkIndices{iFreq};
        stray = setdiff(ci(2:end-1), bounds);
        assert(isempty(stray), 'noise frequency %d: chunk edges cut inside continuous segments: %s', ...
            iFreq, mat2str(stray'));
    end

    % --- 3. The cleaning still works: 60 Hz must be strongly attenuated. ---
    [pRaw,   f] = pwelch_segmented(data,    [], [], 8192, srate, bounds);
    [pClean, ~] = pwelch_segmented(cleaned, [], [], 8192, srate, bounds);
    noiseBin = f > 59.5 & f < 60.5;
    attenDb  = 10*log10(mean(pRaw(noiseBin,:), 'all') / mean(pClean(noiseBin,:), 'all'));
    assert(attenDb > 10, 'expected >10 dB attenuation at 60 Hz, got %.2f dB', attenDb);

    % --- 4. Every chunk's shortest continuous stretch is a FULL segment, so the
    %        pwelch window inside that chunk can never collapse onto a fragment. ---
    for iFreq = 1:numel(analytics.chunkIndices)
        ci = analytics.chunkIndices{iFreq};
        for k = 1:numel(ci)-1
            chunkStart = ci(k);
            chunkEnd   = ci(k+1) - 1;
            inChunk    = bounds(bounds > chunkStart & bounds <= chunkEnd);
            pieces     = diff([chunkStart; inChunk; chunkEnd+1]);
            assert(all(pieces == segLen), ...
                'noise frequency %d, chunk %d contains a partial segment (piece lengths %s)', ...
                iFreq, k, mat2str(pieces'));
        end
    end

    assert(all(isfinite(analytics.proportionRemovedNoise)), ...
        'proportionRemovedNoise must be finite');
end
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_zapline_boundary_aligned_chunks"
```

Expected: FAIL at assertion 1, `analyticsResults must record chunkIndices`. (If it instead fails at assertion 2, that is also a valid RED — it means the pre-change chunker is already cutting inside segments, which is the bug.)

- [ ] **Step 3: Replace the adaptive-chunk block**

In `vendor/zapline-plus/clean_data_with_zapline_plus.m`, replace lines **411-476 inclusive** — from `disp('Using adaptive chunk length!')` down to and including the `end` on line 476 that closes the `if chunkLength ~= 0 ... else` — with the block below, which supplies its own closing `end`. Leave the `if chunkLength ~= 0` fixed-chunk branch above untouched, and leave `nChunks = length(chunkIndices)-1;` (line 478) untouched.

```matlab
        disp('Using adaptive chunk length!')
        %% find chunk indices

        % RELAX-GEDAI local change (VENDORING.md C8): upstream anchors the narrowband
        % filter, the covariance epochs and the chunk edges at sample 1 and lets all
        % three run straight across discontinuities. On boundary-heavy recordings that
        % slices a stretch of continuous data into a fragment, which then becomes the
        % shortest segment in its chunk and collapses that chunk's pwelch window.
        segments = zapline_boundary_segments(size(data,1), boundaries);

        % bandpass() applies a linear-phase FIR single-pass and compensates the delay;
        % it is NOT zero-phase, so filtfilt is not a substitute. Design it once, then
        % apply the identical function inside each segment with reflection padding.
        narrowBand = [noisefreq-detectionWinsize/2 noisefreq+detectionWinsize/2];
        [~, narrowFilt] = bandpass(data(:,1), narrowBand, srate);
        narrowOrder = filtord(narrowFilt);
        data_narrowfilt = zapline_filter_by_segments(data, segments, ...
            @(x) bandpass(x, narrowBand, srate), narrowOrder, narrowOrder + 1);

        epochSamples        = round(segmentLength*srate);
        minEpochSamples     = max(epochSamples, narrowOrder + 1);
        wholeSegmentSamples = max(round(minChunkLength/2*srate), epochSamples);

        [epochs, epochSegment] = zapline_covariance_epochs(segments, ...
            epochSamples, wholeSegmentSamples, minEpochSamples);

        nEpochs = size(epochs,1);
        fprintf('%.0f covariance epochs across %.0f continuous segments.\n', ...
            nEpochs, size(segments,1));

        %% compute covmatrices
        covarianceMatrices = zeros(size(data_narrowfilt,2),size(data_narrowfilt,2),nEpochs);

        for iEpoch = 1:nEpochs

            covarianceMatrices(:,:,iEpoch) = ...
                cov(data_narrowfilt(epochs(iEpoch,1):epochs(iEpoch,2),:));

        end

        %% find distances
        distances = zeros(max(nEpochs-1,0),1);

        for iEpoch = 2:nEpochs

            distances(iEpoch-1) = ...
                sum(pdist(covarianceMatrices(:,:,iEpoch)-covarianceMatrices(:,:,iEpoch-1)))/2;

        end

        %% find peaks
        % 'MinPeakDistance', minChunkLength is deliberately NOT passed: it counts epoch
        % INDICES, and boundary-aligned epochs are no longer uniformly segmentLength
        % long, so it would silently mean minChunkLength * (mean epoch length) seconds.
        % zapline_chunk_indices enforces minChunkLength in samples instead.
        pks  = [];
        locs = [];
        if numel(distances) >= 3
            [~,~,~,proms] = findpeaks(distances);
            if ~isempty(proms)
                [pks,locs] = findpeaks(distances,'MinPeakProminence',quantile(proms,prominenceQuantile));
            end
        end

        %% create final chunk indices
        chunkIndices = zapline_chunk_indices(size(data,1), epochs, epochSegment, segments, ...
            locs, pks, round(minChunkLength*srate), wholeSegmentSamples);
    end
```

- [ ] **Step 4: Record the chunk layout in the analytics**

Three edits, all in `vendor/zapline-plus/clean_data_with_zapline_plus.m`.

(a) Immediately after the `[... resFoundNoise] = deal([]);` statement that ends near line 280, add:

```matlab
% RELAX-GEDAI local change (VENDORING.md C9): the chunk layout is recorded so a run's
% boundary alignment can be audited after the fact. Cell, because the number of chunks
% differs per noise frequency.
resChunkIndices = {};
```

(b) In the per-noisefreq accumulator block (with `resFoundNoise(i_noisefreq,...)` near line 972), add:

```matlab
    resChunkIndices{i_noisefreq} = chunkIndices;
```

(c) Beside `analyticsResults.NremoveFinal = resNremoveFinal;` near line 1065, add:

```matlab
analyticsResults.chunkIndices = resChunkIndices;
```

- [ ] **Step 5: Run the end-to-end test to verify it passes**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); test_zapline_boundary_aligned_chunks"
```

Expected: PASS.

- [ ] **Step 6: Run the full suite**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"
```

Expected: every test file passes. `tests/test_zapline_boundary_crash.m` (the nfft-grid regression) and `tests/test_pwelch_segmented.m` must both still pass — if either breaks, stop and report rather than editing them.

- [ ] **Step 7: Update `vendor/zapline-plus/VENDORING.md`**

Follow the existing C1-C7 table/section format exactly. Add:

- **C8 — `clean_data_with_zapline_plus.m`, adaptive-chunk block.** Replaced the sample-1-anchored narrowband filter, covariance tiling and chunk-edge arithmetic with boundary-aligned equivalents built on `zapline_boundary_segments`, `zapline_covariance_epochs`, `zapline_chunk_indices` and `zapline_filter_by_segments`. Chunk edges are now the first sample of the following epoch (upstream used the last sample of the current one). `'MinPeakDistance', minChunkLength` was removed from `findpeaks` because it counts epoch indices; `minChunkLength` is now enforced in samples inside `zapline_chunk_indices`.
- **C9 — `clean_data_with_zapline_plus.m`, analytics.** Added `analyticsResults.chunkIndices`, a cell of the chunk layout per noise frequency.

Also add to the "files with no upstream counterpart" inventory, alongside `pwelch_segmented.m`: `zapline_boundary_segments.m`, `zapline_covariance_epochs.m`, `zapline_chunk_indices.m`, `zapline_filter_by_segments.m`. State explicitly that `zapline_filter_by_segments.m` calls `RELAX_reflect_pad` from the repo root, so the vendored tree is not standalone.

- [ ] **Step 8: Update `CLAUDE.md`**

In the "Boundary/discontinuity handling" section, extend the sentence about the vendored zapline-plus to record that its **adaptive chunker** is now boundary-aligned: covariance epochs never straddle a boundary, an inter-boundary segment shorter than `minChunkLength/2` is never subdivided, chunk edges either fall on a boundary or leave at least `minChunkLength/2` on both sides, and the narrowband detection filter runs per segment. Note that this holds only for the vendored copy reached via `RELAX_add_vendored_zapline()`; the EEGLAB plugin copy remains sample-1-anchored.

- [ ] **Step 9: Run the full suite one final time and commit**

```bash
matlab -nodisplay -nosplash -batch "addpath(pwd); addpath(fullfile(pwd,'tests')); run_all_tests"

git add vendor/zapline-plus/clean_data_with_zapline_plus.m \
        vendor/zapline-plus/VENDORING.md \
        CLAUDE.md \
        tests/test_zapline_boundary_aligned_chunks.m
git commit -m "Align zapline's adaptive chunker to data boundaries"
```

**Do not** `git add -A`. The working tree holds unrelated uncommitted user edits.

---

## Manual verification on real data

Not a subagent step — this is a ~60 minute cluster job for the controller or the user to run and interpret after Task 4 lands.

```bash
sbatch --array=<m208-index> run_relax-gedai_slurm.sh
```

In the resulting log, check:

1. `N covariance epochs across M continuous segments.` — for m208_dlpfc_day1 expect M = 61 and N = 61 (every segment is 7.93 s, below the 15 s threshold, so exactly one epoch each).
2. `X chunks will be created.` — expect several chunks, each at least 30 s, not 1 and not 61.
3. No `zapline_filter_by_segments:shortSegment` warnings.
4. `Removing noise at <f>Hz` still resolves a frequency off the coarse grid.
5. Compare `proportion of removed power` against the pre-change run's `-0.12351`, and the `Nremove is larger than 1/5th of the components, using that (51)!` escalations. If the negative proportion and the 51-of-257 escalation persist, that confirms they are driven by the separate stopping-test grid issue listed under Out of scope, not by chunk fragmentation.

## Plan validation already performed

Every code block in this plan was extracted, assembled and executed before the plan was handed over, so the implementer should expect it to work as written rather than treat it as pseudocode:

- `test_zapline_chunk_segmentation`, `test_zapline_chunk_indices` and `test_zapline_filter_by_segments` were run against the Task 1-3 implementations exactly as printed here: **3/3 PASS**.
- The Task 4 wiring block was spliced into a scratch copy of `clean_data_with_zapline_plus.m` at lines 411-476, the three analytics edits were applied, and `test_zapline_boundary_aligned_chunks` was run end to end: **PASS**, with every chunk edge landing on a boundary and 60 Hz attenuated by well over the asserted 10 dB.

If a step here fails during implementation, the likely cause is a transcription difference, not a design error — re-read the block before changing it.

## Out of scope

Recorded so they are not silently absorbed into this change:

- **Whole-data spectrum grid.** The three whole-data `pwelch_segmented` calls (lines 621/623/625) share `specNfft` with the per-chunk calls, leaving roughly 3 grid bins inside `detailedFreqBoundsUpper` against a 0.5% threshold, which quantises the "cleaning too weak" test to 0/33/67/100%. Upstream had ~52 bins there via a 300 s window. Deliberately deferred to a separate plan. **Confirmed independent of this change** during plan validation: with the boundary-aligned chunker in place and every chunk edge on a boundary, the fixture still logged `100% of frequency samples above thresh` and escalated `Cleaning too weak!` twice — while simultaneously reporting `proportion of removed power at noise frequency: 0.99928`. The stopping test fires on data it has already cleaned almost perfectly, so chunk fragmentation is not what drives it.
- **`pwelch_segmented`'s window is still the minimum segment length.** This plan makes that minimum a *complete* segment rather than a fragment, which is what mattered; decoupling the window from the minimum entirely remains open.
- **`VENDORING.md`'s C7 entry is stale.** Commit `8a87fd7` reverted `nt_smooth.m` from `filtfilt` back to `filter` and added `nodelayflag=1` at both `nt_zapline_plus.m` call sites, but C7 still documents the change in the old direction and never mentions `nodelayflag`. Reconciling that is the user's call, not this plan's.
- **The `noisefreqs='line'` 49-50/59-60 vs 49-51/59-61 upstream bug** (U2 in VENDORING.md). Dead code in RELAX.
