# Vendored NoiseTools (partial)

Source: `~/TMS_response_variability/NoiseTools` (the tree used by the TMS_response_variability
repo), copied 2026-09-14. Upstream project: http://audition.ens.fr/adc/NoiseTools
(public mirror: https://github.com/deron2/NoiseTools).

## This is a PARTIAL mirror, unlike vendor/zapline-plus

`vendor/zapline-plus` is a complete byte-for-byte mirror because RELAX shadows an
installed EEGLAB plugin with it, and a partial shadow would mix two versions. Nothing
shadows NoiseTools - it is not an EEGLAB plugin and is not otherwise installed - so
only the runtime closure of `nt_detrend` is vendored: 10 files instead of 138, and
6.2 MB of demo data and unrelated toolboxes left behind.

## Closure

`RELAX_repair_steps.m` calls exactly one NoiseTools entry point, `nt_detrend`. Its
runtime closure is:

    nt_detrend -> nt_regw, nt_demean
    nt_regw    -> nt_demean, nt_fold, nt_unfold
    nt_demean  -> nt_demean2, nt_fold, nt_unfold, nt_vecadd, nt_vecmult
    nt_fold    -> nt_greetings, nt_unfold
    nt_unfold  -> nt_fold, nt_greetings
    nt_vecmult -> nt_fold, nt_unfold
    nt_greetings -> nt_version

`nt_version` is not in the plan's original 9-file list. It was discovered during the
Step 2 closure check below: `nt_greetings.m` calls it bare (`nt_version` with no
parens, on the first invocation only - the result is cached in a persistent
variable) to build a banner string. It is a leaf (calls nothing else) and was copied
in, making this a 10-file closure.

Two names the naive grep-based check below also flags are NOT missing dependencies,
just false positives of matching `nt_[a-z0-9_]+` against non-call tokens:

  - `nt_detrend_helper` is a LOCAL SUBFUNCTION at `nt_detrend.m:163`, not a file -
    it can never appear as an entry in `present.txt` because `ls *.m` only sees
    top-level file names, so it will always show up in the `comm -23` diff.
  - `nt_greeted` is a `persistent` LOCAL VARIABLE declared at `nt_greetings.m:4`,
    not a function call at all - the regex can't distinguish a bare identifier
    used as a variable from one used as a niladic function call.

Three more names referenced from `nt_detrend.m` are genuinely unreachable at runtime
for how `RELAX_repair_steps` calls this function (always as `y = nt_detrend(x, ...)`,
capturing the output):

  - `nt_detrend2` and `nt_read_data` appear only inside the `if 0` demo block from
    `nt_detrend.m:269` onward.
  - `nt_imagescc` appears twice: once in that same `if 0` block, and once at
    `nt_detrend.m:157` inside `if ~nargout` (the plotting branch taken only when
    `nt_detrend` is called as a bare statement with no output captured). Any caller
    that does `y = nt_detrend(...)` never takes this branch.

To verify the closure after any change, run from this directory:

    grep -ohE '\bnt_[a-z0-9_]+' *.m | sort -u > /tmp/called.txt
    ls *.m | sed 's/\.m$//' | sort > /tmp/present.txt
    comm -23 /tmp/called.txt /tmp/present.txt

That must print exactly `nt_detrend2`, `nt_detrend_helper`, `nt_greeted`,
`nt_imagescc`, `nt_read_data` and nothing else. If any other name appears, that is a
genuine missing dependency - copy it from the source tree and re-run this check.

## Why the local tree rather than the public mirror

The public mirror's `nt_regw.m` has its all-zero-weight guard commented out, so a
degenerate fit silently returns `b = 0` and propagates `NaN` through the detrended
output. The local tree raises `error('weights all zero')` instead. A loud failure is
required here: `RELAX_repair_steps` catches it per sub-segment and steps down its
treatment ladder, which is only possible if the failure is signalled at all.

## Customizations

NONE. All ten files are byte-identical copies of the source tree. `RELAX_repair_steps`
works around `nt_detrend`'s rough edges from the outside (try/catch around each fit,
never handing it an all-zero-weight window) rather than patching it, so a later re-sync
is a straight overwrite.

## Re-sync procedure

1. Copy the files again from the source tree (or from upstream) - re-derive the file
   list from the closure check below rather than assuming it's still these ten.
2. Re-run the closure check above.
3. Run `tests/test_RELAX_repair_steps.m` - in particular the spectral-fidelity case,
   which fails if `nt_detrend`'s transfer function changed.
4. Update the checksum table below.

## Checksums at time of vendoring

| file | lines | sha256 |
| --- | --- | --- |
| `nt_demean2.m` | 57 | `18dfa324ab194ce2` |
| `nt_demean.m` | 73 | `43b6f3dccc90971e` |
| `nt_detrend.m` | 319 | `c962eabea1b8b5e3` |
| `nt_fold.m` | 28 | `e1b0bde503d51706` |
| `nt_greetings.m` | 16 | `abc03a2a77cb6f8c` |
| `nt_regw.m` | 120 | `818e19bfdcf85a5a` |
| `nt_unfold.m` | 21 | `b3ff85c24177cbc7` |
| `nt_vecadd.m` | 42 | `fa4604fda3d4aa1e` |
| `nt_vecmult.m` | 53 | `ad19b85d4b1e66bd` |
| `nt_version.m` | 2 | `fa8eedf8e153bbe4` |
