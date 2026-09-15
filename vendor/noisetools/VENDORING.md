# Vendored NoiseTools (partial)

Source: `~/TMS_response_variability/NoiseTools` (the tree used by the TMS_response_variability
repo), copied 2026-09-14. Upstream project: http://audition.ens.fr/adc/NoiseTools
(public mirror: https://github.com/deron2/NoiseTools).

## This is a PARTIAL mirror, unlike vendor/zapline-plus

`vendor/zapline-plus` is a complete byte-for-byte mirror because RELAX shadows an
installed EEGLAB plugin with it, and a partial shadow would mix two versions. Nothing
shadows NoiseTools - it is not an EEGLAB plugin and is not otherwise installed - so
only the runtime closure of `nt_detrend` is vendored, plus one harmless orphan file
carried over from the original plan (`nt_demean2` - see "Closure" below): 10 files
instead of 138, and 6.2 MB of demo data and unrelated toolboxes left behind.

## Closure

`RELAX_repair_steps.m` calls exactly one NoiseTools entry point, `nt_detrend`. Its
runtime closure, re-derived by reading every call site directly (excluding commented-out
lines and the unreachable `if 0` demo blocks - see "Fix round 1" note below), is:

    nt_detrend   -> nt_greetings, nt_regw
    nt_regw      -> nt_unfold, nt_demean, nt_fold
    nt_demean    -> nt_unfold, nt_vecadd, nt_vecmult, nt_fold
    nt_fold      -> nt_unfold            (only when the input is already 3-D)
    nt_unfold    -> (nothing)
    nt_vecadd    -> nt_unfold, nt_fold
    nt_vecmult   -> nt_unfold, nt_fold
    nt_greetings -> nt_version

`nt_detrend` calls `nt_greetings` directly at `nt_detrend.m:20`, and again from its
local subfunction `nt_detrend_helper` at `nt_detrend.m:197`; `nt_detrend_helper` calls
`nt_regw` at `nt_detrend.m:245`. Both calls are folded into the `nt_detrend` node above
since `nt_detrend_helper` is not a separate file (see the grep-false-positives note
below).

`nt_fold.m:11` and `nt_unfold.m:10` each contain a commented-out `%nt_greetings;` line -
dead code from an earlier version, not a live dependency. Neither function calls
`nt_greetings`, or anything else, at runtime: `nt_unfold` calls no other NoiseTools
function at all, and `nt_fold`'s only call is the conditional `nt_unfold` above.

`nt_demean2.m` is vendored but is **not actually reached by this closure**. Direct
inspection of `nt_demean.m` shows it never references `nt_demean2` anywhere - the
dependency runs the other way: `nt_demean2` calls `nt_demean` (`nt_demean2.m:52,55`,
to demean each page of a 4-D array) and otherwise only recurses into itself
(`nt_demean2.m:29,32,34`). Nothing in this closure calls `nt_demean2`. It was included
in the original 9-file list on the mistaken assumption that `nt_demean` called it; that
assumption was wrong, but the file is left in place (this is a documentation-only fix -
see "Fix round 1" below) since it is harmless: an uncalled file changes nothing about
`nt_detrend`'s behavior, and removing it is a decision for a future vendoring pass, not
this one.

`nt_version` is not in the plan's original 9-file list either. It was discovered during
the Step 2 closure check below: `nt_greetings.m` calls it bare (`nt_version` with no
parens, on the first invocation only - the result is cached in a persistent variable) to
build a banner string. It is a leaf (calls nothing else) and was copied in, making this
a 10-file vendored set (9 in the true closure above, plus the harmless orphan
`nt_demean2`).

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
   list from the closure check below rather than assuming it's still these ten. Note
   that the check only flags *missing* files (a name called but not present); it says
   nothing about *unreached* ones like `nt_demean2`, so it won't catch a similar orphan
   on its own - cross-check against the hand-derived graph in "Closure" above too.
2. Re-run the closure check above.
3. Run `tests/test_RELAX_repair_steps.m` - in particular the spectral-fidelity case,
   which fails if `nt_detrend`'s transfer function changed.
4. Update the checksum table below.

## Fix round 1 (2026-09-14)

The initial version of the "Closure" section above was copied verbatim from the
implementation plan, which built it by grepping for `nt_*` tokens without excluding
commented-out lines or the `if 0` demo blocks. That produced four errors, caught in
review and corrected here: two false edges from commented-out `%nt_greetings;` calls
(`nt_fold -> nt_greetings`, `nt_unfold -> nt_fold, nt_greetings` - `nt_unfold` in fact
calls nothing), one imprecise edge (`nt_detrend -> nt_demean` was asserted directly, but
the only direct call site is in the unreachable `if 0` block - `nt_demean` is reached
only transitively, via `nt_regw`), and one omitted real edge (`nt_detrend -> nt_greetings`,
called bare at `nt_detrend.m:20` and `nt_detrend.m:197`). While re-deriving the graph to
fix those four, direct inspection also turned up a fifth, unreported error: the original
diagram's `nt_demean -> nt_demean2` edge is backwards and the true edge
(`nt_demean2 -> nt_demean`) is unreachable from `nt_detrend` - see the `nt_demean2`
paragraph above. The vendored file set itself was not changed; the `.m` files remain
byte-identical to the source tree.

## Checksums at time of vendoring

| file | lines | sha256 |
| --- | --- | --- |
| `nt_demean2.m` (vendored but unreached - see "Closure" above) | 57 | `18dfa324ab194ce2` |
| `nt_demean.m` | 73 | `43b6f3dccc90971e` |
| `nt_detrend.m` | 319 | `c962eabea1b8b5e3` |
| `nt_fold.m` | 28 | `e1b0bde503d51706` |
| `nt_greetings.m` | 16 | `abc03a2a77cb6f8c` |
| `nt_regw.m` | 120 | `818e19bfdcf85a5a` |
| `nt_unfold.m` | 21 | `b3ff85c24177cbc7` |
| `nt_vecadd.m` | 42 | `fa4604fda3d4aa1e` |
| `nt_vecmult.m` | 53 | `ad19b85d4b1e66bd` |
| `nt_version.m` | 2 | `fa8eedf8e153bbe4` |
