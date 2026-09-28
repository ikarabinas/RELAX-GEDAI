# Parked tests

Not discovered by `run_all_tests` (it only globs `tests/test_*.m`). Parked 2026-09-24 by the user:

- `test_RELAX_detect_reference_pops.m`, `test_RELAX_bridged_channels.m` — the reference-pop and
  bridged-channel features are in limbo and still under development. Both failed at the time
  ("total amplitude -235.9"; "The cutoff frequencies must be within the interval of (0,1)", after
  an uncommitted 2026-09-18 edit switched `RELAX_bridged_channels`'s `Filter` default to 1).
  Move them back to `tests/` when work on those features resumes.

`detect_steps_plateau_hold_sections.m` (sections 2d/2e, the plateau- and hold-drift gates) was
deleted 2026-09-25 with the gates themselves (plan Task 4D). The same change removed the
`plateauDrift`/`holdDrift` assertions from `test_RELAX_detect_reference_pops.m`.
