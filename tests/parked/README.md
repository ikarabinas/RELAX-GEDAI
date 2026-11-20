# Parked tests

Not discovered by `run_all_tests` (it only globs `tests/test_*.m`). Parked 2026-09-24 by the user:

- `test_RELAX_bridged_channels.m` — the bridged-channel feature is in limbo and still under
  development. It failed at the time ("The cutoff frequencies must be within the interval of
  (0,1)", after an uncommitted 2026-09-18 edit switched `RELAX_bridged_channels`'s `Filter`
  default to 1). Move it back to `tests/` when work on that feature resumes.
- `test_RELAX_detect_reference_pops.m` was parked alongside it and returned to `tests/` on
  2026-11-20, passing; reference-pop repair is on in the slurm entry script.

Parked 2026-11-20:

- `test_GEDAI_alpha_eyeoc.m` — written ahead of `GEDAI_alpha_eyeoc_summary` /
  `GEDAI_alpha_eyeoc`, which do not exist yet ("Undefined function 'GEDAI_alpha_eyeoc_summary'").
  Move it back to `tests/` once they do.

`detect_steps_plateau_hold_sections.m` (sections 2d/2e, the plateau- and hold-drift gates) was
deleted 2026-09-25 with the gates themselves (plan Task 4D). The same change removed the
`plateauDrift`/`holdDrift` assertions from `test_RELAX_detect_reference_pops.m`.
