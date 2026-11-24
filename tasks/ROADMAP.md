# probmed roadmap

**Reviewed 2026-11-23.** This roadmap records the package's current priorities.
[.STATUS](../.STATUS) owns the current state, and the [task checklist](todo.md)
links to the acceptance rules and evidence for each phase.

## Current position

- **Package.** `dev` integrates through `main` for releases. Version
  [0.4.1](https://github.com/Data-Wise/probmed/releases/tag/v0.4.1) is the
  latest GitHub release. Subsequent changes on `dev` are development work.
  medfit provides model extraction, and probmed owns its P_med computation
  and bootstrap. The package release and pmed submission are on hold.
- **Phase 7.** I5 is complete and the revised pmed Results wording was
  approved. I6 failed the prespecified coverage and width gates at s6,
  N = 1,000, so **I7 is stopped**. The separate pmed exact-Gaussian follow-up
  passed its own six-cell rule. [Issue #47](https://github.com/Data-Wise/probmed/issues/47)
  is the only open probmed issue as of this review and remains held with I7.
  See the [Phase 7 plan](plan-mc-bootstrap.md).
- **Phase 8 remains NO-GO.** Retain the chi-square interval with one degree
  of freedom and the `near_double_null` diagnostic. The offline T5 audit
  closed [#50](https://github.com/Data-Wise/probmed/issues/50) after
  [PR #70](https://github.com/Data-Wise/probmed/pull/70) (`dbb4292`).
  No candidate passed both gates. `MC21_t_05` was inconclusive, and
  `MC21_t_all` was excluded from interval evidence. See the
  [Phase 8 plan](plan-mbco-null.md).

## Recent closeouts

| Item | Disposition |
| --- | --- |
| [#38](https://github.com/Data-Wise/probmed/issues/38) | Closed after the Gauge and Sobol manuscript fixes and PDF verification in pmed-modern at `d29215f` and `dc8b25b`. |
| [PR #71](https://github.com/Data-Wise/probmed/pull/71) | Merged into `dev` as `bc16bfd`. The fixed-treatment wording and scope of `pmed_sensitivity()` were verified locally. See the [claims audit plan](PLAN-research-workflow-enhancements-2026-10-20.md). |

## Next work and holds

1. **Available now.** Reconcile the remaining [mediationverse compatibility
   documentation item](todo-mediationverse-compat.md)
   against the README ecosystem table and `.STATUS`, using the existing
   compatibility evidence.
2. **After an authorized site deployment.** Verify the corrected vignettes
   and sensitivity help on the public pkgdown site.
3. **Held.** Keep I7 and #47 stopped under the I6 decision. Keep Phase 8 at
   NO-GO and the package release and pmed submission on hold. The available
   documentation work uses existing evidence and requires no new simulations.
