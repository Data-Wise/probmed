# probmed roadmap

**Reviewed 2026-11-24.** This roadmap records the package's current priorities.
[.STATUS](../.STATUS) owns the current state, and the [task checklist](todo.md)
links to the acceptance rules and evidence for each phase.

## Current position

- **Package.** `dev` integrates through `main` for releases. Version
  [0.4.1](https://github.com/Data-Wise/probmed/releases/tag/v0.4.1) is the
  latest GitHub release. The 0.4.2 candidate (`d8d7e9c`) passed local
  preflight. All five hosted jobs failed before probmed tests could run
  because `medfit 0.3.2` failed to load with `S7 1.0.0`.
  [Hosted run](https://github.com/Data-Wise/probmed/actions/runs/45563290371)
  medfit provides model extraction, and probmed owns its P_med computation
  and bootstrap.
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
  `MC21_t_all` was excluded from interval evidence. The
  [Phase 8 research follow-up](plan-phase8-follow-up.md) has started.
  The first all-target Gaussian rule controls size, but reaching the
  ideal-origin target gives empty sets and wider strong-path intervals.

## Recent closeouts

| Item | Disposition |
| --- | --- |
| [#38](https://github.com/Data-Wise/probmed/issues/38) | Closed after the Gauge and Sobol manuscript fixes and PDF verification in pmed-modern at `d29215f` and `dc8b25b`. |
| [PR #71](https://github.com/Data-Wise/probmed/pull/71) | Merged into `dev` as `bc16bfd`. The vignette and sensitivity help corrections are live and verified from `gh-pages` commit `9d4ba50`. See the [claims audit plan](PLAN-research-workflow-enhancements-2026-10-20.md). |
| [mediationverse compatibility T6](todo-mediationverse-compat.md) | Closed using the saved test evidence. The README records the existing medfit floor and the missing-value guard's release. |

## Next work and holds

For the next package task, choose and document a compatible medfit source
for hosted CI and for GitHub installs. The local `medfit 0.20.0` and
`S7 1.0.0` install and plug-in fit passed. Verify the chosen dependency
setup with a clean install and the full package checks, then run the
hosted matrix on the updated candidate. See the
[compatibility plan](plan-mediationverse-compat.md#2026-11-24-release-preflight-follow-up).

I7 and #47 remain stopped under the I6 decision. Phase 8 remains NO-GO.
Package publishing and pmed submission remain on hold.
