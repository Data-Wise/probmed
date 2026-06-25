# P2 Submission Deliverables — Wave Plan

**Date:** 2026-06-25
**Scope:** Sequence and de-risk the three open deliverables tracked by pmed-modern
issues [#4](https://github.com/Data-Wise/pmed-modern/issues/4),
[#5](https://github.com/Data-Wise/pmed-modern/issues/5),
[#6](https://github.com/Data-Wise/pmed-modern/issues/6) for the Paper 2
(incremental mediated elasticity) JASA submission.
**Status:** approved structure; implementation pending.

---

## Context

Two P2 checklist items closed via probmed PR #10 (g-score EIF term II + `incr_pmed()`
promotion). Three remain open in `pmed-modern/02-incremental-pmed/SUBMISSION.md`:

| Line | Item | Issue | State |
|------|------|-------|-------|
| 20 | M–Y unmeasured-confounding sensitivity analysis | #4 | open, tractable |
| 21 | `incr_pmed()` vignette | #5 | open, documentation |
| 19 | Continuous / multi-valued exposure extension | #6 (tracker) | **unstarted**, theory-heavy |

**Key grounding finding.** `probmed::pmed_sensitivity()` (`R/sensitivity.R`) already
accepts `indirect`/`total` components directly and its docstring certifies the
additive-offset bias model is *exact for the Paper 2 incremental share* (numerator is a
difference of cross-world means). So #4 is **wiring + reporting**, not new theory — this
is why it leads.

**Cross-repo constraint.** This session is rooted in `probmed`. All `pmed-modern` writes
(manuscript prose, figures, SUBMISSION.md checkbox flips) are **deferred and
permission-gated** per the `no-cross-repo-work-without-permission` rule. Each wave below
splits its probmed (code) part from its pmed-modern (manuscript) part explicitly.

---

## Wave 1 — Sensitivity (unblocks the submission claim)

**Issue #4. Highest leverage, lowest risk.**

### probmed (code — this repo, actionable now)
- Compute, for the incremental curve, the per-δ `(med, tot)` components already produced
  by `incr_pmed()` (`R/incremental-pmed.R`: `med`, `tot` columns of `@curve`).
- Feed them into `pmed_sensitivity(indirect = med, total = tot, threshold = 0)` to get
  the **tipping bias** that zeroes the incremental share, per δ (or at a representative δ
  such as δ = 1).
- Decide surface: either (a) a thin convenience wrapper / method that maps an
  `IncrPmedResult` → a per-δ sensitivity table, or (b) document the direct-component call
  pattern. Prefer (a) only if it stays small; otherwise (b) avoids API bloat.
- Add a focused testthat case: components round-trip; tipping bias `= -med`; threshold
  bias `= thr*tot - med`; finite/positive where expected.

### pmed-modern (manuscript — deferred, permission-gated)
- Sensitivity row/figure for the JOBS II application across the δ-sweep.
- Prose in `manuscript/incremental-pmed.qmd`; re-render `.tex`/`.pdf`.
- Flip SUBMISSION.md line 20.

**Done when:** probmed exposes the incremental sensitivity path + test (suite green);
manuscript sensitivity result drafted (gated).

---

## Wave 2 — Vignette (documentation; reuses Wave 1)

**Issue #5. Depends on Wave 1 only for the sensitivity snippet.**

### probmed (code — this repo)
- New `vignettes/incremental-pmed.qmd` covering:
  - the δ-sweep and how to read the curve (flat ⇒ no A×M interaction; bends otherwise),
  - the g-cross-fit one-step estimator and pointwise CIs (full EIF incl. term II),
  - a short sensitivity callout reusing the Wave 1 path.
- Follow existing vignette conventions (`vignettes/gauge-residual.qmd` as template;
  native pipe, `package::fn()` namespacing).
- Ensure it builds under the project's quarto-vignette toolchain (watch the Windows
  quarto-vignette CI note in memory `ci-triggers-and-release-to-main`).

### pmed-modern (manuscript — deferred, permission-gated)
- Reference the vignette as the software illustration; flip the SUBMISSION.md line 21
  remainder.

**Done when:** vignette drafted and builds locally; referenced from P2 (gated).

---

## Wave 3 — Scope decision (decision gate, not implementation)

**Issue #6 tracker. The live work is the continuous / multi-valued exposure extension
(SUBMISSION.md line 19).**

- This wave is a **decision**, not code: full extension (modified treatment policies for
  continuous/multi-valued exposure) **vs.** a scope-limiting statement in the paper that
  defers it.
- Resolve explicitly before the two newly-checked boxes create the impression the
  checklist is complete (#6's stated purpose).
- If "full extension" is chosen, it becomes its own brainstorm → spec → plan cycle
  (new estimand + EIF, not a wiring task). If "scope statement", it is a short
  manuscript-prose edit (gated, pmed-modern).
- #6 stays open as the tracker either way.

**Done when:** scope call made and recorded; #6 updated accordingly.

---

## Dependencies & sequencing

```
Wave 1 (sensitivity, probmed)  ──┐
                                 ├──► Wave 2 (vignette, reuses sensitivity snippet)
Wave 3 (scope decision) ── independent, gated on user's scope call
```

- 1 → 2: vignette's sensitivity callout reuses Wave 1; core curve demo is independent, so
  Wave 2 can start in parallel and only the callout waits.
- 3 is independent and may be resolved at any time; it gates no other wave.

## Constraints (apply to every wave)

- **No pmed-modern writes without explicit per-instance permission** — manuscript,
  figures, and SUBMISSION.md checkbox flips are reported and wait for a "go".
- probmed code work follows the multi-branch workflow (feature branch → dev via PR).
- New vignette/test files are fine on a `feature/*` branch, not on `dev`.
