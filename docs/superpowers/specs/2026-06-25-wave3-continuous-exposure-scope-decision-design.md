# Wave 3 — Continuous/Multi-Valued Exposure: Scope Decision

**Date:** 2026-06-25
**Topic:** Issue #6 (Data-Wise/pmed-modern) last open item — continuous / multi-valued
exposure extension for the incremental mediated elasticity `P_med^δ`.
**Decision owner:** D. Tofighi · **Recorded by:** Claude (brainstorming skill)

## Decision

**Hybrid: scope-limiting statement for the Paper-2 submission now; the full
continuous-exposure method becomes a separately-tracked Paper-3 line.**
Cross-repo work is dispatched via **GitHub issues** on `Data-Wise/pmed-modern`,
not direct edits (standing no-cross-repo-edit rule).

## Rationale (the grill)

The incremental propensity-score intervention is **intrinsically binary**. Its
engine is the odds-ratio tilt

```
q_δ(C) = δ g / (δ g + 1 − g),   g = P(A = 1 | C)
```

which requires `g` to be a Bernoulli probability. There is no `δ` that tilts a
continuous `A` this way. The correct object for continuous/multi-valued exposure
is a **modified treatment policy (MTP)** / shift intervention `A → A + δ`
(Díaz & van der Laan 2012/2018; Haneuse & Rotnitzky 2013) — a *different*
estimand, EIF, mediation decomposition, identification set, and simulation suite.

Therefore continuous-exposure incremental `P_med` is **Paper-3-magnitude
methodology, not a Wave-3 patch to Paper 2**. Paper 2's identity is binary-exposure
incremental elasticity; bolting on a half-built second estimand mid-submission
delays a ready paper and invites reviewer scope-creep.

Rejected alternatives:
- **"Method must be general"** — would gate a ready submission on a new paper's worth of theory.
- **Bare scope note** — an unexplained "binary only" reads as a gap; reviewers punish
  unjustified scope limits. A *justified* boundary (binary IPS is the natural incremental
  object; continuous needs MTPs, named as future work) reads as command of the literature.

## Deliverables

| # | Deliverable | Location | Channel | Status |
|---|-------------|----------|---------|--------|
| 1 | Scope paragraph (Discussion/Limitations) + SUBMISSION.md line-19 flip | `pmed-modern/02-incremental-pmed/` | **GitHub issue** | filed |
| 2 | Paper-3 stub spec — MTP-mediation estimand (shift intervention, EIF to derive, sim plan) | research repo | **GitHub issue** | filed |
| 3 | Roxygen "Scope (exposure type)" note — binary-by-construction | `probmed/R/incremental-pmed.R` | direct (this repo) | done (this branch) |

Note: issue #6's other two open items (M–Y sensitivity, vignette) were **already
closed** by Wave 1 (`incr_sensitivity`) and Wave 2 (`incremental-pmed.qmd`). The
scope statement is the single remaining item between current state and a complete
P2 checklist.

## Content for the scope paragraph (drop-in)

> The incremental intervention studied here is defined for a binary exposure: the
> odds-ratio tilt `q_δ(C) = δg/(δg+1−g)` requires the propensity `g` to be a
> Bernoulli probability. Continuous and multi-valued exposures fall outside this
> construction and call instead for modified treatment policies — stochastic shift
> interventions of the form `A → A + δ` — which define a distinct mediated-elasticity
> estimand with its own efficient influence function. We treat that extension as
> future work rather than a special case of the present estimator.

## Out of scope (explicitly)

- No continuous-exposure code in `probmed` for Paper 2.
- No MTP EIF derivation in this cycle (it is the Paper-3 deliverable, tracked by issue #2 above).
