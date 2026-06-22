# wasserstein_pmed() — integration spec (planning → R/ promotion)

Prototype for **pmed-modern Paper 4** ($P_{med}^{W}$, transport scale). Mirrors the gauge/sobol
graduation pattern: verified prototype lives here in `planning/pmed-extensions/`; promotion to
`probmed/R/` (S7 + NAMESPACE + tests + docs) is a **Claude Code** task (per Cowork CLAUDE.md, package
R/ coding belongs in Claude Code via craft/rforge — not Cowork).

## Branch
`feature/wasserstein-pmed` (worktree off `dev`). Do **not** open the PR until integrating; rebase on
`dev`; never commit to `main` (protected).

## What the prototype delivers (VERIFIED, d=1)
- `.w2_1d()` — 1-D quadratic Wasserstein via quantile transport (Brenier).
- `.potIF()` — Kantorovich-potential influence function of $W_2^2$ (verified vs Gaussian closed form,
  `test-wasserstein-pmed.R`; and FACT5 in the manuscript repo).
- `.pmedW_engine()` — NIE$^W$/NDE$^W$/TE$^W$/synergy/$P_{med}^W$ + **EIF-based SE** (transport layer);
  EIF-SE ≈ empirical SD verified (FACT6 / test 4).
- `.corner_laws_gcomp()` — MC parametric **g-computation** of the four corner laws (continuous 1-D Y).
- `wasserstein_pmed()` S7 generic + `data.frame` method; `WassersteinPmedResult` S7 class;
  `method = "eif" | "bootstrap"`.
- `test-wasserstein-pmed.R` — 8 passing tests: potIF vs closed form, kill-test, reduction, EIF-SE≈empirical,
  boundary endpoints.

## Promotion checklist (Claude Code)
1. **Corner-law EIF (the real DR layer).** Replace the g-comp point estimator's plug-in SE with the
   **cross-fit triply-robust corner-LAW EIF**: reuse the `.gauge_fit` cross-fit nuisance scaffold from
   `gauge_pmed.R` (π, q→density ratio, μ, η) but retain the *conditional law* of the corner pseudo-outcome,
   then compose with `.potIF` (chain) + ratio identity. This removes the empirical-$W_2$ positive bias seen
   in the reduction test and gives valid DR inference. Derivation: vault `DERIVATION-PT2.md` §A.
2. **Boundary guard** at $\mathrm{NIE}^W=0$ (non-regular): implement the regular sub-parameter test +
   split-guard; mark the CI PROVISIONAL when not rejecting (watch the Sobol A-12 over-rejection lesson).
   Spec: manuscript repo `04-wasserstein-pmed/BOUNDARY-SPEC.md`; vault `DERIVATION-PT4.md`.
3. **Entropic OT for d>1** (Sinkhorn): `.potIF` → entropic potentials; $\varepsilon_n$ schedule (PT3).
4. **S7/R/**: move class to `R/classes.R`, generic to `R/generics.R`, method to `R/methods-wasserstein.R`;
   `@export`; roxygen; `print/summary/plot` methods; `R CMD check --as-cran` clean; coverage >90%.
5. **medsim**: add corner-law DGP + `variance_share`-style estimand wiring if the sim grid runs through medsim.

## Reuse map (do NOT duplicate)
- corner nuisances + cross-fit: `gauge_pmed.R::.gauge_fit` (π, q, μ, η).
- S7 result-class pattern: `GaugePmedResult` in `gauge_pmed.R`.
- kill-test pattern: `pmed_killtest.R`.
- manuscript + closed-form truth + sims: `~/projects/research/pmed-modern/04-wasserstein-pmed/`.
