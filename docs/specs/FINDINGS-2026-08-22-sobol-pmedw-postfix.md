# Findings: `sobol_pmed()` and `pmedW_dr()` after the corner-EIF weight fix

**Date:** 2026-08-22 (night), PLAN-2026-08-22-postfix-followups items 1-2.
**Why:** both estimators build on the corner machinery fixed in PR #34 (`a73efc8`);
their documentation and their one validation each were written against the buggy
weights.

## 1. `sobol_pmed()` — six cells against exact truth

**Script:** `sobol_shipped_check.R` (session scratchpad; `.sobol_fit()` called as
`sobol_pmed()` calls it: `boundary_test = "split"`, `procedure = "B"`, K = 5,
`pd = pm = 0.5`), 250 datasets per cell for the analytic arm, the first 100 also
with `se_method = "bootstrap"`, B = 100. Truth = `sobol_from_theta` algebra on
`truth(s, tau, FALSE)$theta` from `inst/sim/phi_decomposition.R`; oracle se =
the same delta-method IF on `phi_oracle` (true nuisances). Collated table:
`sobol_postfix_summary.csv` (promote with the code change). Monte-Carlo
half-width on a coverage: ~0.03 (250) / ~0.04 (100).

| cell | n | `P` true | `Delta_m` true | bias | empSD / oracle SD | se/empSD (CV) | oracle se/empSD | boundary rate |
|---|---|---|---|---|---|---|---|---|
| s=1, tau=0 | 800 | 0.414 | 0.42 | +0.007 | 0.115 / 0.112 | 0.97 (0.17) | 0.94 | 0 |
| s=1, tau=0 | 3000 | 0.414 | 0.42 | 0.000 | 0.055 / 0.055 | 1.05 (0.09) | 1.03 | 0 |
| s=1, tau=0.8 | 800 | 0.419 | 0.66 | -0.001 | 0.082 / 0.080 | 0.99 (0.11) | 0.97 | 0 |
| s=1, tau=0.8 | 3000 | 0.419 | 0.66 | -0.002 | 0.041 / 0.041 | 1.02 (0.06) | 1.01 | 0 |
| s=0.5, tau=0.4 (intermediate) | 800 | 0.421 | 0.24 | +0.002 | 0.184 / 0.182 | 0.96 (0.28) | 0.92 | 0.50 |
| s=0.2, tau=0.4 (near-null) | 800 | 0.417 | 0.089 | +0.020 | 0.337 / 0.333 | 0.91 (1.44) | 0.92 | 0.98 |

The estimator is at oracle efficiency in every cell (empSD within 3% of the
oracle's) and its analytic se matches the oracle's ratio. The near-null se CV of
1.44 is the degeneracy of `phi_P = 2 c_m Delta_m phi_{Delta_m}` as `Delta_m -> 0`
— the estimand's, shared with the oracle.

**Coverage of the four interval constructions (analytic se):**

| cell | Wald for `P` | oracle Wald | Procedure A (gated) | **Procedure B1 (default)** | B2 (split) | B1 / Wald width | B2 / B1 width |
|---|---|---|---|---|---|---|---|
| s=1, tau=0, n=800 | 0.912 | 0.908 | 0.912 | **0.956** | 0.932 | 1.06 | 1.42 |
| s=1, tau=0, n=3000 | 0.944 | 0.948 | 0.944 | **0.956** | 0.972 | 1.06 | 1.43 |
| s=1, tau=0.8, n=800 | 0.928 | 0.924 | 0.928 | **0.996** | 0.968 | 1.38 | 1.41 |
| s=1, tau=0.8, n=3000 | 0.940 | 0.940 | 0.940 | **0.988** | 0.980 | 1.38 | 1.44 |
| intermediate | 0.868 | 0.876 | 0.904 | **0.968** | 0.944 | 1.17 | 1.42 |
| near-null | 0.740 | 0.740 | 0.916 | **0.952** | 0.988 | 1.42 | 1.67 |

- **The Wald interval for `P` fails exactly where the oracle's does** (0.87 /
  0.74 at the transition / boundary, 0.91-0.94 in the regular cells): the
  square map `delta -> c_m delta^2 / V_T`, not the se. Nothing to fix in the se.
- **Procedure B1, the default, is valid everywhere** (0.952-0.996), conservative
  under a strong interaction (tau = 0.8: 0.99-1.00 at 1.38x the Wald width) and
  nominal at the boundary (0.952 with 98% of draws gated).
- **Procedure A under-covers at the transition** (0.904 / 0.916): the
  Leeb-Potscher pre-test effect the docs describe. Docs right.
- **B2 adds nothing** (1.4x wider than B1, coverage not better in any cell).
  Docs right ("experimental ... did not improve coverage").

**The `Delta_m` scale, where the docs' near-null story lives:**

| cell | se(`Delta_m`)/empSD (CV) | `Delta_m` Wald coverage |
|---|---|---|
| s=1, tau=0.8, n=800 | 0.99 (0.08) | 0.952 |
| intermediate | 1.00 (0.08) | 0.960 |
| near-null | 0.99 (0.08) | 0.960 |

The documented "residual analytic-se bias: se(`Delta_m`) ~0.8x its sampling SD
near the null, so the default Procedure-B interval covers ~0.85; remedy
`se_method = "bootstrap"` (valid, ~1.25x conservative, coverage ~0.97)" and the
"fold-split Monte-Carlo variance is ~80% of `Var(Delta_m_hat)` near the null;
remedy `reps > 1`" (both in `?SobolPmedResult`, "Coverage near the boundary",
items (2)) **were the weight bug**. Post-fix the `Delta_m` se is exact at every
effect size with CV 0.08, and `reps > 1` has nothing left to remove (empSD is
at the oracle's).

**Bootstrap arm (B = 100, 100 datasets):** bootstrap se / empSD 0.94 / 1.05 /
0.95 / 1.02 in the regular cells, **0.90 intermediate, 0.84 near-null** (CV
0.08-0.21); Wald-with-bootstrap-se 0.93 / 0.94 / 0.94 / 0.95 / 0.90 / 0.88; B1
with bootstrap se 0.96 / 0.94 / 1.00 / 1.00 / 0.98 / 0.97. It is now slightly
*anti*-conservative where the docs call it conservative, and never better than
the analytic default.

### What this changes in `sobol_pmed` / `SobolPmedResult`

1. `?SobolPmedResult` "Coverage near the boundary — (2) The near-null `Delta_m`
   standard error": replace the three "pinned facts" (fold-split 80%, normal
   shape, 0.8x residual bias) with: the analytic se of `Delta_m` is exact
   post-fix (table above); the only near-null pathology is the Wald interval
   for `P` itself, which Procedure B already avoids. Item (1) on Procedure A
   stands.
2. "Practical guidance": drop "pass `reps > 1` ... and `se_method = "bootstrap"`
   for a valid interval"; the analytic Procedure B default is valid everywhere
   measured. Keep `reps` and the bootstrap as options.
3. `se_method` parameter docs (class and generic) and the boundary `message()`
   text in `.sobol_fit()` ("the analytic se_Dm is anti-conservative here ...
   pass se_method = 'bootstrap'"): rewrite; the message should say Procedure B
   is reported and is valid at the boundary.
4. `ci_B2` stays "experimental, did not improve coverage" — confirmed.
5. Promote the script and CSV to `inst/sim/`.

Code + Rd change: needs a worktree.

## 2. `pmedW_dr()` — pending

Script `pmedw_shipped_check.R`: point estimate pre- vs post-fix against a
2e6-draw simulated truth; `wasserstein_pmed(method = "dr")` bootstrap coverage
on a subset. Results to be appended.
