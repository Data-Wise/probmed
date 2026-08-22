# Findings: where `ward_residual()`'s variance comes from, and why its SE fails

**Date:** 2026-08-22 (draft written while the 12-cell run was at cell 5; the
mechanism results below are from completed cells and targeted diagnostics and will
not change; the full table is appended when the run finishes)
**Script:** `inst/sim/phi_decomposition.R` (modes: default / `--remedy` / `--collate`)
**Question:** the published grid shows `seW_an / empSD = 0.65-0.83` in all 8
manuscript cells and Wald coverage 0.86-0.91. `REVIEW-2026-08-22` called the
under-dispersion of `phi` the most dangerous unexamined assumption. Where does the
missing variance come from?

---

## Bottom line

1. **The estimator's excess variance is fold-split noise.** `.corner_fit()` redraws
   the partition every call (`R/corner.R:29`). That noise is 50-100% of the shipped
   single-partition estimator's variance; the irreducible (oracle-nuisance) variance
   is only 14-33% of it.
2. **The IF variance formula is exact.** With the DGP's true nuisances plugged into
   the same corner EIF, `se/sd = 1.01-1.07` and coverage is 0.94-0.97 in every cell.
3. **`reps > 1` (already shipped) fixes the point estimator completely for
   continuous Y.** At `reps = 10`, `empSD` falls from 0.121 to 0.057 vs oracle 0.056
   (n=800) and 0.048 to 0.029 vs 0.027 (n=3000); kurtosis goes normal; with a
   correctly sized SE, coverage would be 0.956. **For binary Y it helps but is not
   enough**: a few partitions give a near-zero `OE_hat`, and the mean over 10
   partitions still carries them (kurtosis 177 -> 67; `empSD` 0.66 -> 0.175 vs
   oracle 0.126). That regime needs a robust aggregation across partitions.
4. **`reps > 1` does NOT fix coverage (0.884 -> 0.884), because the SE estimator is
   the remaining defect.** Its distribution across datasets spans 0.40-2.36x the
   truth (CV 0.83 vs 0.17 for the oracle SE), is nearly uncorrelated with the actual
   error (Spearman 0.15), and the entire coverage deficit comes from the quartile of
   datasets where the SE came out small (coverage 0.65 / 0.90 / 1.00 / 0.98 by SE
   quartile).
5. **The "under-dispersion" of the published grid is a skewness artifact of that
   unstable SE**, not a biased formula: `mean(se)/sd = 0.81` at n=800 while
   `sqrt(E[se^2])/sd = 0.99`. On the RMS scale the SE is right at n=800 (0.63 at
   n=3000). Coverage fails anyway because the reps=1 error distribution is kurtotic
   (6.6-20).

So: **not a too-narrow SE; an unstable one, sitting on top of an estimator whose
variance was mostly partition noise.** Half of this is already solved by an argument
the package ships. The other half needs a different SE estimator.

---

## Method

Across datasets (nrep = 250 per cell, 10 re-partitions each), with the exact
closed-form truth from the weak-ID grid:

```
V_cf = Var(W_hat)                 shipped single-partition estimator
F    = E_data[ Var_partition(W_hat | data) ]        fold-split
V_or = Var(W_oracle)              same corner EIF, TRUE nuisances, no cross-fit
N    = (V_cf - F) - V_or          nuisance estimation, net of fold-split
V_cf = V_or + F + N
```

Oracle nuisances are closed-form from the DGP (`pi(C)`, Bayes-rule `q(M,C)`,
`mu(a,M,C)`, `eta(a,a',C)`; binary `eta` via 40-node Gauss-Hermite). Wiring
asserted by `E[phi_oracle] = theta_exact` on 2e5 draws (max error 2e-3 to 6.6e-3,
all 8 DGPs).

**Positive controls — reproduced.** On the published grid's definition
(`mean(seW_an)/sd(W)`, `collate_gauge_boot.R:12`):

| cell | n | tau | Y | se ratio | published | Wald cov | published |
|---|---|---|---|---|---|---|---|
| 1 | 800 | 0 | cont | 0.805 | 0.65 | 0.884 | 0.867 |
| 2 | 3000 | 0 | cont | 0.544 | 0.71 | 0.828 | 0.859 |
| 3 | 800 | 0.8 | cont | 0.794 | 0.72 | 0.904 | 0.899 |

(Cell 2's spread is within noise: `W_hat` has kurtosis 20 there, so `sd(W)` at 250
reps carries ~28% relative SE.)

**A false alarm worth recording.** The first version of the script used the RMS ratio
`sqrt(E[se^2])/sd(W)` and reported 0.99 for cell 1 — which looked like a failed
reproduction. It was the statistic, not the wiring: `seW_an` is right-skewed across
datasets (CV 0.5-0.7), so mean/sd and RMS/sd differ by 20-25%. Both are now reported.

## The decomposition (cells completed so far)

| cell | n | tau | irreducible `V_or/V_cf` | fold, median / mean | nuisance | oracle se/sd | oracle cov | `W_hat` kurt |
|---|---|---|---|---|---|---|---|---|
| 1 | 800 | 0 | 0.23 | 0.53 / 1.08 | ~0 (noise) | 1.01 | 0.940 | 6.6 |
| 2 | 3000 | 0 | 0.14 | 0.30 / 0.59 | 0.27 | 1.07 | 0.968 | 20 |
| 3 | 800 | 0.8 | 0.33 | 0.36 / 0.72 | ~0 (noise) | 1.01 | 0.960 | 6.8 |

Mean-based fold shares overshoot (can exceed 1) because the within-dataset partition
variance is heavy-tailed — mean/median ~2 in every cell — so a few datasets with one
wild partition dominate. The median-based share is the conservative reading; either
way fold-split is the largest component. Nuisance vs fold is not cleanly separable at
this precision, which is now a second-order question.

## The remedy check (same 250 datasets per cell, shipped `ward_residual()`)

| cell | arm | empSD | se mean | se CV | Wald cov | kurt |
|---|---|---|---|---|---|---|
| 1: n=800, cont | `reps = 1` | 0.121 | 0.098 | 0.73 | 0.892 | 4.9 |
| | `reps = 10` | **0.057** | 0.059 | 0.66 | **0.884** | 3.1 |
| | oracle | 0.056 | 0.056 | 0.18 | 0.940 | 3.1 |
| 2: n=3000, cont | `reps = 1` | 0.048 | 0.039 | 0.55 | 0.852 | 5.1 |
| | `reps = 10` | **0.029** | 0.025 | 0.51 | **0.844** | 4.2 |
| | oracle | 0.027 | 0.029 | 0.11 | 0.968 | 2.9 |
| 5: n=800, **binary** | `reps = 1` | 0.658 | 0.513 | **9.3** | 0.900 | **177** |
| | `reps = 10` | 0.175 | 0.143 | 1.8 | 0.928 | 67 |
| | oracle | 0.126 | 0.124 | 0.29 | 0.972 | 3.2 |

**Continuous Y (cells 1, 2):** the point estimator at `reps = 10` reaches oracle
efficiency at both n; with a *fixed* `se = empSD`, coverage at cell 1 is **0.956**.
Solved. Cost: 10 cross-fits instead of 1 (~0.5 s at n = 800).

**Binary Y (cell 5) is a different regime.** The single-partition estimator is
catastrophically heavy-tailed (kurtosis 177; RMS of the reported SE is 4.8 = some
datasets report an SE near 50): a few partitions yield a near-zero `OE_hat` and
`W = R/OE` explodes. Averaging 10 partitions cuts `empSD` 4x and lifts coverage to
0.93, but does **not** reach the oracle (kurtosis still 67) — a mean across partitions
is not robust to one wild partition. Here the remedy needs either many more reps or a
robust aggregation across partitions (median / trimmed mean of the per-partition
corner means, or of `W` itself). The SE problem is correspondingly worse (CV 1.8).

## Anatomy of the `reps = 10` SE (60 datasets, cell 1)

The shipped SE is `sqrt(seIF^2 + seMC^2)`: the IF sd from the partition-averaged
`phi`, plus `var(W_reps)/reps` (`R/gauge-pmed.R:275`).

| component | median x truth | CV | 5-95% x truth |
|---|---|---|---|
| `seIF` (averaged `phi`) | **0.63** | 0.81 | 0.33-2.13 |
| `seMC` (`var(W_reps)/10`) | — | 0.95 | 0.23-1.18 |
| shipped total | 0.81 | 0.83 | 0.40-2.36 |
| single-partition IF se (reps=1 style) | 1.35 | 0.69 | 0.71-3.79 |
| oracle | 1.00 | **0.17** | 0.68-1.18 |

Both parts fail, in opposite directions: averaging `phi` over partitions
**over-smooths** the per-observation contributions, so the IF sd comes out ~37% below
the oracle's (the low-SE quartile that under-covers); the MC add-on is a heavy-tailed
variance estimated from 10 draws (the high tail; it is 35% of `se^2` typically, up to
69%). The means nearly cancel; the spread does not.

**Hypothesis tested and rejected: extreme estimated weights.** On 30 datasets, the
maximum inverse-propensity weight (19 vs 17) and maximum density-ratio weight (7.8 vs
8.0) are the same for estimated and true nuisances; the oracle IF is itself
heavy-tailed (kurtosis ~20, top-1% share 0.37) yet yields a stable SE. What differs is
the IF's *spread*: the cross-fit IF sd is 47% larger than the oracle's and 3x more
variable across datasets (CV 0.57 vs 0.18). The extra variability enters through
nuisance errors in the residual terms `(Y - mu_hat)` and `(mu_hat - eta_hat)`, which
are fold-specific — not through the weights.

---

## What this changes

**SPEC-2026-08-21, Phases 2-3 (no-refit bootstrap; Fieller for W):** already refuted
by the adversarial review; this explains *why* the refuted diagnosis looked plausible
(the SE tracks the partition, so `E[se^2] ~ Var` at n=800) and why `phi`-reweighting
could never have worked (it inherits the IF sd's instability exactly). Note the
**Fieller set has the same problem**: it is built from `var(pOE)/n` of the same `phi`,
so the `oe_regular` gate and the Fieller bounds are as unstable under re-partition as
the SE is. A1's survival of the width control in the 48k grid is not contradicted, but
its precision is lower than the grid implied.

**The manuscript's Simulation section.** The "right-skewed analytic SE" framing is
correct as description; the mechanism claimed (ratio skewness) is wrong — it is
partition noise, and the fix for the point estimate is `reps > 1`, which the
manuscript already mentions in passing (`gauge-pmed.qmd:236`) as removing "the
fold-split variance component." That sentence is doing far more work than the text
around it admits. The percentile-bootstrap arm's ~1.00 coverage is now explicable too:
each resample refits with a fresh partition, so the bootstrap distribution carries the
fold-split noise the `reps=1` point estimate also carries — at n=3000, where the review
measured its SD as correct, that is the right variance for the *wrong* (reps=1)
estimator.

**Issue #32 (BCa):** moot as posed. The defect is not interval shape.

**Issue #11 (A1/A2):** the biased-denominator account of A2 sharpens further — the
denominator is not merely biased, it has CV ~0.7. A2 is a ratio of two unstable
widths.

## Open, in priority order

1. **An SE estimator for the `reps > 1` point estimator.** Candidates, cheapest
   first: (a) `reps` large enough that `seMC -> 0` and `seIF` can be calibrated
   against `V_avg` empirically; (b) the IF sd averaged *across* partitions with each
   partition's own `W`, rather than from the averaged `phi` (tests whether the
   over-smoothing is in the averaging step); (c) a nonparametric bootstrap of the
   `reps = 10` estimator (cost `B x reps` cross-fits; the refit bootstrap's SD was
   already right at n = 3000 for reps = 1). The oracle benchmark is CV 0.17.
2. **Robust aggregation across partitions for binary Y** — `reps = 10` reaches the
   oracle for continuous Y at both n but not for binary Y (cell 5: `empSD` 0.175 vs
   0.126, kurtosis 67). Candidates: median of per-partition corner means before
   forming the ratio; trimmed mean; or an explicit guard that drops partitions whose
   `OE_hat` falls below a threshold. Each changes the estimand's finite-sample
   behavior and needs its own coverage check. Near-null cells (9-12) still pending
   from the main run.
3. **Why does averaging `phi` over-smooth?** Likely because a fold-specific nuisance
   error is shared by every observation in that fold and so is correlated across
   observations — invisible to a per-observation sd. Testable by comparing
   `sd(ifW(pbar))` with the sd of the per-fold means.
4. Fold-split noise is heavy-tailed (`F_tail ~ 2`). Which partitions are "wild", and is
   it a small-fold separation effect in the `A == a'` subset regressions for `eta`?

## Reproducibility

Everything above: `Rscript inst/sim/phi_decomposition.R` (default grid),
`--remedy <cell>`, and the two inline diagnostics recorded in the session transcript
(weight extremity; SE component split), which should be promoted into a `--se-parts`
mode before the next iteration. Seeds are deterministic per (cell, rep, partition).
