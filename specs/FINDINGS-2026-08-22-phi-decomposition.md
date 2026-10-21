# Findings: where `ward_residual()`'s variance comes from, and why its SE fails

> **SUPERSEDED 2026-08-22 (later the same day) — the "fold-split noise" was a bug.**
> Every package-estimator number below (the `V_cf` column, `F`, `N`, `seW_an`,
> the `reps` remedy rows, the binary "transition regime", the SE-candidate rounds
> in §Open item 1) was produced with `.corner_fit()`'s corner weights built as
> `ifelse(z == 1, p1, 1 - p1)`. `ifelse()` returns a result shaped like its
> *test*; `z` is a scalar; so every row in a fold was weighted by the **first test
> row's** propensity (and `q`). A different first row per partition is a
> different constant weight per partition — that is the "fold-split noise", and
> the "dataset-level nuisance-fit quality" of the SE rounds. The oracle columns
> (`V_or`, exact truth, `phi_oracle`) were computed by this script's own correct
> EIF and stand. With the weights fixed, the single-partition estimator's `empSD`
> is 0.063 vs oracle 0.056 (cell 1) and 0.121 vs 0.122 (cell 5), and the
> shipped analytic Wald interval covers 0.936 / 0.972 at `reps = 1` — see
> `inst/sim/se_shipped_check.R`, `inst/sim/results/se_shipped_postfix.csv`, and
> NEWS 0.3.0.9000. What survives of the bottom line: item 2 (the IF formula is
> exact under true nuisances) and the near-null non-regularity (A1's regime).
> What does not: items 1, 3, 4, the three regimes, the remedy section, and §Open
> items 1-4 as posed. Kept unedited below as the record of how the bug was found.

**Date:** 2026-08-22 (final: all 12 cells at nrep = 250 x 10 partitions, plus
remedy checks on cells 1, 2, 5 and two targeted diagnostics; summary table promoted to
`inst/sim/results/phi_decomp_summary.csv`)
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

6. **Two mechanisms, separated cleanly by `oe_snr`.** In the 8 strong-ID cells the
   oracle is exact and the excess variance is cross-fitting noise. In the 2 near-null
   cells (`oe_snr` ~1-1.7) **the oracle SE explodes as well** (se/sd 591 and 333;
   coverage 0.98-1.00 from enormous intervals) and fold-split is 1-3% of the variance:
   the ratio `W = R/OE` is non-regular there and no nuisance quality helps. That is
   A1's regime, and it is the manuscript's own Application-section argument (Fieller
   set unbounded; report the unnormalized effects). The manuscript's *Simulation*
   section, by contrast, lives entirely in the first regime, where the story is
   partition noise.

So: **not a too-narrow SE; an unstable one, sitting on top of an estimator whose
variance was mostly partition noise** — in the regime where the estimand is regular.
Half of this is already solved by an argument the package ships. The other half needs
a different SE estimator. Outside that regime the problem is the estimand, which A1
already flags.

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

| cell | n | tau | Y | se ratio | published | median ratio | Wald cov | published |
|---|---|---|---|---|---|---|---|---|
| 1 | 800 | 0 | cont | 0.805 | 0.65 | 0.648 | 0.884 | 0.867 |
| 2 | 3000 | 0 | cont | 0.544 | 0.71 | 0.466 | 0.828 | 0.859 |
| 3 | 800 | 0.8 | cont | 0.794 | 0.72 | 0.672 | 0.904 | 0.899 |
| 4 | 3000 | 0.8 | cont | 0.732 | 0.83 | 0.645 | 0.880 | 0.913 |
| 5 | 800 | 0 | bin | 1.35* | 0.72 | 0.287 | 0.932 | 0.879 |
| 6 | 3000 | 0 | bin | 0.763 | 0.74 | 0.654 | 0.856 | 0.869 |
| 7 | 800 | 0.8 | bin | 0.707 | 0.74 | 0.535 | 0.888 | 0.892 |
| 8 | 3000 | 0.8 | bin | 0.638 | 0.67 | 0.508 | 0.860 | 0.874 |

Reproduced in 7 of 8 cells, within the noise of a 250-rep `sd(W)` under kurtosis
7-20 (~20-30% relative SE). *Cell 5 is the explosive cell (kurtosis 130): a handful
of datasets report SEs ~50x the empirical SD, so any mean-based ratio is meaningless
there — the median ratio (0.29) is the readable number, and the published 0.72 at
2000 reps was itself a mean over the same kind of tail.

**A false alarm worth recording.** The first version of the script used the RMS ratio
`sqrt(E[se^2])/sd(W)` and reported 0.99 for cell 1 — which looked like a failed
reproduction. It was the statistic, not the wiring: `seW_an` is right-skewed across
datasets (CV 0.5-0.7), so mean/sd and RMS/sd differ by 20-25%. Both are now reported.

## The decomposition — all 12 cells

| cell | regime | n | tau | Y | s | `oe_snr` med | irreducible `V_or/V_cf` | fold (median) | `F_tail` | oracle se/sd | oracle cov | `W_hat` kurt | SE CV |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | manuscript | 800 | 0 | cont | 1 | 8.3 | 0.23 | 0.53 | 2.1 | 1.01 | 0.940 | 6.6 | 0.72 |
| 2 | manuscript | 3000 | 0 | cont | 1 | 17 | 0.14 | 0.30 | 2.0 | 1.07 | 0.968 | 20 | 0.58 |
| 3 | manuscript | 800 | 0.8 | cont | 1 | 9.8 | 0.33 | 0.36 | 2.0 | 1.01 | 0.960 | 6.8 | 0.49 |
| 4 | manuscript | 3000 | 0.8 | cont | 1 | 19.5 | 0.24 | 0.36 | 2.0 | 1.08 | 0.960 | 6.8 | 0.49 |
| 5 | manuscript | 800 | 0 | bin | 1 | 4.5 | **0.05** | 0.11 | 3.1 | 1.02 | 0.972 | **130** | 10.7 |
| 6 | manuscript | 3000 | 0 | bin | 1 | 8.7 | 0.23 | 0.45 | 3.2 | 0.98 | 0.944 | 10 | 0.67 |
| 7 | manuscript | 800 | 0.8 | bin | 1 | 5.3 | 0.16 | 0.34 | 79 | 1.06 | 0.968 | 11 | 1.0 |
| 8 | manuscript | 3000 | 0.8 | bin | 1 | 10.5 | 0.16 | 0.30 | 2.8 | 1.01 | 0.940 | 18 | 0.93 |
| 9 | near-null | 800 | 0.4 | cont | 0.2 | 1.7 | 0.35 | **0.03** | 141 | **591** | 0.996 | 90 | 9.8 |
| 10 | near-null | 800 | 0.4 | bin | 0.2 | 1.0 | 3.4 | **0.01** | 4420 | **333** | 0.980 | 70 | 7.4 |
| 11 | intermediate | 800 | 0.4 | cont | 0.5 | 4.3 | 0.30 | 0.35 | 6.2 | 1.01 | 0.972 | 13 | 1.3 |
| 12 | intermediate | 800 | 0.4 | bin | 0.5 | 2.5 | 0.18 | 0.26 | 48 | 1.52 | 0.988 | 27 | 5.5 |

(`SE CV` = across-dataset CV of the shipped analytic SE. Mean-based fold shares are
omitted: the within-dataset partition variance is heavy-tailed — `F_tail` = mean/median
~2 in the well-behaved cells and 50-4000 in the explosive ones — so the median share is
the only readable one. Full columns in `inst/sim/results/phi_decomp_summary.csv`.)

### Three regimes

**Strong identification, regular (cells 1-4, 6-8, 11; `oe_snr` ≳ 4).** Oracle exact
(se/sd 0.98-1.08, coverage 0.94-0.97). Irreducible variance 14-33% of the shipped
estimator's; fold-split the largest component (30-53% by median). `W_hat` kurtosis
6.6-20; Wald coverage 0.83-0.90. **This is where the manuscript's Simulation section
lives, and the story is partition noise.**

**Near-null, non-regular (cells 9-10; `oe_snr` 1-1.7).** The oracle SE explodes (se/sd
591, 333) and the oracle over-covers (0.98-1.00) through enormous intervals; fold-split
is 1-3% of the variance. Even with true nuisances, `OE_hat` lands near zero in some
datasets and `W = R/OE` blows up. **Nothing about cross-fitting matters here; the
estimand is non-regular.** This is exactly A1's regime (`oe_snr < 2` -> `oe_regular =
FALSE`), and it is the manuscript's own Application-section conclusion: report the
unnormalized effects, the Fieller set is unbounded.

**Transition (cell 5: binary, tau=0, n=800, `oe_snr` 4.5; cell 12, `oe_snr` 2.5).**
Oracle fine, but the cross-fit estimator has 20x the oracle variance (cell 5), and the
median fold share is small (0.11): the explosion is partly *dataset*-level (nuisance
error pushes `OE_hat` toward zero for the whole dataset, every partition) rather than
*partition*-level. That is why `reps = 10` halved `empSD` there but did not reach the
oracle. Binary outcomes with a small true `OE` on the probability scale sit closest to
this edge.

Nuisance vs fold is not cleanly separable at this precision in the regular regime;
that is now a second-order question.

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

1. **An SE estimator for the `reps > 1` point estimator — RESOLVED (2026-08-22,
   three rounds, 250 datasets per cell, `reps = 10`, scratch scripts
   `se_candidates{,2,3}.R` in the session scratchpad — promote into `inst/sim/` with
   the code change).** Score = median(SE)/empSD(W), CV of the SE across
   datasets, Wald coverage; the oracle-nuisance SE is the benchmark.

   *Round 1 (cell 1): every partition-aggregation candidate fails.* Shipped total
   0.81 / CV 0.62 / cov 0.896; per-partition IF mean 1.33 / 0.57 / 0.984; per-partition
   median 1.23 / 0.61 / 0.968; fold-clustered sandwich, median over partitions 0.96 /
   0.54 / 0.948; oracle 0.92 / 0.18 / 0.940. The best candidate gets the median and
   the coverage right only because its per-dataset errors (0.5-2.4x truth) balance —
   CV 0.54-0.66 across the whole family. So the instability is a property of the
   *dataset's* nuisance fits, not of how partitions are combined. Side result: the
   fold-clustered SE is *smaller* than the iid one (1.19 vs 1.33), i.e. within-fold IF
   contributions are negatively correlated — the shared-fold-error over-smoothing story
   in item 3 below is wrong in its simple form.

   *Round 2 (cell 1): the IF variance from full-sample nuisance fits works.* Keep the
   cross-fit `reps = 10` point estimate; for the SE only, refit the three nuisance GLMs
   (and the `eta` projections) on all n, form the corner `phi` in-sample, and take
   `sd(ifW(phi_full, W_crossfit)) / sqrt(n)`. Result 0.88 / **0.20** / 0.936 vs oracle
   0.92 / 0.18 / 0.940. Winsorizing the IF contributions (1/99%) is ruled out: it drops
   even the *oracle* to 0.72 / cov 0.85, so the IF tails are signal.

   *Round 3: holds outside cell 1.* Intermediate (s = 0.5, tau = 0.4, cont):
   1.00 / 0.23 / 0.940 vs oracle 1.00 / 0.25 / 0.964. Binary strong-ID (cell 5):
   0.91 / 0.31 / 0.948 vs oracle 0.92 / 0.31 / 0.972. In both, the shipped SE has CV
   1.7-1.9 (its `seMC` add-on explodes on unstable partitions). The full-sample SE's CV
   equals the oracle's in all three cells.

   Mechanism: the cross-fit IF evaluates each observation under a 4/5-sample nuisance
   fit, and the IF *variance* inherits that fit's noise — dataset-specific, so no
   partition aggregation removes it; the full-sample fit's noise is at the oracle
   rate. The point estimate keeps cross-fitting (where the bias protection lives).
   Caveat: in-sample nuisances can understate the variance if the nuisance models
   overfit — negligible for these GLMs at n = 800 (0.88-1.00 vs oracle 0.92-1.00);
   re-check if flexible learners are ever plugged into `.corner_fit`.

   Shipping it is a code change in `R/corner.R` / `R/gauge-pmed.R` (new `se_method`
   value or the default for `reps > 1`) with a coverage test against the saved tables.
   Pending: the nonparametric bootstrap of the `reps = 10` estimator (B = 50, 60
   datasets, cell 1) as a cost comparison — it is B x reps = 500 cross-fits per
   dataset against one extra GLM fit for the full-sample IF.

   Discrepancy to re-check under item 2: in this 250-dataset run the `reps = 10`
   point estimator's `empSD` in cell 5 was 0.125 (oracle 0.122), not the 0.175 the
   remedy run reported — different seeds and rep count; one of the two is a tail
   event.
2. **Robust aggregation across partitions for binary Y** — `reps = 10` reaches the
   oracle for continuous Y at both n but not for binary Y (cell 5: `empSD` 0.175 vs
   0.126, kurtosis 67). Candidates: median of per-partition corner means before
   forming the ratio; trimmed mean; or an explicit guard that drops partitions whose
   `OE_hat` falls below a threshold. Each changes the estimand's finite-sample
   behavior and needs its own coverage check. Near-null cells (9-12) still pending
   from the main run.
3. **Why does averaging `phi` over-smooth?** — CLOSED by item 1. The per-observation
   sd is not hiding positively correlated shared-fold error (the fold-clustered SE is
   smaller than the iid SE, so the within-fold correlation is negative); the
   under-dispersion is fold-fit noise entering the IF itself, and it disappears when
   the IF is formed from full-sample nuisance fits.
4. Fold-split noise is heavy-tailed (`F_tail ~ 2`). Which partitions are "wild", and is
   it a small-fold separation effect in the `A == a'` subset regressions for `eta`?

## Reproducibility

Everything above: `Rscript inst/sim/phi_decomposition.R` (default grid),
`--remedy <cell>`, and the two inline diagnostics recorded in the session transcript
(weight extremity; SE component split), which should be promoted into a `--se-parts`
mode before the next iteration. Seeds are deterministic per (cell, rep, partition).
