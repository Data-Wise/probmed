# T5 feasibility: a calibrated MBCO interval near the composite null (2026-10-14)

**Status:** NO-GO. Keep the χ²₁ inversion and the `near_double_null` label.
This records the decision made at Phase 8 Checkpoint B (`tasks/plan-mbco-null.md`).

Provenance: an external feasibility audit (DeepSeek) produced the first draft.
It was then reviewed and corrected against the lattice data and the papers'
abstracts. The corrections are listed at the end.

## 1. Main argument: the conservativeness is not confined to the null test

An interval for P_med inverts a family of tests, one for each p*. A
composite-null method changes only the test at p* = 0.5 (ab = 0). The lattice
(`inst/sim/results/calibration/mbco_null/scores.csv`, χ²₁ rule, pooled over
N ∈ {100, 1000} and with or without the covariate) shows overcoverage at
points where the true P_med is **not** 0.5. There the interval accepts the
truth through the test at the true p* ≠ 0.5, which a null-only method leaves
unchanged:

| max(h_a, h_b), both paths nonzero | coverage range |
|---|---|
| 0.5 | 0.9955–0.9965 |
| 1 | 0.986–0.9975 |
| 2 | 0.9615–0.9935 |
| 4 | 0.9415–0.9755 |

So a better null test cannot remove the conservativeness at (h_a, h_b) = (1, 1)
or (2, 2). Grafting a sharper p* = 0.5 test onto the χ²₁ family has two
effects. It changes only whether 0.5 is inside the set, and it can cut a hole
there: 0.5 is excluded while its neighbors stay in. A genuinely less
conservative interval near the origin would need a test family that is
uniformly valid over every p* near 0.5. None of the candidates below supplies one.

## 2. Candidate methods, at the level checked

Each claim below was checked against the arXiv abstract only, not the full
text.

| Method | What the abstract establishes | Interval for all p*? |
|---|---|---|
| Miles & Chambaz, arXiv 2107.07575 | Tests the composite null "if either regression coefficient is zero"; one approach is a closed-form test, minimax optimal for local power; preserves level α | No; the null only |
| van Garderen & van Giersbergen, arXiv 2012.11342 | Tests H0: θ1θ2 = 0; a similar test exists only if the size is the reciprocal of an integer, and is unique but has objectionable properties; proposes a powerful, easy-to-use test | No; the null only |
| He, Song & Xu, JRSS-B 2024 (arXiv 2310.11630) | Adaptive bootstrap framework for composite nulls in the mediation pathway; type I error control under the composite null, with improved power | No; the null only |

The abstracts make no claims about covariates, finite-N behavior or
computational cost. This report makes none either.

**Covariates are not an obstacle to these methods.** In the single-mediator
model the likelihood factors into the M model and the Y | M model. The
information matrix is therefore block-diagonal, so â and b̂ are asymptotically
independent with or without covariates. The real open questions for adopting
any of them are finite-N behavior (N = 100) and plugging in estimated standard
errors. These matter only for a better **test** of mediation, not for the
interval.

## 3. Why the T4 rules failed, and why the same applies here

At p* = 0.5 the profile LR equals min(LR_a0, LR_b0). A single null with the
other path strong has min(LR) ≈ LR of the null branch, which is ~χ²₁. Any rule
that rejects when min(LR) exceeds a cutoff below the χ²₁ value therefore
over-rejects at single nulls. Off the null, LR(truth) is no longer
min(branch). The T4 switch and blend rules fixed the origin but undercovered
at nearby alternatives (minimum coverage 0.84–0.92; see the Phase 8 addendum
in `CALIBRATION-2026-10-14.md`). A cutoff tuned to the null test does not
transfer to the test family.

## 4. Decision and what would reopen it

- **Shipped:** the χ²₁ interval, which never undercovered on the lattice
  (minimum 0.9415), and the `near_double_null` label with a print note.
- **Possible separate feature (not T5):** a composite-null **test** of mediation
  (for example Miles–Chambaz), reported beside the interval, never by editing
  the interval. Its own gate would be size ≤ α at every null lattice point, and
  power reported at the non-null points.
- **What would reopen T5:** a published construction giving uniformly valid tests
  of P_med = p* for every p* near 0.5. It must pass the existing lattice
  acceptance rule: minimum coverage ≥ 0.94 at both N, with and without the
  covariate, and origin ≤ 0.97.

**Harness gap:** the plan called for storing the branch Wald z values. The
stored columns are `lr_true, lr_a0, lr_b0, ci_lower, ci_upper`. |z| = √LR can
be recovered, but the sign cannot. Any future null-test prototype needs the
signed z values added to `calib_mbco_null.R`.

## Corrections to the external draft

1. **Removed** its prototype: "reject when min(lr_a0, lr_b0) > c_α". That rule
   is the joint-significance test, not Miles–Chambaz's minimax test. With
   c_α below the χ²₁ cutoff it has size > α at single nulls (§3).
2. **Reversed** its claim that covariates violate a diagonal-covariance
   assumption. The likelihood factors, so â and b̂ are asymptotically
   independent (§2).
3. **Cut** the claims that go beyond the abstracts: diagonal covariance,
   "computationally intensive", "validity at N = 100 questionable",
   implementation details.
4. **Corrected** its statement that the per-replication statistics were
   unavailable. They are committed in
   `inst/sim/results/calibration/mbco_null/replications.csv.gz`.
5. **Added** §1, the lattice evidence, as the main argument.

## Offline reproducibility checkpoint (2026-11-23)

[Issue #50](https://github.com/Data-Wise/probmed/issues/50#issuecomment-7123567442) closed after [PR #70](https://github.com/Data-Wise/probmed/pull/70) merged as [`dbb4292`](https://github.com/Data-Wise/probmed/commit/dbb429292bbf5bc633906840f19076e2dedc57d5).
The [audit spec](../specs/SPEC-2026-10-15-t5-prototype-audit.md) defines the
prototype scoring and the frozen pointwise Monte Carlo gates. The first
checkpoint uses the committed 100-cell lattice and its T4 score table.

```sh
lattice=inst/sim/results/calibration/mbco_null
Rscript inst/sim/calibration/calib_mbco_null_t5_baseline.R \
  "$lattice/replications.csv.gz" "$lattice/scores.csv" /tmp/t5-baseline.csv
```

All 400 baseline rows reproduce, with identical acceptance and replication
counts and exact 95% binomial interval endpoints within `1e-12` of the saved
values. Each cell has 2,000 replications. The chi-square rule's minimum
coverage at N = 1,000 without the covariate is 0.9440.

| Baseline rule | Minimum coverage gate | Origin ceiling | Overall |
|---|---|---|---|
| `chi2` | INCONCLUSIVE | FAIL | FAIL |
| `switch_crit` | FAIL | INCONCLUSIVE | FAIL |
| `switch_logN` | FAIL | PASS | FAIL |
| `smooth` | FAIL | FAIL | FAIL |

Branch-specific reconstruction gives 157 LR-versus-absolute-t ordering
disagreements among 200,000 rows. In task 1, replication 558, LR selects
the mediator branch and absolute t selects the outcome branch. An offline
mapping that uses the t order must therefore take its degrees of freedom
from the outcome regression for that row.

Exclude `MC21_t_all` from interval evidence. Section 3.1 of
[Miles and Chambaz (v4)](https://arxiv.org/html/2107.07575v4) defines
the rejection region for a fixed composite null in a normal experiment.
The result supplies no finite-sample inversion over candidate P_med
values. Reconstruction of the absolute t statistics resolves the branch
mapping and leaves that theoretical gap open.

## Offline prototype scores (2026-11-23)

The original October 15 pasted scorer was recovered and rerun on the frozen
replications. The committed scorer reproduced the acceptance counts and
replication counts for all 1,500 retained rule-cell rows, including all 400
T4 baseline rows. Exact 95% Clopper–Pearson endpoints matched within
`1e-12`. The saved table has 15 rules × 100 cells, with 2,000 replications
per cell. It omits `MC21_t_all` under the spec's exclusion clause.

```sh
lattice=inst/sim/results/calibration/mbco_null
Rscript inst/sim/calibration/calib_mbco_null_t5_score.R \
  "$lattice/replications.csv.gz" /tmp/t5_prototype_scores.csv
cmp /tmp/t5_prototype_scores.csv "$lattice/t5_prototype_scores.csv"
```

The `_05` rules use the prototype decision when the true P_med is 0.5
and the baseline χ²₁ decision at the 64 alternative cells. Their null-cell
values are nonrejection rates. The alternative-cell values for these
rules describe inclusion of the true P_med under the baseline decision.
The `_all` rules apply their cutoff to the profile LR at the true P_med
throughout the lattice. These are research diagnostics. The two gates
use the pointwise exact binomial intervals in the score table.

| Rule | Null min | Alternative min | Origin range | Minimum gate | Origin gate | Overall |
|---|---:|---:|---:|---|---|---|
| `chi2` | 0.9470 | 0.9415 | 0.9960–0.9985 | Inconclusive | Fail | Fail |
| `switch_crit` | 0.9190 | 0.8630 | 0.9625–0.9655 | Fail | Inconclusive | Fail |
| `switch_logN` | 0.8600 | 0.8355 | 0.9510–0.9575 | Fail | Pass | Fail |
| `smooth` | 0.9470 | 0.9080 | 0.9960–0.9985 | Fail | Fail | Fail |
| `MC21_t_05` | 0.9380 | 0.9415 | 0.9440–0.9560 | Inconclusive | Pass | Inconclusive |
| `MC21_LR_05` | 0.9390 | 0.9415 | 0.9415–0.9565 | Inconclusive | Pass | Inconclusive |
| `MC21_LR_all` | 0.9390 | 0.7900 | 0.9415–0.9565 | Fail | Pass | Fail |
| `vGvG_strip_05` | 0.9415 | 0.9415 | 0.9435–0.9520 | Inconclusive | Pass | Inconclusive |
| `vGvG_strip_all` | 0.9415 | 0.7990 | 0.9435–0.9520 | Fail | Pass | Fail |
| `HSX24_logN_05` | 0.8200 | 0.9415 | 0.9435–0.9515 | Fail | Pass | Fail |
| `HSX24_logN_all` | 0.8200 | 0.7835 | 0.9435–0.9515 | Fail | Pass | Fail |
| `HSX24_2logN_05` | 0.7965 | 0.9415 | 0.9420–0.9515 | Fail | Pass | Fail |
| `HSX24_2logN_all` | 0.7965 | 0.7625 | 0.9420–0.9515 | Fail | Pass | Fail |
| `HSX24_plugin_05` | 0.9430 | 0.9415 | 0.9915–0.9955 | Inconclusive | Fail | Fail |
| `HSX24_plugin_all` | 0.9430 | 0.9005 | 0.9915–0.9955 | Fail | Fail | Fail |

`MC21_t_05` directly compares the folded Student-t probability bands
using each branch's own residual degrees of freedom. Its worst null
cell was (h_a, h_b) = (2, 0), N = 1,000, without a covariate. It accepted
1,876 of 2,000 replications (0.9380, exact 95% CI [0.9265, 0.9482]).
The interval crosses the 0.94 floor, so the gate is inconclusive.
The LR analogue substitutes folded normal probabilities evaluated at
the square roots of the branch LRs. `MC21_LR_all` uses the squared normal
quantile at the lower edge of the larger branch's probability band as
the profile LR cutoff. Theorem 1 and Corollary 1 apply to the independent
normal experiment at a unit-fraction level. Theorem 3 concerns rejection
regions approaching the composite-null set. Interval inversion and the
finite-sample Student-t substitution fall outside those results.
[Miles and Chambaz, arXiv v4, §3.1](https://arxiv.org/html/2107.07575v4).

The `vGvG_strip` prototype sets the strip width so that its added rejection
probability at the independent standard-normal origin is 0.05 − 0.05²,
with the larger absolute statistic restricted to [0.25, 1.96), using the
unrounded normal quantile as the upper endpoint. The authors' nearly
similar test uses an optimized boundary. Their separate LR optimality
theorem restricts the comparison to the coherent invariant class.
[Nearly similar test, §4](https://arxiv.org/pdf/2012.11342)
[LR optimality, Theorem 1](https://www.mdpi.com/2073-8994/14/1/178).

The `HSX24_*` rules are analytical joint-significance analogues. They use
the squared 95th quantile of the smaller of two independent absolute
normal variables, one centered and the other shifted. The threshold
variants set the shift to zero until the larger square-root branch LR
exceeds √(2 log N) or 2√(log N), respectively. The plug-in variant uses that
larger value directly. The scorer interpolates the squared quantiles on
the original 0 to 8 grid with spacing 0.01. The authors' Theorem 2 treats
product-of-coefficients bootstrap consistency under the assumptions of
Theorem 1, fixed targets with local perturbations, and a diverging
threshold smaller in order than √N. Appendix B treats the resampled
joint-significance statistic. The lattice results here apply to the
analytical analogues. [He, Song, and Xu, §3.2 and Appendix B](https://arxiv.org/html/2310.11630).

No rule passed both gates. Phase 8 remains NO-GO with the χ²₁ interval
and `near_double_null` label. I7, the package release, and pmed
submission remain on hold. All scoring used saved replications.

## Sources

- Miles, C. H., & Chambaz, A. (2021). Optimal tests of the composite null hypothesis arising in mediation analysis. arXiv:2107.07575.
- van Garderen, K. J., & van Giersbergen, N. (2022). A nearly similar powerful test for mediation. arXiv:2012.11342.
- He, Y., Song, P. X.-K., & Xu, G. (2024). Adaptive bootstrap tests for composite null hypotheses in the mediation pathway analysis. *JRSS-B, 86*(2), 411–434. arXiv:2310.11630.
