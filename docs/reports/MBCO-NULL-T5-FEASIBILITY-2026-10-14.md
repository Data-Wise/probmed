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

## Sources

- Miles, C. H., & Chambaz, A. (2021). Optimal tests of the composite null hypothesis arising in mediation analysis. arXiv:2107.07575.
- van Garderen, K. J., & van Giersbergen, N. (2022). A nearly similar powerful test for mediation. arXiv:2012.11342.
- He, Y., Song, P. X.-K., & Xu, G. (2024). Adaptive bootstrap tests for composite null hypotheses in the mediation pathway analysis. *JRSS-B, 86*(2), 411–434. arXiv:2310.11630.
