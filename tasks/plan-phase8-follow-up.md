# Phase 8 research plan for the MBCO interval

**Research started 2026-11-25. Phase 8 remains NO-GO.** Retain the χ²₁ interval
and `near_double_null`. T5 is closed. I7, package publishing, and pmed
submission remain on hold.

The goal is a less conservative 95% confidence interval for `P_med`
in the Gaussian single-mediator model near the double null. On the saved
lattice, the current rule includes the truth in 99.60–99.85% of origin
replications. No earlier prototype cleared both gates.
`MC21_t_05` was inconclusive and `MC21_t_all` is excluded from interval
evidence. Each candidate interval needs a test at every target value.
[Saved audit](../docs/reports/MBCO-NULL-T5-FEASIBILITY-2026-10-14.md#offline-prototype-scores-2026-11-23)

## 1. Verify the Gaussian benchmark

- [x] Archive the sharp 95% proof and exact-arithmetic certificate.
  Checked the large-mean geometry, compact-region certificate, and
  strong-path limit. Both certificate implementations reproduced the
  bound from a clean directory.
  [Proof and replay](../docs/reports/phase8-gaussian-benchmark/README.md)
- [x] Record the benchmark's scope. In the independent, unit-variance
  Gaussian shift experiment, the proof gives a worst-case rejection
  probability of 0.05 for the χ²₁ cutoff at every product target. Its
  inversion is the product image of a disk. At the origin it covers
  99.75%.

The benchmark is ready for the test-family comparison. Check the exact
95% statement against the literature before any novelty claim.

## 2. Develop a test family and check its inversion

- [x] Extend the fixed-target rule to the axes and both product signs.
  Prove the uniform 0.05 size bound for every target. Its inversion is
  a single interval when nonempty, with exact origin empty-set
  probability equal to its radial rejection budget.
- [x] Compare the rule to χ²₁ and χ²₂ by deterministic integration and
  derive its strong-path limit. Reaching 97% coverage at the ideal
  origin costs at least 2.978% empty inversions there and 24.41% wider
  intervals in the strong-path limit. [Proof and checks](../docs/reports/phase8-gaussian-benchmark/test-family.md)
- [ ] Investigate whether the χ²₁ rule's local size slack can reduce
  origin coverage while keeping the inversion nonempty and preserving
  χ²₁ strong-path widths. Establish a uniform size bound and verify the
  inversion before moving a revised candidate to fitted-model work.

Advance when there is a full test family with a uniform size argument
and a practically useful gain over χ²₁. Define the gain before screening.
A better test of `P_med = 0.5` can be evaluated separately under
its own size and power criteria. It supplies a decision for one target.
[Miles and Chambaz](https://arxiv.org/html/2107.07575v4)

## 3. Connect the rule to the fitted mediation model

The [fitted-model note](../docs/reports/phase8-gaussian-benchmark/fitted-model.md)
proves the local bound for the exact χ²₁ profile LR on bounded local
path regions, under Gaussian regressions with common controls.
Stronger-path sequences and numerical inversion are open. Any adaptive
test family will need its own transfer argument.

- [ ] Derive the local `P_med` constraint with estimated
  residual variances and standardized path estimates. Compare the
  constrained and unconstrained profile likelihoods with their Gaussian
  approximations. State the uniformity domain and address both paths
  weak, either path stronger, and covariate adjustment.
- [ ] Check inversion against both constrained likelihood branches, using
  deterministic fixtures. Record the treatment of fitting failures and
  endpoint accuracy. For a fitted-data claim, state the assumptions and
  whether the guarantee is asymptotic or finite sample.

Advance to an interval evaluation only when the regression argument and
numerical verification support the stated scope. If the argument stops
at the ideal experiment, retain the result as theory and keep NO-GO.

## 4. Screen the candidate and make a decision

- [ ] Freeze one candidate and its width comparison before scoring. The
  saved 100-cell lattice supports truth-inclusion scoring for rules that
  use the recorded likelihood ratios and metadata. It contains baseline
  endpoints for 50 replications per cell. If the rule needs extra inputs
  such as signed path estimates, or candidate widths need extra data,
  write a minimal evidence proposal for author approval.
- [ ] Apply the frozen T5 gates. Require a pointwise 95% Clopper–Pearson
  lower bound of at least 0.94 in each of the 100 cells and an upper bound
  of at most 0.97 in each of the four origin cells. Compare paired widths
  and report inversion failures. A bound that crosses its threshold
  leaves the gate inconclusive. These are lattice screening criteria.
  [Scoring specification](../docs/specs/SPEC-2026-10-15-t5-prototype-audit.md)
- [ ] Bring the proof, implementation checks, coverage, and width results
  to the author for a new decision. Any failed or inconclusive gate keeps
  Phase 8 at NO-GO. Start with proof review and deterministic checks. New
  simulations, a package API change, and manuscript work each require a
  separate decision.
