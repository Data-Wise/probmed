# Plan: Monte Carlo error in simulation-based bootstrap intervals for P_med

Created 2026-10-14 after the pmed K_BOOT = 50k rerun failed its null-coverage gate and the
Codex adversarial review (verdict: needs-attention). Covers the pmed manuscript simulation and
probmed issue #47; one method decision serves both.

## What is established, and how firmly

| Claim | Evidence | Strength |
|---|---|---|
| Fresh per-resample draws inflate bootstrap width and coverage | probmed T20 (exact vs MC: 0.977 → 0.945); pmed width ratios predicted 1.10/1.15, observed 1.10/1.145 | Verified (positive control) |
| The manuscript's MBCO width advantage is that artifact | K = 5k and K = 50k reruns: bootstrap/MBCO width 0.97–1.01 in all Gaussian cells | Verified, two runs |
| The percentile bootstrap is valid at the null | Exact per-resample P_med, s6, N = 1000: 0.951 [0.941, 0.960] | **One cell only** |
| CRN undercovers at the null from slope noise | K = 5k 0.916 / 0.933, K = 50k 0.9315; mechanism derived, not measured | **Hypothesis** |
| CRN is fine away from the null | probmed T22 (all non-null) PASS; pmed s1/s2/s4/s7/s8 ≈ 0.94–0.95 | Probable; s8 and s2 slightly low |

## Open questions

1. Is the exact-bootstrap validity at the null general (all N, both null branches s6 and s8), or specific to one cell?
2. Is the CRN failure the slope noise as derived? Its prediction can fail: per-dataset width sd (CRN vs exact) should scale like 1/√(K/√N).
3. Which Monte Carlo control is valid everywhere, including binary outcomes with no closed form?
4. Are the reported coverages computed over all attempted replications? `analyze.R` uses `na.rm = TRUE`, which drops failures instead of counting them as misses.
5. Does MBCO itself have issues the bootstrap comparison hid (s5 bias 0.004 → 0.92 at N = 1000)?

## Tasks

### I1: Audit the reported numbers from raw replications (S, local, ~15 min)
- [x] For the fresh, crn5k and crn50k runs: per cell, attempted vs finite replications, and coverage with failures counted as misses
- [x] Flag any cell where the two coverage definitions differ by more than 0.005
**Can fail:** yes, if failures are being silently dropped. **Output:** a table in the report. **Done 2026-10-15:** 0 drops in 210 cells; see `docs/reports/I1-COVERAGE-AUDIT-2026-10-15.md`.

### I2: Exact-bootstrap reference grid for the Gaussian mechanisms (M, local or hopper, ~2 h)
- [ ] Closed-form P_med per resample for s1, s2, s6, s7 and s8 × all 5 N, 2000 reps, same datasets as the crn50k run (reuse its seeds)
- [ ] Also exact MBCO already on file; compare
- [ ] Positive control: exact must reproduce 0.951 at s6 N = 1000
**Answers Q1.** This is the ground truth the candidate controls are judged against.

### I3: Test the slope-noise mechanism (S, local, ~1 h)
- [ ] s6 and s1 at N ∈ {250, 1000}, K ∈ {5k, 20k, 50k, 200k}, 500 reps: the sd of (CRN width / exact width) per dataset
- [ ] Prediction: relative sd ≈ sd(W) / (E[W] √n_flip), so it grows with N and shrinks like K^{-1/2}; much smaller at s1 than s6
- [ ] If the scaling doesn't match, the mechanism is wrong; record that and rely on I4's empirical comparison
**Answers Q2.**

### I4: Head-to-head comparison of Monte Carlo controls (M–L, hopper, ~3–4 h)
Same datasets, same resamples, for each candidate:
- (a) fresh draws (original)
- (b) CRN + control variate, K = 50k (current)
- (c) **fresh draws + variance correction**: shrink the deviations by √(1 − σ̂²_MC / var(P*)), with σ̂²_MC = var(δ_k)/K estimated from the draws themselves, so it also works for binary outcomes with ties
- (d) fresh draws, large K (e.g. 200k), brute force
- (e) exact, for the Gaussian cells (reference from I2)

Cells: s1, s3, s4, s6, s7, s8 × N ∈ {50, 250, 1000}, 2000 reps.
- [ ] Coverage over all attempted reps, mean width, width sd, and the per-dataset agreement with (e)
- [ ] For binary s3/s4 (no exact): use (d) at K = 1e6 as the reference
- [ ] Check (c) where the bootstrap distribution is skewed (s3 small N, s7 near the boundary); the correction assumes additive independent noise and scales the quantiles linearly
**Answers Q3.** Prespecify the decision rule before running: choose the simplest control whose coverage CI overlaps (e) in every Gaussian cell and whose widths are within 2% of (e).

### I5: MBCO side check (S, local, ~30 min)
- [ ] s5 at N = 1000: confirm the 0.92 is the 0.004 misspecification bias (bias/SE vs coverage shortfall), for both methods
- [ ] Note in the manuscript as a misspecification result, not an MC one
**Answers Q5.**

### I6: Authoritative rerun with the chosen control (M, hopper, ~2 h)
- [ ] Full 400-task grid with the I4 winner; manifest committed before submission (control, K, seed base 20261017, gate)
- [ ] Gate: every Gaussian cell's coverage CI overlaps the I2 exact reference; widths within 2% of exact
- [ ] Coverage computed over all attempted replications (I1)

### I7: Manuscript and package updates (M, after I6 passes)
- [ ] pmed.qmd + supplement: rewrite the MBCO-versus-bootstrap comparison (no width advantage), the large-N coverage text, line 432's MC statement, and the method description
- [ ] probmed #47: apply the same control to `model_bootstrap` if (b) loses to (c), and add null cells to the calibration
- [ ] Push pmed `draft` only with the user's OK

## Checkpoints
- **After I1–I3 (one evening, local):** review with the user. If I3 refutes the mechanism, reconsider I4's candidate list.
- **After I4:** the user picks the control from the comparison table before I6.
- **After I6:** gate passes → I7; gate fails → stop and report.

## Risks
| Risk | Mitigation |
|---|---|
| Variance correction distorts skewed bootstrap distributions | I4 checks skewed cells; fall back to (d) there |
| Seeds differ across runs, so cell comparisons mix MC and sampling noise | I2 and I4 reuse one dataset set per cell |
| hopper queue (16 concurrent) stretches I4/I6 | Pilot via sbatch first; tasks are ~7 min each |
| Deciding the method after seeing the results | Prespecified decision rule in I4 and gate in I6 |
