# Plan: Monte Carlo error in simulation-based bootstrap intervals for P_med

Created 2026-10-14 after the pmed K_BOOT = 50k rerun failed its null-coverage gate and the
Codex adversarial review (verdict: needs-attention). Covers the pmed manuscript simulation and
probmed issue #47; one method decision serves both.

## What is established, and how firmly

| Claim | Evidence | Strength |
|---|---|---|
| Fresh per-resample draws inflate bootstrap width and coverage | probmed T20 (exact vs MC: 0.977 → 0.945); pmed width ratios predicted 1.10/1.15, observed 1.10/1.145 | Verified (positive control) |
| The manuscript's MBCO width advantage is that artifact | K = 5k and K = 50k reruns: bootstrap/MBCO width 0.97–1.01 in all Gaussian cells | Verified, two runs |
| Exact Gaussian percentile-bootstrap coverage met the historical pmed gate in six null cells at N >= 250 | I2/I6 references; completed pmed 6,000-per-cell run reported 94.38%–95.32% and six two-sided 90% exact-binomial intervals inside its frozen 93%–97% band; s6/s7 truths are 0.5 | Verified for that reported Gaussian grid and historical gate; not a prospective pass of the later 94% minimum, binary support, or general model support |
| Finite-K Monte Carlo noise broadens intervals as sampling uncertainty shrinks | I3 scaling check; I4 and I6 s6/N=1000 misses; I6 observed d/e width 1.029588 vs approximation 1.029691 | Verified in the tested cells |
| CRN may be adequate away from the null | probmed T22 (non-null cells) PASS; pmed non-null s1/s2/s4/s8 ≈ 0.94–0.95 | Probable in tested non-null cells; s8 and s2 slightly low; s7 is a null branch |

## Open questions

1. Does the exact bootstrap meet a separately prespecified 94% minimum in all six s6/s7 null cells? The completed pmed run used an earlier 93%–97% rule; s8 is not a null branch.
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
- [x] Closed-form P_med per resample for s1, s2, s6, s7 and s8 × all 5 N, 2000 reps, same datasets/resample RNG stream as crn50k (seed base 20261016; 8-worker chunks)
- [x] Compare exact bootstrap with saved CRN/MBCO: width ratios exact/CRN 0.9988–1.0003, exact/MBCO 0.9658–0.9981; MBCO coverage 0.9400–0.9815, max paired coverage difference 0.0125
- [x] Positive control: s6 N = 1000 coverage 0.944 [exact-binomial 95% CI 0.9330, 0.9537], containing 0.951; no equality requirement.
**Partial evidence for Q1.** At 2,000 replications per cell, the tested s6/s7 null cells from N = 250 onward were near 95%; s6 overcovers at N = 50 and 100. This predates the later 94% decision and is not its formal confirmation. See `docs/reports/I2-EXACT-REFERENCE-2026-10-15.md` and the 50,000-row result grid. **Done 2026-10-15.**

### I3: Test the slope-noise mechanism (S, local, ~1 h)
- [x] Produce an exact comparator subset for s6 and s1 at N ∈ {250, 1000} (500 matched datasets/seeds), then measure the sd of (CRN width / exact width) per dataset for K ∈ {5k, 20k, 50k, 200k}; this lets I3 run before the full I2 grid.
- [x] Prediction: relative sd ≈ sd(W) / (E[W] √n_flip), so it grows with N and shrinks like K^{-1/2}; much smaller at s1 than s6
- [x] The scaling matched the prediction in these four cells; historical array rows were not reproduced because its per-worker RNG streams were not replayed.
**Answers Q2.**
**Done 2026-10-15:** 500 datasets per cell and 1,000 bootstrap resamples per dataset; fresh internally paired run. The user directed “do next” at the I1–I3 checkpoint on 2026-10-15; I2 proceeded.

- [x] After I1–I3 checkpoint: user directed proceeding to I2 (2026-10-15).

### I4: Head-to-head comparison of Monte Carlo controls (M–L, hopper) — complete 2026-10-19
Same datasets, same resamples, for each candidate:
- (a) fresh draws (original)
- (b) CRN + control variate, K = 50k (current)
- (c) **fresh draws + variance correction**: shrink the deviations by √(1 − σ̂²_MC / var(P*)), with σ̂²_MC = var(δ_k)/K estimated from the draws themselves, so it also works for binary outcomes with ties
- (d) fresh draws, large K (e.g. 200k), brute force
- (e) exact, for the Gaussian cells (reference from I2)

Cells: s1, s3, s4, s6, s7, s8 × N ∈ {50, 250, 1000}, 2000 reps.
- [x] Coverage over all attempted reps, mean width, width sd, and the per-dataset agreement with (e)
- [x] For binary s3/s4 (no exact), use (d) at K = 1e6 as the reference
- [x] Check (c) where the bootstrap distribution is skewed (s3 small N, s7 near the boundary)
**Answers Q3.** Prespecification approved by the user on 2026-10-15; do not inspect I4 comparison results before the manifest and pilot are recorded.

- Coverage equivalence margin: ±0.01 absolute coverage (±1 percentage point) for each matched cell. For each candidate/reference pair, compute the replication-level coverage-indicator difference and form a paired 90% bootstrap CI by resampling outer replication IDs (10,000 bootstrap resamples, preserving pairing). This is a two-one-sided-test equivalence rule at α = .05; pass only when the full CI lies inside [−0.01, 0.01]. Apply against I2 exact in Gaussian cells and against fresh draws at K = 1e6 in binary s3/s4 cells.
- Width gate: candidate/reference mean interval-width ratio must be in [0.98, 1.02] in every cell. Use the same references as above.
- Rationale: ±1 percentage point caps the allowed coverage difference from the reference; the 2% width tolerance is already specified in this plan. Requiring every cell to pass avoids averaging away a weak mechanism or sample-size condition.
- Selection: among candidates passing both gates in every applicable cell, choose the simplest. If none passes, report no-go and retain the current method. I6 remains an independent-seed validation against the I2 exact Gaussian reference. The first full 200-replication chunk reproduced all 200 I2 Gaussian point estimates and interval endpoints exactly (max absolute difference 0); the one-replication smoke test was invalid for RNG replay because it changed the `mclapply` worker topology.

**Result:** Hopper job 4339014 completed 180/180 tasks with no missing active intervals or failed resamples. No candidate passed every applicable cell: b and c passed 16/18, d passed 11/12 Gaussian cells, and a passed 2/18. The frozen rule therefore returned no-go and retained current b. See `docs/reports/I4-MC-CONTROLS-2026-10-19.md`.

### I5: MBCO side check (S, local, ~30 min)
- [ ] s5 at N = 1000: confirm the 0.92 is the 0.004 misspecification bias (bias/SE vs coverage shortfall), for both methods
- [ ] Note in the manuscript as a misspecification result, not an MC one
**Answers Q5.**

### I6: Authoritative rerun with the chosen control (M, hopper) — complete; FAIL 2026-10-19
I4 ended in a no-go, so I6 did not run "the I4 winner". The user selected fresh K=200,000 candidate d for Gaussian cells and current b for report-only binary cells. The committed manifest is `inst/sim/results/calibration/i6_authoritative_manifest_2026-10-19.csv`.
- [x] Decision (a) 2026-10-19: use d (fresh K=200k) for Gaussian cells. Alternatives were: keep b and accept a documented large-N null undercoverage (mechanism in CRN50K-RERUN), or use d (K=200k) for Gaussian cells, or make Gaussian P_med exact. Record the evidence that would justify proceeding with b
- [x] Decision (b) 2026-10-19: binary cells report-only, run with b and no gate
- [x] Reference is in-run exact candidate e on the same datasets and bootstrap indices; seed base 20261017 is independent of I2/I4
- [x] Apply the unchanged Gaussian gate: paired 90% coverage-difference CI strictly inside (-0.01, +0.01) and width ratio 0.98-1.02 in every cell
- [x] Count coverage over all attempted replications using closed-form truths; s6/s7 are exactly 0.5
- [x] Record within-bootstrap failed resamples separately; all active `boot_failed` counts were zero
- [x] Pilot the committed `inst/sim/hopper/i6/runner.R`, then run Hopper job 4341753: 300/300 tasks completed and 60,000 rows validated

**Result:** d passed both gates in 19/20 Gaussian cells. It failed s6/N=1000: coverage 0.9515 vs exact 0.9425, paired difference +0.0090 with 90% CI [0.0055, 0.0125], and width ratio 1.029588. Under the frozen all-cell rule, I6 fails and work stops before I7. Binary s3/s4 results are descriptive only; s3 undercoverage at low N must be disclosed if reused. See `docs/reports/I6-GATE-2026-10-19.md` and `docs/reports/I6-S6-N1000-DIAGNOSTIC-2026-10-19.md`.

### Cross-project exact-Gaussian follow-up — historical result and optional new decision

- [x] pmed reports a separate 6,000-replication-per-cell exact-Gaussian run with seed base 20261021, accepted under its previously frozen 93%–97% six-cell rule; see its `docs/reports/EDIT-REPORT-pmed-exact-gaussian-2026-10-19.md`. This does **not** reverse I6.
- [ ] Reconcile its manifest, six cell counts/bounds, hashes, and saved chunks before reusing its results beyond the scoped pmed statement. At least one reported cell has 94.38% observed coverage, below the later rule's 94.52% threshold; even a post hoc application cannot pass all six. Do not recast this as a failed *prespecified* 94% test.
- [ ] If a new prospective minimum-94% claim is needed, obtain separate compute authorization and freeze a distinct seed family, manifest, all-attempt analysis, and stop rules before a new run. See `docs/PROPOSAL-exact-gaussian-followup-2026-10-19.md`; 6,000 per cell may legitimately yield an inconclusive result.
- [ ] Keep medfit's generic fitting/extraction/bootstrap contract, RMediation's path-product CI/MBCO, and probmed's tie-adjusted P_med inference distinct. Exact-kernel package integration needs its own decision and tests; neither dependency is a P_med oracle.

### I7: Manuscript and package updates (M, only after I6 passes) — stopped; I6 failed
- [ ] pmed.qmd + supplement: rewrite the MBCO-versus-bootstrap comparison (no width advantage), the large-N coverage text, line 432's MC statement, and the method description
- [ ] probmed #47: apply the same control to `model_bootstrap` if (b) loses to (c), and add null cells to the calibration
- [ ] Push pmed `draft` only with the user's OK

## Checkpoints
- **After I1–I3 (one evening, local):** review with the user. If I3 refutes the mechanism, reconsider I4's candidate list.
- **Before I4:** agree on the paired coverage-equivalence margin and rationale before inspecting comparison results.
- **After I4:** the user picks the control from the comparison table before I6.
- **After I6:** gate passes → I7; gate fails → stop and report.

## Risks
| Risk | Mitigation |
|---|---|
| Variance correction distorts skewed bootstrap distributions | I4 checks skewed cells; fall back to (d) there |
| Seeds differ across runs, so cell comparisons mix MC and sampling noise | I2 and I4 reuse one dataset set per cell |
| hopper queue (16 concurrent) stretches I4/I6 | Pilot via sbatch first; tasks are ~7 min each |
| Deciding the method after seeing the results | Prespecified decision rule in I4 and gate in I6 |
