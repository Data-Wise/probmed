# Plan: MBCO at the composite null, and the T22 claim (Codex review, 2026-10-14)

## Findings being addressed
- **F1 (high).** MBCO 95% intervals cover 0.997–0.998 at a = b = 0 and about 0.994 near
  (0.05, 0.05). Cause: at p* = 0.5 the profile LR is the minimum of the two branch
  statistics (a = 0 and b = 0). With roughly independent χ²₁ branches, P(reject) ≈ 0.05²,
  so coverage ≈ 0.9975, which matches what was observed. `R/compute-mbco.R:115` uses
  qchisq(level, 1) everywhere.
- **F2 (medium).** The T22 addendum says coverage is nominal at N ≥ 250, but boot-1-N1000 is
  INCONCLUSIVE (0.9345 [0.9228, 0.9449]) under the frozen rule.

## Architecture decisions
- **Measure before calibrating.** An adaptive critical value can fix the double null and break
  nearby points. That is the non-uniformity problem in Miles & Chambaz (2021),
  van Garderen & van Giersbergen and He, Song & Xu (2024). No calibration ships unless a
  lattice harness shows coverage ≥ 0.94 at every point near the origin.
- **Calibration can be evaluated offline.** Coverage of a test inversion equals acceptance of
  the truth. So the harness records, per replication, LR(true p*), LR(a = 0) and LR(b = 0).
  Any cutoff rule c(·) is then scored as P(LR_true ≤ c), with no refitting. Widths come from
  a subsample with full inversion.
- **Label first, calibrate only if it is proven safe.** The diagnostic flag and documentation
  ship regardless; the calibration is a separate, opt-in, gated task.

## Dependency graph
```
T1 (F2 report fix) ── independent
T2 harness (records LR stats) ──► T4 calibration candidates (offline) ──► T5 opt-in API (gated)
T3 diagnostic flag + docs ── independent of T2; T5 reuses its trigger
```

## Tasks

### T1: Correct the T22 claim (XS)
**Acceptance:** the CALIBRATION addendum and NEWS say boot-1-N1000 is INCONCLUSIVE, and the
supported claim excludes it; no other number changes.
**Verify:** grep finds no remaining "nominal at N ≥ 250 for every configuration".

### T2: Composite-null coverage harness (M)
`inst/sim/calibration/calib_mbco_null.R`: a local-alternative lattice
(a, b) = (h_a, h_b)/√N with h ∈ {0, 0.5, 1, 2, 4}² (25 points) × N ∈ {100, 1000},
2000 replications per point. It records LR_true, LR_a0, LR_b0 and the branch Wald z values.
Interval width comes from the first 200 replications per point with full inversion.
**Acceptance:**
- coverage at (0, 0) reproduces the calibration (about 0.997, as a positive control);
- coverage at single nulls with h_other = 4 is about 0.95;
- output is a tidy CSV with a manifest.
**Verify:** unit tests in `test-mbco-null-harness.R` on the helper that computes LR stats
(known answer at the MLE gives LR_true = 0 for p* = P̂), plus a smoke run.

### T3: Near-double-null diagnostic and labeling (S)
- `@diagnostics$pmed$near_double_null` is TRUE when both branch LR statistics (H0: a = 0 and
  H0: b = 0) are below qchisq(level, 1), i.e. neither path is individually significant.
- `print`/`summary` add one note saying the interval is conservative there.
- Docs (`pmed()` method docs, NEWS, calibration report) state the 0.997 coverage at the
  double null.

**Acceptance, TDD (tests written first and failing):**
- the flag is TRUE on simulated a = b = 0 data and FALSE when either path is strong;
- the note appears in print only when the flag is TRUE.

**Verify:** full suite plus R CMD check.

### Checkpoint A (after T1–T3)
Suite green, harness smoke run reproduces the positive controls, and a PR is opened for T1 + T3.

### T4: Offline evaluation of calibrations (M; T2 on hopper)
Candidates, scored from the T2 records:
- (a) qchisq(level, 1) (status quo);
- (b) a pretest-switched cutoff: qchisq(1 − √α, 1) when both branch statistics fall below a
  threshold λ_N, otherwise qchisq(level, 1), with λ_N ∈ {crit, log N} following He–Song–Xu
  thresholding;
- (c) a smooth interpolation between (a) and (b) on max(LR_a0, LR_b0).

**Prespecified acceptance:** min coverage over the lattice ≥ 0.94 at both N, coverage at
(0, 0) ≤ 0.97, and mean width reduction reported.

### Checkpoint B
The user decides from the T4 table whether any candidate ships.

### T5: Opt-in calibrated MBCO (M, only if T4 passes)
A `null_calibration = c("none", "adaptive")` argument, default "none", with tests that
reproduce the T4 coverage on a fixed seed. The calibration report gets an addendum.

## Risks
| Risk | Mitigation |
|---|---|
| A calibration undercovers near (not at) the null | Lattice in √N-local coordinates; acceptance uses the minimum over the lattice |
| The branch statistics are not independent with covariates | T2 also runs with one covariate |
| Harness cost | 3 profile evaluations per replication (no inversion), widths on a subsample; local pilot, hopper for the full run |

## Sources
- Miles & Chambaz, optimal tests of the composite null in mediation: https://arxiv.org/pdf/2107.07575
- van Garderen & van Giersbergen, nearly similar test for mediation: https://arxiv.org/pdf/2012.11342
- He, Song & Xu (JRSS-B 2024), adaptive bootstrap for composite nulls: https://arxiv.org/pdf/2310.11630
