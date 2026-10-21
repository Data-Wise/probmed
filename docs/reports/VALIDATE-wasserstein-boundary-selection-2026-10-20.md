# Boundary point estimation and Fieller calibration: fresh validation

**Decision: neither candidate passes. Keep PRs #54 and #55 DRAFT and blocked.**
No production adoption, diagnostic relaxation, merge, required-status switch,
release, or submission. Full exact-commit local software CI remains pending.

## Prespecification and provenance

Design: [frozen plan](PLAN-wasserstein-boundary-selection-2026-10-20.md), commit
3303ccb. Source frozen at f27e7960114f9e8fb01e710b6c39aa03fc6bd0a6 before any fresh
datasets. Seven unchanged cells, 2000 independent datasets each; seed
20561020 + cell*10000 + replicate. Independent audit confirms no overlap with
the three preceding studies. 14000 datasets, 44000 paired method rows, zero
failed fits; elapsed 317.84 seconds, R 4.6.1.

The numerical precision check initially failed: 32768 common Gaussian draws
gave SE/reference ratio 1.102004 on a rare-event fixture, beyond the frozen 10%
tolerance. The documented pre-simulation amendment doubled draws to 65536;
the worst ratio became 1.038229 against an independent 262144-draw reference.
Thresholds, estimators, tolerances, seeds and diagnostics were not retuned.

Raw [replicates](data/wasserstein-selection-2026-10-20/replicates.csv),
[summary](data/wasserstein-selection-2026-10-20/summary.csv), and
[source hashes/session](data/wasserstein-selection-2026-10-20/provenance.txt).
Replicate SHA256:
`1b4e271e1c1bc69185d9da31a393e0a5c3e701c0381d3caff4ba0c9790bece21`.

## Results: interval calibration is not point-estimation calibration

Frozen gates remain: zero failures; absolute bias/empirical SD < .10;
mean SE/empirical SD in [.90,1.10]; coverage in [.92,.98]. No rounding at gates.

All three Fieller methods share the same intervals. Their observed coverage
passes in every binary cell. Binomial MCSE is .00372–.00514; this is empirical
finite-sample evidence, not a theorem of uniform or exact coverage.

| Binary cell | Fieller coverage | k4 bias/SD; SE/SD | Select95 bias/SD; SE/SD | Select995 bias/SD; SE/SD |
|---|---:|---:|---:|---:|
| Regular n=500 | .9475 | .0636; 1.0658 | .4423; .7355 | .8929; .4945 |
| Regular n=1000 | .9500 | .0267; 1.0338 | .0779; 1.2818 | .2974; .8351 |
| Zero direct, P=1 | .9440 | 3.2082; 2.8478 | .2430; 1.0680 | .0761; 1.9622 |
| Zero indirect, P=0 | .9450 | 2.6219; 2.6005 | .4601; 1.0860 | .7488; .9649 |
| Weak direct | .9715 | 1.0407; 2.7565 | 1.8730; 1.1012 | 7.8773; 2.5141 |

The k4/Fieller combination passes both regular binary cells and fails all three
stress cells. Both Gaussian raw/Wald controls pass (coverage .9555/.9595).
Neither selection candidate passes any binary cell's complete diagnostics.
Raw/Wald binary n=500 again fails bias/SD (.1192), and raw boundary coverage
is .4815/.5545/.8055. All differences above are paired on the same datasets.

## Why selection failed

At regular binary n=500, primary selection snaps to an axis in 22.95% of fits;
the stronger selector does so in 52.70%. This induces bias where the true score
is interior. At weak direct, selected rates are 94.95% and 99.65%; the stronger
candidate's small empirical SD (.02872) around a badly displaced estimate makes
its bias/SD particularly poor. This is not an interval-coverage failure.

For zero indirect, both contrasts are classified weak in 15.45% and 41.95% of
fits. The prespecified fallback avoids an arbitrary double-zero score but
retains substantial positive endpoint bias. For zero direct, stronger selection
does reduce bias/SD below .10, but its Gaussian-proxy SE is almost twice the
empirical SD: passing one diagnostic cannot substitute for the full contract.

The primary selector's estimates lie inside all its intervals. Strong selection
is incompatible with the 95% inversion on some datasets: point-in-interval
rates range from .7025 to .9670 across binary cells. We did not enlarge intervals
after seeing this. Its more aggressive threshold is therefore not a drop-in
replacement for the point estimate accompanying this interval.

## What is supported, and what is not

The union-of-lines chi-square-1 inversion has an independent statistical
justification: under a regular joint-contrast CLT, its rejection event is a
subset of rejection for the true line. It is not the previous 95% joint ellipse
with a cosmetic label change. Fresh coverage supports this limited interval
candidate; fitted contrast covariance still gives asymptotic, not exact normal
inference. The double-zero population ratio remains undefined and is not a
coverage target. Nearby weak-signal uniform validity is not established.

Hard thresholding and the constrained Gaussian-experiment SE are rejected as a
complete repair under these frozen cells. These results do not rule out every
boundary estimator, nor validate model-selected proxy uncertainty generally.

Sources read for design: von Luxburg and Franz (2009), sections 2–3,
https://www3.stat.sinica.edu.tw/sstest/oldpdf/A19n312.pdf (ratio-test geometry);
Pötscher and Schneider (2010), sections 1–3,
https://arxiv.org/pdf/0806.1652 (thresholding/nonuniformity cautions). Neither
source proves validity for this fitted nonlinear ratio.

## Verification and next research boundary

`Rscript scripts/test-wasserstein-boundary-selection.R` passes scalar/vectorized
agreement, axis selection, constrained-center orthogonality, independent proxy
precision, scale invariance, and independent linear-null endpoint equations.
`Rscript scripts/audit-wasserstein-boundary-selection.R` independently reproduces
all 22 summaries, seed disjointness, Fieller inclusion for every binary row,
and selected-axis/center calculations. First 10 fits per binary cell also agree
with production point/SE oracles within 1e-9/1e-8. These are research checks,
not a passing full-package CI claim.

Next: derive a boundary point-estimation risk/bias trade-off before proposing a
new smooth constrained estimator and its uncertainty. Preserve all regular and
weak-signal cells; do not escalate thresholds to chase exact-null bias. A new
candidate needs a frozen design, independent numerical tests, and fresh seeds.
Both PRs remain blocked until the unchanged complete diagnostics pass.
