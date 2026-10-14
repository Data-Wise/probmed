# pmed manuscript: coverage claims checked against the probmed calibration (2026-10-14)

This was a read-only review of `~/projects/research/pmed` at `a2c0b74` (`pmed.qmd`,
`code/simulation_functions.R`, `code/config.R`), and nothing in pmed was edited.
The calibration it is checked against is `docs/reports/CALIBRATION-2026-10-14.md`
(T17, T20, T22).

## Main finding: the manuscript's percentile-bootstrap intervals carry Monte Carlo inflation

`pmed_point()` (`code/simulation_functions.R:27`) draws **K = 5,000 independent
cross-world pairs** per estimate (`code/config.R:12`), with fresh draws in every
one of the 1,000 bootstrap resamples. This is the same defect T20 confirmed in
probmed. The Monte Carlo sd is √(p(1−p)/K), about 0.007 at P_med = 0.55. Each
resample adds that variance to the bootstrap spread, so the intervals are too
wide, and more so as N grows.

The predicted inflation matches the manuscript's own numbers. The prediction
takes the MBCO width (closed-form, no Monte Carlo) as the noise-free reference:
SE = width / 3.92 and inflation = √(SE² + MC²) / SE.

| Cell | MC sd | Predicted bootstrap/MBCO width ratio | Manuscript ratio | Predicted coverage | Manuscript coverage |
|---|---|---|---|---|---|
| Mech. 1, N = 250 | 0.0070 | 1.10 | 0.066 / 0.060 = **1.10** | 0.969 | **0.97** |
| Mech. 7, N = 1,000 | 0.0042 | 1.15 | 0.0332 / 0.0290 = **1.145** | 0.976 | 0.95–0.98 |

## Claims affected

| Line (`pmed.qmd`) | Claim | Status |
|---|---|---|
| 432 | Monte Carlo error is "of order O_p((nK)^{-1/2}), small relative to the n^{-1/2} sampling error" | **Incorrect for the code as run.** K pairs are drawn in total, not per observation, so the MC error is O(K^{-1/2}), about 0.007. That is comparable to the sampling SE at N ≥ 500 (SE ≈ 0.0074 for Mech. 7 at N = 1,000). |
| 227, 542 | MBCO intervals are "consistently tighter than the percentile bootstrap… an advantage that widens in larger samples"; MBCO 13% narrower near the boundary at N = 1,000 | **The width advantage is explained entirely by the bootstrap's MC inflation** (predicted 1.10 and 1.15, observed 1.10 and 1.145). The comparison does not show that MBCO is more efficient. |
| 491, tbl-performance | Bootstrap coverage 0.94–0.97 (non-null) at N = 250; Mech. 1 and 5 at 0.97 | Values are accurate, but the reported 0.97s are MC overcoverage (predicted 0.969), not property of the method. |
| 544 | Rare binary (Mech. 3) coverage 0.97 at N = 500 and 0.99 at N = 1,000; Mech. 4 "becomes conservative at the larger sample sizes" | The large-N conservativeness matches MC inflation; the manuscript gives no explanation for it. |
| 491, 542 | Bootstrap is conservative at the null (0.99) | Probably partly MC inflation and partly the null boundary. This can't be separated without a rerun. |
| 544 | Rare-binary undercoverage at N = 50 (0.77), linked to separation | **Consistent** with the probmed calibration: separation and fit failures dominate rare outcomes at small N. |
| 542 | MBCO coverage close to nominal at both single-null branches | **Consistent.** MBCO is closed-form here, so the MC issue doesn't affect it, and probmed's MBCO passes at single nulls. The double null a = b = 0 (probmed: conservative, about 0.997) is not in the manuscript's design. |
| 469 | "with the parametric bootstrap, the returned point estimate is the mean of the bootstrap estimates" | **Still accurate** for `parametric_bootstrap` (`R/compute-bootstrap.R:75`). It does not describe `model_bootstrap`, which returns the plug-in. The paper doesn't mention `model_bootstrap`. |

## Options for the manuscript (author's decision)

1. **Rerun the bootstrap half of the simulation, recommended.** Use common random
   numbers with control-variate centering, as in probmed PR #46, or a much larger
   K for the point estimate. The cost is about that of the original run. Then
   restate the MBCO-versus-bootstrap comparison from the new widths. The likely
   outcome is that both intervals are near nominal at similar widths, so the
   width-advantage claim goes.
2. **Keep the results but reframe.** Correct line 432's error order, attribute
   the bootstrap's large-N conservativeness and width to K, and drop the MBCO
   width-advantage claims at lines 227 and 542.
3. Increase K to about 100,000 in the bootstrap resamples (MC sd 0.0016). This
   is simpler than common random numbers but costs about 20× the compute.

Option 1 or 2 is needed before submission, because line 542's "consistently
tighter" claim is currently an artifact.
