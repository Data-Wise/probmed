# I6 width-failure investigation (2026-10-20)

## Decision boundary

I6 remains **FAIL** under its frozen all-cell gate. This work diagnoses the
Gaussian fresh-draw candidate `d` at `K = 200,000`; it neither changes that
gate nor authorizes I7, package-method claims, pmed edits, or release. The
existing [s6 diagnostic](../docs/reports/I6-S6-N1000-DIAGNOSTIC-2026-10-19.md)
already found a finite-`K` explanation. This plan makes that explanation
reproducible across the retained Gaussian grid and tests its limits.

## ADEMP

- **Aim:** Determine whether the s6/`N = 1000` width excess is consistent with
  finite-`K` score noise, and look for pairing, failure, or endpoint-asymmetry
  alternatives in the saved I6 outputs.
- **DGP:** Reuse the frozen I6 `s1,s6,s7,s8` Gaussian cells at
  `N = 50,100,250,500,1000`, 2,000 outer replications per cell. No new data
  generation or seed change.
- **Estimand:** The same tie-adjusted Gaussian `P_med` and closed-form truths
  as the I6 gate; this is a method-comparison diagnostic, not a new estimand.
- **Methods:** Pair candidate `d = Binomial(200000,p)/200000` with in-run exact
  `e` on each saved dataset and bootstrap-index stream. Do not pair I6 to I2
  or the separate pmed run.
- **Performance:** Report mean-width ratio, paired outer-replication 95% CI,
  a finite-`K` variance prediction, endpoint shifts, coverage-flip counts,
  and all interval/failure/key checks. The diagnostic CI is descriptive, not
  the prespecified coverage-equivalence test.

## Implementation and checks

1. Read the 300 retained `i6_*.csv` parts into an isolated scratch directory.
   Check 60,000 unique `(scenario,N,rep)` keys, 2,000 rows in each of 30 cells,
   `boot_attempted = 2000`, finite ordered Gaussian intervals, and zero failed
   `d` resamples. Preserve the source and runner identities in the report.
2. Run `inst/sim/i6_width_investigation.R` with a fixed diagnostic bootstrap
   seed. For each of the 20 Gaussian cells, predict
   `sqrt(1 + mean[p(1-p)] / (K * (mean_width_e / (2 z_0.975))^2))`, using the
   saved exact bootstrap medians as `p` proxies. Compare to observed `d/e`.
   The normal-spread approximation is explanatory, not an endpoint proof.
3. Test the diagnostic code with small known-answer paired fixtures, including
   duplicate keys, failed draws, and a one-sided endpoint shift. Reconcile
   s6/`N = 1000` against the frozen I6 report before interpreting other cells.
4. Report where prediction and observation differ and whether the saved rows
   suggest a competing mechanism. Do not select a new `K` from these results.

**Stop/next decision:** Preserve I6 FAIL regardless of descriptive agreement.
If a different Gaussian method is proposed, review its implementation and
predeclare a fresh validation plan before any run. The binary report-only
limitations remain separate.
