<!-- Draft comment for Data-Wise/probmed#38. Not posted; posting needs the author's yes. -->

**Gauge Table 1 post-fix rerun: done.** Collated on hopper 2026-08-30 (jobs 4311422 /
4311424, 64/64 parts, nsim = 2000 per cell); committed as
`inst/sim/results/gauge_boot_coverage_nsim2000_postfix.csv` beside the pre-fix file.

With the corrected corner-EIF weights (PR #34) every cell is calibrated:

| Quantity (8 cells) | Pre-fix | Post-fix |
|---|---|---|
| Wald coverage of W | 0.859–0.913 | 0.942–0.961 |
| Percentile coverage of W | 1.000 | 0.940–0.953 |
| Wald coverage of P_med | 0.899–0.934 | 0.943–0.961 |
| Percentile coverage of P_med | 0.997–1.000 | 0.941–0.956 |
| Analytic SE / empirical SD | 0.65–0.83 | 0.98–1.02 |
| Bootstrap SE / empirical SD | 0.97–16.6 | 0.98–1.51 |

Manuscript consequence for Table 1: the "Wald under-covers at ~0.86–0.91 /
percentile ~1.00" reading came from the weight bug and should be replaced by the
post-fix numbers. The sobol near-null correction (the other half of this issue) is
already reflected in PR #37.

MCSE at nsim = 2000 is about 0.005 for coverage near 0.95.
