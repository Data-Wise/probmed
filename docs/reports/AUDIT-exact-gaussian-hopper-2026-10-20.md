# Frozen Hopper source audit: exact-Gaussian pmed follow-up

**Decision:** The completed run's original six-cell gate is **PASS**, independently reproduced from its saved outer-replication rows. This audit verifies the *saved-data decision and provenance checks*, not the mathematical kernel or a new package method. The separate probmed I6 result remains **FAIL** and I7 remains held.

## Run identity and independent procedure

- Hopper job `4342657`, directory `/users/dtofighi/probmed_exact_20261019_full_e`.
- `code/manifest.csv`: 600 tasks × 200 outer replications = 120,000 rows over 20 cells; 2,000 nonparametric bootstrap draws per outer data set; seed base `20261021`.
- Run source files (`runner.R`, `kernel.R`, `config.R`, `analyze.R`, `submit_exact.sh`) have SHA-256 values matching every manifest row. The manifest SHA-256 is `2b392f36f88608b90260576a9565bfcc64feaf6a436abcb24ed11adbf4b5b872`; analyzer SHA-256 is `fba709c2e4f330f7b21a95645ddec4a64172e5328329b992e108f51089414c55`.
- `sha256sum --check --quiet part_hashes.sha256` returned exit 0 for the 600 saved parts. The independent script also checked all 600 per-part hash sidecars against that list, the complete manifest task windows, expected per-file schema, bootstrap counts, finite/bounded complete intervals, all `(scenario,N,rep)` keys, and the 20-cell summary.
- [Independent base-R audit script](audit_exact_gaussian_hopper.R) ran on Hopper under the run's `r/4.4.0-ytj2` module, reading the parts without sourcing `code/analyze.R` or modifying the remote run. It rebuilt every recomputable summary metric checked—including both 90% and 95% exact-binomial bounds—from 120,000 unique rows. All 20 cells had zero incomplete intervals and zero failed bootstrap fits. The saved `summary.csv` SHA-256 is `02f026751affb9321083e02c20ea7bde36e2dfd53ff6226489c886da6862b997`, the same as the earlier pmed report. `DECISION.txt` and `slurm-4342657.out` agree with the independent calculation.

The read-only verification commands were `sha256sum --check --quiet /users/dtofighi/probmed_exact_20261019_full_e/part_hashes.sha256` and, from this checkout, `ssh -o BatchMode=yes hopper 'source /etc/profile.d/modules.sh && module load r/4.4.0-ytj2 && Rscript --vanilla - /users/dtofighi/probmed_exact_20261019_full_e' < docs/reports/audit_exact_gaussian_hopper.R`. The auditor had two local setup corrections: Hopper's login shell needed the recorded R module loaded, and the saved hash-list uses one space between hash and path. Neither correction changed remote files or the gate.

## Six original acceptance cells

Each cell's closed-form null truth is 0.5. Coverage counts include all 6,000 attempts; an incomplete interval would be a miss. The frozen gate asks whether each **two-sided 90% exact-binomial interval for the coverage of a 95% percentile-bootstrap interval** is strictly inside 93%–97%.

| Scenario | N | Covered / attempted | Coverage | 90% exact bounds | Mean width | Mean point estimate | MCSE | Gate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| s6 | 250 | 5,719 / 6,000 | 95.317% | 94.843%–95.757% | 0.038005 | 0.499872 | 0.273% | PASS |
| s6 | 500 | 5,711 / 6,000 | 95.183% | 94.704%–95.630% | 0.025850 | 0.499947 | 0.276% | PASS |
| s6 | 1,000 | 5,663 / 6,000 | 94.383% | 93.870%–94.865% | 0.017866 | 0.499992 | 0.297% | PASS |
| s7 | 250 | 5,706 / 6,000 | 95.100% | 94.617%–95.550% | 0.053243 | 0.500321 | 0.279% | PASS |
| s7 | 500 | 5,674 / 6,000 | 94.567% | 94.061%–95.040% | 0.037105 | 0.500010 | 0.293% | PASS |
| s7 | 1,000 | 5,690 / 6,000 | 94.833% | 94.339%–95.295% | 0.026097 | 0.499930 | 0.286% | PASS |

The weakest lower limit is 93.870% (s6, N=1,000), still strictly above 93%. In the out-of-gate s6, N=50 cell, observed coverage was 98.133%; it remains a separate overcoverage disclosure. Width and mean point estimates were rebuilt and matched the saved summary, but neither is a retroactive gate.

## Limits and next work

This audit establishes that the current saved chunks, sidecars, manifest, summary, and decision are internally consistent. The current manifest/source hashes do not by themselves prove historical prespecification or guard against an earlier coordinated rewrite of all files; the prospective plan and job provenance are separate evidence. The audit did **not** independently derive/test the exact-Gaussian kernel against quadrature, replay original RNG streams, validate binary or other models, or test probmed/medsim package adapters. Those remain distinct tasks in the [follow-up plan](../PROPOSAL-exact-gaussian-followup-2026-10-19.md). Do not broaden the pmed manuscript's six-cell statement or lift its submission hold from this result.
