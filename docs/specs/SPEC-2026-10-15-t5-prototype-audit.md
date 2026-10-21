# SPEC — Reproducible audit of Gemini T5 composite-null prototypes

- **Date:** 2026-10-15
- **Tracker:** [probmed #50](https://github.com/Data-Wise/probmed/issues/50)
- **Capability map:** `reproducible-scoring` → `mc21-finite-sample` → `evidence-disposition`
- **Target:** `dev` via a feature branch and PR, following package branch policy
- **Current decision:** NO-GO; retain the chi-square(1) MBCO interval and `near_double_null` label
- **Scope:** make the proposed T5 scoring auditable, repair or exclude the finite-sample MC21 t cutoff mapping, and document evidence without reopening calibration or changing the package API

## Objective

Audit the 16 Gemini T5 prototype rules against the frozen lattice of 100 cells with 2,000 replications per cell (200,000 rows total). Reconstruct the finite-sample branch t statistics from the stored branch likelihood ratios, verify the MC21 branch/df mapping, and report results and theorem scope precisely enough to support the existing NO-GO decision.

Success means the scoring script reproduces all 400 baseline rows from the committed T4 score table, emits the expected 1,600 rule-by-cell rows with exact binomial intervals, and exposes whether any prototype passes the existing frozen decision rule. The prototypes remain research diagnostics; this work does not ship an adaptive interval or claim uniform validity.

## Assumptions

1. The committed replication data and T4 score table are the authoritative frozen inputs; no simulation rerun is required.
2. The stored branch likelihood ratios are Gaussian nested-model likelihood-ratio statistics, LR = N log(1 + t^2 / df). Reconstruct them as `t_a = sqrt(df_m * (exp(lr_a0 / N) - 1))` and `t_b = sqrt(df_y * (exp(lr_b0 / N) - 1))`, where `df_m = N - 2 - as.integer(cov)` and `df_y = N - 3 - as.integer(cov)`. Reproduce the 157/200,000 LR-versus-t branch-order disagreements from the committed data before applying the correction. If the method's derivation does not support this finite-sample mapping, exclude `MC21_t_all` and document why.
3. Preserve the T5 point thresholds and use exact 95% Clopper–Pearson intervals to classify Monte Carlo uncertainty under the explicit decision function below. Do not substitute `calib_collate.R`'s separate [0.93, 0.97] rule.
4. The external tracker for implementation tasks is issue #50. Existing `tasks/plan.md` and `tasks/todo.md` are active for other work and must not be overwritten.

## Commands

After implementation, the primary offline reproduction command will be:

```sh
Rscript inst/sim/calibration/calib_mbco_null_t5_score.R inst/sim/results/calibration/mbco_null/replications.csv.gz /tmp/t5_prototype_scores.csv
```

Focused verification:

```r
devtools::test(filter = "mbco-null-score")
```

Required package check for code changes:

```r
rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "error")
```

Do not run the 100-cell simulation harness or submit a full compute job for this specification.

## Project Structure

- `inst/sim/calibration/calib_mbco_null_t5_score.R` — T5 CLI; sources and reuses scoring helpers from `calib_mbco_null_score.R` without changing T4 `score_rules()` semantics
- `inst/sim/calibration/calib_mbco_null_score.R` — authoritative T4 baseline scoring helpers
- `inst/sim/results/calibration/mbco_null/` — frozen 2,000-replication-per-cell input, baseline T4 score table, and reproducible T5 score output
- `tests/testthat/test-mbco-null-score.R` — fixture and frozen-data regression coverage for baseline and T5 scorer
- `docs/reports/MBCO-NULL-T5-FEASIBILITY-2026-10-14.md` — evidence, corrected summaries, theorem-scope notes, and final disposition
- `docs/specs/` — this specification

Do not alter the raw replication data or baseline score table.

## Requirements and Acceptance Criteria

### Capability: reproducible-scoring

- The T5 CLI sources/reuses `calib_mbco_null_score.R` and accepts the frozen compressed replication CSV and an output path.
- It writes exactly 1,600 rows: 16 named rules × 100 lattice cells, with 2,000 replications per cell, accepted count, coverage, and exact 95% Clopper–Pearson intervals.
- For `chi2`, `switch_crit`, `switch_logN`, and `smooth`, accepted counts and valid replication counts match the 400 baseline rows exactly; coverage is derived from those counts, and interval columns match within absolute tolerance `1e-12`.
- The report identifies the input data, output table, command, rule definitions, and summary method.
- The script and output are present in version control; no dependency on the user's pasted attachment or temporary files remains.

### Capability: mc21-finite-sample

- Reconstruct finite-sample statistics from each branch LR using `t = sqrt(df * (exp(LR / N) - 1))`, with `df_m = N - 2 - as.integer(cov)` and `df_y = N - 3 - as.integer(cov)`; select the branch and df by `min(|t_a|, |t_b|)`, not LR ordering.
- Tests cover both df values, an LR/t ordering disagreement, and the t-cutoff-to-LR mapping for the selected branch.
- Reproduce the 157/200,000 LR-versus-t ordering disagreements from the committed data and verify that all are handled by t ordering. If the primary method does not justify the finite-sample mapping, explicitly exclude `MC21_t_all`.
- `MC21_t_05` remains separately scored from its direct rejection indicator; the `MC21_t_all` mapping must not alter it.

### Capability: evidence-disposition

- Reproduce all per-cell minima, origin summaries, null/alternative summaries, and uncertainty intervals directly from the saved score output.
- Distinguish a composite-null test at `p*=0.5` from a confidence interval requiring inversion over candidate values of `p*`.
- Label the HSX rules as analytical joint-significance analogues; do not call them the authors' adaptive bootstrap.
- Verify quotations, version, theorem, and assumptions against the cited primary sources. Do not infer full-null or confidence-set validity from a theorem about a narrower test class or a fixed target/local experiment.
- Preserve the frozen T5 gates: minimum coverage target 0.94 across the 100 cells, and origin coverage ceiling 0.97 across the four origin cells. For each cell, use its exact 95% Clopper–Pearson interval, without a familywise adjustment. For the minimum-coverage gate, classify a cell PASS if its lower bound is at least 0.94, FAIL if its upper bound is below 0.94, and INCONCLUSIVE otherwise; the gate is FAIL if any cell fails, PASS if all pass, and INCONCLUSIVE otherwise. For the origin ceiling, classify a cell PASS if its upper bound is at most 0.97, FAIL if its lower bound is above 0.97, and INCONCLUSIVE otherwise; aggregate the four origin cells by the same rule. Overall disposition is FAIL if either gate fails, PASS only if both pass, and INCONCLUSIVE otherwise. This is a pointwise Monte Carlo screening rule, not a familywise confidence guarantee. The [0.93, 0.97] interval rule in `calib_collate.R` is a separate procedure and does not apply here.
- If no prototype passes, retain the chi-square(1) interval and `near_double_null` label. Do not add a calibration API or revise manuscript claims.

## Code Style

Use base R for the standalone T5 CLI, explicit namespaces for non-base functions, deterministic ordering of rules and cells, and no new package dependencies. Source and reuse `calib_mbco_null_score.R` helpers rather than duplicate or alter T4 `score_rules()` semantics.

Example branch selection:

```r
df_m <- N - 2L - as.integer(cov)
df_y <- N - 3L - as.integer(cov)
t_a <- sqrt(df_m * (exp(lr_a0 / N) - 1))
t_b <- sqrt(df_y * (exp(lr_b0 / N) - 1))
branch_is_a0 <- t_a <= t_b
branch_df <- ifelse(branch_is_a0, df_m, df_y)
```

Tests must include a stored-data case where LR and reconstructed t order the branches differently.

## Testing Strategy

- Extend `tests/testthat/test-mbco-null-score.R` with synthetic cases that verify the reconstruction formula, branch selection, and cutoff mapping.
- Add a frozen-data regression check for the 1,600 output rows, 100 unique cells × 16 rules, exactly 2,000 replications per cell, baseline accepted-count/replication-count equality, interval tolerance `1e-12`, and the reproduced 157 ordering disagreements.
- Check deterministic output and no missing or duplicate rule-cell keys.
- Review the report against computed output and source theorem statements; a successful test run does not establish theorem applicability or full-interval validity.
- Run the focused test first, then the package check required by repository policy. No full simulation is part of this work.

## Boundaries

- **Always:** preserve the current NO-GO and package submission hold; keep raw simulation inputs unchanged; separate test validity from interval coverage; report Monte Carlo uncertainty.
- **Ask first:** launching a full simulation/Hopper job; changing acceptance thresholds; adding an adaptive interval API; changing manuscript language or submission status.
- **Never:** claim uniform composite-null or confidence-set validity without supporting proof; describe analytical HSX analogues as bootstrap implementations; modify unrelated active plans or task lists.

## Success Criteria

- A clean checkout can run the documented T5 command and regenerate the committed score table.
- All 400 T4 baseline accepted counts and replication counts reproduce exactly; confidence interval columns agree within `1e-12`.
- The 157 LR-versus-t ordering disagreements are reproduced and handled by the documented formula, or `MC21_t_all` is explicitly excluded.
- Per-cell and overall PASS/FAIL/INCONCLUSIVE labels follow the specified gate logic.
- The feasibility report still records NO-GO unless a prototype passes both frozen gates.
- Focused tests and required package checks pass; no full simulation, API change, or submission action occurs.

## Open Questions

- Does the primary MC21 construction justify applying this finite-sample t-scale mapping to the interval inversion? Resolve this before retaining `MC21_t_all`; if not established, exclude it and report its scores only as an unvalidated diagnostic.
