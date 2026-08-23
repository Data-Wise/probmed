# Handoff: gauge manuscript claims to correct after the corner-EIF weight fix

**For:** `~/projects/research/pmed-modern/01-gauge-pmed/manuscript/gauge-pmed.qmd`
(cross-repo; this file only lists what to change and what is already measured —
the manuscript edit and the grid rerun are the author's).
**Source of the numbers:** probmed PR #34 (`a73efc8`, the fix), PR #35
(`8fbfa43`, bootstrap arm + gates re-measured), `docs/specs/FINDINGS-2026-08-22-postfix-bootstrap-gates.md`,
`inst/sim/results/se_shipped_postfix.csv`, `inst/sim/results/boot_gates_postfix.csv`.

## What was wrong

Every number in Table 1 and every coverage sentence in the abstract, §Results,
§Discussion and the Supplement was produced by an estimator whose inverse-
probability weights were one constant per cross-fit fold: `ifelse(z == 1, p1,
1 - p1)` with scalar `z` returned the first fold row's propensity (and
mediator-density proxy) for every row, from the first gauge commit until
2026-08-22. The point estimates were not visibly biased (triply-robust EIF,
correctly specified outcome / projection models), but the influence function was
not the efficient one and its variance depended on which row landed first in each
fold. Table 1 was produced by the hopper boot grid (`inst/sim/hopper/run_gauge_boot_grid.R`,
job 4239663) running the *installed* probmed in `~/Rlib/4.4-gauge`, which carries the
bug (so does `4.4-weakid`). The fixed dev package is installed as `~/Rlib/4.4-postfix`
(known-answer probe matches local dev to 1e-10; the old libs give the pre-fix value)
and the same 64-task grid was relaunched 2026-08-22 21:53 MDT (jobs 4311422 + 4311424,
`~/gauge_boot_postfix/`; README "Post-fix rerun"). The local standalone
(`inst/sim/gauge_coverage_standalone.R`, continuous Y only, a different n grid) had the
same bug and is fixed, verified identical to the package (`inst/sim/check_standalone_equiv.R`).

## Claims to change, with what replaces them

| where (qmd line, 2026-08-22 copy) | claim | post-fix measurement (n = 800, 250 datasets, exact truth) | action |
|---|---|---|---|
| abstract 24-26; l. 83-84; l. 217-225; Table 1 l. 250-261; Fig. (b) l. 263 | analytic coverage of W 0.86-0.91, anti-conservative; percentile restores "nominal-or-above" (1.00 everywhere) | Wald **0.936** (cont, tau = 0) / **0.972** (bin, tau = 0); se/empSD 0.89 / 0.99 with dispersion matching the oracle's; percentile **0.905 / 0.930** at the same width as Wald, bootstrap SD unstable for binary Y (CV 2.2) | rerun Table 1 with the fixed standalone (all 8 cells, 2000 reps); rewrite the sentences: Wald calibrated, percentile no better |
| l. 222 "estimated SE under-estimating the empirical SD throughout" | 30% se shortfall | se/empSD median 0.89-0.99 | delete; was the bug |
| l. 232-235, 422-424 | ratios ⇒ report the tail-aware percentile interval | no longer supported | make the analytic Wald the reported interval; bootstrap optional |
| l. 242, 458 | `reps > 1` "removes the fold-split component of the variance" | the component was the bug; post-fix empSD 0.063 (reps 1) vs 0.061 (reps 10), oracle 0.056 | keep `reps` as a minor efficiency option or drop the paragraph |
| l. 449, 472 (ii) | weak-identification flag = Wald-vs-percentile width ratio ≥ 3 | flagged draws cover like unflagged at every threshold 1.2-3; tracks bootstrap instability, which the OE-regularity gate already detects; 2/98 beyond it | demote to a one-line bootstrap-instability note or remove; keep (i) the OE-regularity gate |
| §Application l. 280 (near-null, Fieller) | Wald interval uninformative near `OE = 0`; Fieller set unbounded | stands: near-null cell Wald 0.945 overall, `oe_regular = FALSE` draws 0.99 at 8.7x \|W\|, oracle se explodes too | keep |
| §Simulation design l. 460-468 (MCSE) | n_sim ≥ 1000 to separate conservative from nominal | unaffected | keep |
| Supplement: bootstrap validity via Lin & Han (2026) | already corrected 2026-08-22 (refit scheme not covered) | unchanged | keep the correction |

## What still needs the manuscript's own grid

The package measurements are 250 datasets per cell at n = 800 (and 200 x B =
200 for the bootstrap arm). Table 1 needs all 8 manuscript cells at 2000 reps,
n = 800 and 3000, both Y types, both arms: rerun `inst/sim/hopper/` (SLURM array;
pilot with `--array=1-1` first, per the hopper memory), then
`inst/sim/collect_gauge.R`. Expected: Wald ~0.94-0.97 everywhere regular;
percentile ~0.90-0.93; bias unchanged (it was never the weights' problem).

## Probmed-side references for the text

- `?ward_residual` `se_method` and `reps`; `?GaugePmedResult` `weak_id`,
  `oe_regular` (post-fix wording, PR #35).
- NEWS 0.3.0.9000: "Bug fix: corner-EIF inverse-probability weights were one
  constant per fold" and "`weak_id` no longer warns".
- Issues #11 (closed: A1 stands, A2 uninformative) and #32 (closed: BCa not
  planned) carry the public record of the corrections.
