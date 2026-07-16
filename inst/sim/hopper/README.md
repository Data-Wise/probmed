# Hopper (HPC) gauge coverage grid

The bootstrap arm of the gauge coverage study is CPU-expensive (each rep costs
`B` bootstrap refits x cross-fitting), so it runs as a SLURM array on **hopper**
(UNM CARC) rather than locally. These are the scripts that produced
`../results/gauge_boot_coverage_nsim2000.csv`; they previously lived only in
`~/gauge_boot/` on the cluster, making the committed results unreproducible.

## Files

| File | Role |
|---|---|
| `run_gauge_boot_grid.R` | The array worker. One `$SLURM_ARRAY_TASK_ID` (1..64) -> one (cell, chunk) pair; writes `parts/part_cell<ci>_chunk<ch>.rds`. |
| `submit_gauge_boot.sh` | `sbatch` for the 64-task array (`--array=1-64`, 12h, 4G). |
| `collate_gauge_boot.R` | Aggregates all `parts/*.rds` -> per-cell coverage + MCSE; writes the CSV and `gauge_boot_raw.rds`. |
| `submit_collate.sh` | `sbatch` for the collate step. |

## Design

- **Grid:** 8 cells = `n` {800, 3000} x `tint` {0, 0.8} x `binY` {FALSE, TRUE}.
- **Reps:** 8 chunks x 250 = **nsim = 2000** per cell (16,000 reps total),
  `B = 999` bootstrap resamples.
- **Truth:** Monte-Carlo at `N = 2e6` (not closed-form -- the binary-outcome
  cells have no closed form, unlike the local analytic grid).
- **Array split:** cell = `((aid-1) %/% 8) + 1`, chunk = `((aid-1) %% 8) + 1`.

## Relationship to the local grid

This is **not** the same design as the in-repo local study
(`../gauge_coverage.R`, 24 cells: n {500,1000,2000,4000} x tau {0,0.2,0.8} x
se_method, nsim=1000, closed-form truth, continuous Y only). The two are
complementary:

| | Local (`analytic_coverage_nsim1000.csv`) | Hopper (`gauge_boot_coverage_nsim2000.csv`) |
|---|---|---|
| Cells | 24 (4 n x 3 tau x 2 se_method) | 8 (2 n x 2 tint x 2 binY) |
| nsim | 1000 | 2000 |
| Truth | closed-form | Monte-Carlo (N=2e6) |
| Binary Y | no | **yes** |
| Arms reported | analytic only (committed half) | analytic **and** percentile, per row |

## Two guards worth knowing about

1. **Stale-package guard** (`run_gauge_boot_grid.R`, top): aborts the job unless
   `se_method="bootstrap"` genuinely differs from `"analytic"`. Without it, a
   wrong `R_LIBS` ordering silently makes the bootstrap arm a no-op and the
   "percentile" coverage would really be the analytic CI.
2. **Collate-time guard** (`collate_gauge_boot.R`): re-checks the same invariant
   across all reps and writes `STALE_WARNING.txt` if bootstrap ~= analytic.
   The committed run passed: `bootstrap==analytic in 0.0% of reps; median
   seW_bt/seW_an = 1.82`.

## Reproducing

```bash
# on hopper, with the gauge feature-branch probmed installed into $HOME/Rlib/4.4-gauge
sbatch submit_gauge_boot.sh     # 64-task array
sbatch submit_collate.sh        # after the array completes
```

Note both `.sh` files hardcode `/users/dtofighi/...` paths and
`module load r/4.4.0-ytj2` + `R_LIBS=$HOME/Rlib/4.4-gauge:$HOME/Rlib/4.4-a15`;
adjust for another account/cluster.
