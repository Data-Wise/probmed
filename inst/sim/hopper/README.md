# Hopper (HPC) gauge grids

This directory holds two independent SLURM studies:

| Grid | Scripts | Question |
|---|---|---|
| **Coverage** (done 2026-06-23) | `run_gauge_boot_grid.R`, `collate_gauge_boot.R`, `submit_gauge_boot.sh`, `submit_collate.sh` | What is the CI coverage of `W` / `P_med`? |
| **#11 threshold validation** (not yet run) | `run_weakid_validation.R`, `collate_weakid_validation.R`, `submit_weakid_validation.sh`, `submit_weakid_collate.sh` | Are `weak_id_ratio_threshold = 3` / `oe_snr_threshold = 2` the right defaults? |

---

## #11 threshold-validation grid (`run_weakid_validation.R`)

**Status: scripts ready, NOT submitted.** 192 tasks consume a real share of the
CPU allocation and queue behind other work -- submit deliberately.

**Why a second grid.** The coverage grid varies `n x tint x binY`, which lands
`oe_snr` at the *extremes* (~1 or ~12). A threshold arbitrates in the *middle*
(`oe_snr` ~1.5-5), which no coverage cell samples. This grid sweeps the A-effect
**scale `s`**, the axis that moves `oe_snr` monotonically through the crossover.

**Design:** `s` {0.05,0.10,0.15,0.20,0.30,0.50} x `n` {800,3000} x `binY` {F,T}
= **24 cells** x 8 chunks x 250 reps = **nsim 2000/cell** (48,000 reps).
`B = 200` -- deliberately `ward_residual()`'s **default**, since the threshold
governs a user-facing flag and must be calibrated under the conditions users
actually run (a probe showed `B=999` raises the ratio only ~6-8%).

**Records the shipped fields verbatim** (`weak_id`, `weak_id_ratio`, `oe_snr`,
`oe_regular`) rather than a re-derived proxy, so it validates exactly what users
get. It also stores `wid_wald`/`wid_pct`, so the threshold sweep can be redone
without re-running the grid.

**The analysis (`collate_weakid_validation.R`).** Percentile coverage is ~1.00
everywhere, so the flag *cannot* be validated against percentile under-coverage.
Its real claim is "W's CI is least trustworthy here", and the CI users get **by
default is the analytic Wald** one (~0.86-0.91). So the test is an operating
characteristic: **is Wald coverage materially worse when the flag fires?** The
collator sweeps candidate thresholds and reports `covW | flagged` vs
`covW | unflagged`; a validated threshold leaves the unflagged group ~nominal and
the flagged group clearly degraded. It also reports sensitivity/FPR at the
shipped defaults.

**Indications going in** (`../pilot/weak_id_pilot.R`, 120 draws, n=1500,
continuous Y) -- these are *why the grid is needed*, not results it may assume:

| per-draw | pilot | exact 95% CI |
|---|---|---|
| false alarms (`flag \| oe_snr >= 3`) | 0/58 | `[0.00, 0.06]` |
| detection (`flag \| oe_snr <= 1.2`) | 12/28 = 0.43 | `[0.24, 0.63]` |

So the flag looks specific but insensitive -- but the detection interval spans
0.24 to 0.63, which is not a measurement of anything, and 34 further draws sit
between the regimes and are summarised by neither row. The firmer part is the
mechanism: `weak_id_ratio`'s draw-to-draw SD approaches its mean in the weak
regime, which necessarily costs detection. A1 (`oe_regular`) and A2 (`weak_id`)
are *observed* to fire on different draws, suggesting they are complementary
rather than redundant. The grid's job is to turn every one of these hedges into
a number -- at nsim=2000/cell the same rates carry a CI of roughly +/-0.02.

### Why the existing coverage grid cannot substitute for this one

Tempting shortcut, ruled out empirically. The coverage grid's raw data
(`~/gauge_boot/gauge_boot_raw.rds`, 16,000 reps) stores **no CI widths** -- only
`seW_an`, `seW_bt`, and `divW_WaldVsPct` (a max endpoint *displacement*, not a
ratio). The Wald width is recoverable (`2*z*seW_an`, symmetric by construction)
but the **percentile width is not**: it is quantile-based, and `seW_bt` is a
replicate SD, not the quantile spread -- and the skew that separates them is
precisely the phenomenon under study.

You can still run the threshold sweep on that data using
`se_ratio = seW_bt/seW_an` as a stand-in, and it looks encouraging (separation
peaks near 3). **Do not trust that transfer.** A direct check (64 reps spanning
the identification range) found `se_ratio` correlates with the shipped
`weak_id_ratio` at only **Pearson r ~ 0.41 / Spearman ~ 0.53**, on a different
scale (medians 2.39 vs 1.66), agreeing on the ">= 3" call just ~72% of the time.
A threshold calibrated on `se_ratio` therefore says little about the same
threshold on the shipped width ratio. Hence this grid, which records the shipped
fields verbatim.

### Running it

The sbatch files invoke `$HOME/weakid_val/<script>.R`, so the scripts must be
**deployed to the cluster first** — `~/weakid_val/` does not exist until you make
it. Skipping this queues 192 tasks that each fail instantly with "file not found".

```bash
# 1. deploy (from a checkout of this package, on your workstation)
ssh hopper 'mkdir -p ~/weakid_val/{parts,logs}'
scp inst/sim/hopper/run_weakid_validation.R \
    inst/sim/hopper/collate_weakid_validation.R \
    inst/sim/hopper/submit_weakid_validation.sh \
    inst/sim/hopper/submit_weakid_collate.sh   hopper:~/weakid_val/

# 2. install probmed >= PR #23 onto R_LIBS (the weak-ID fields must exist);
#    both the runner and the collator abort loudly if they do not.

# 3. submit
ssh hopper 'cd ~/weakid_val && sbatch submit_weakid_validation.sh'   # 192-task array
ssh hopper 'cd ~/weakid_val && sbatch submit_weakid_collate.sh'      # AFTER it completes
```

Sanity-check one task before committing the whole array:

```bash
ssh hopper 'cd ~/weakid_val && SLURM_ARRAY_TASK_ID=1 Rscript run_weakid_validation.R'
```

Outputs land in `~/weakid_val/parts/`; the collator writes
`weakid_validation_cells.csv`, `weakid_threshold_sweep.csv`,
`oe_snr_threshold_sweep.csv`, and `weakid_validation_raw.rds` to `~/weakid_val/`.
Copy the CSVs back into `../results/` and commit them **together with any script
change**, or the numbers become unreproducible (the mistake this directory exists
to prevent — see the coverage grid's history below).

---

## Gauge coverage grid

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
