# Hopper (HPC) gauge grids

This directory holds two independent SLURM studies:

| Grid | Scripts | Question |
|---|---|---|
| **Coverage** (done 2026-06-23) | `run_gauge_boot_grid.R`, `collate_gauge_boot.R`, `submit_gauge_boot.sh`, `submit_collate.sh` | What is the CI coverage of `W` / `P_med`? |
| **#11 threshold validation** (submitted 2026-07-16, job `4277259`; not yet collated) | `run_weakid_validation.R`, `collate_weakid_validation.R`, `submit_weakid_validation.sh`, `submit_weakid_collate.sh` | Are `weak_id_ratio_threshold = 3` / `oe_snr_threshold = 2` the right defaults? |

---

## #11 threshold-validation grid (`run_weakid_validation.R`)

**Status: COMPLETE 2026-07-16** (job `4277259`, 192/192 tasks, 0 failures; collate
`4277455`). CSVs in `../results/`: `weakid_validation_cells.csv`,
`weakid_threshold_sweep.csv`, `oe_snr_threshold_sweep.csv`; raw per-rep data
(`weakid_validation_raw.rds`, 2 MB) stays on hopper in `~/weakid_val/`.
**Findings below** ("What the grid found") — read them before citing either
threshold; the results validate *less* than this section's original test
anticipated, and one gate's semantics inverted.

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
`covW | unflagged`. *(Post-run caveats: this criterion turned out to be the wrong
test for A1 — see the findings — and the collator's pooled conditioning
understates the A2 contrast ~3x relative to within-cell; fixing it to report
pooled + within-cell + per-cell is an open item.)*

### What the grid found (2026-07-16, adversarially re-verified 2026-07-19)

All numbers recomputed independently from the raw 48,000 reps by a second
analysis attacking the first; where they disagreed, the numbers here are the
corrected ones.

**A2 (`weak_id`, ratio >= 3): direction real, but ~97% width-tautology.**
Flagged draws have worse Wald coverage in **24/24 cells** (each individually
> 2xSE; sign-test p = 6e-8; pooled contrast -0.016, within-cell -0.030..-0.073
depending on cell weighting). But `weak_id_ratio = wid_pct/wid_wald` shares
`wid_wald` with the outcome: in a regression of `covW_an` on flag + cell fixed
effects, the flag coefficient (-0.051) collapses to **-0.002 (SE 0.003)** under a
flexible `log(wid_wald)` control — while controlling `log(wid_pct)` instead makes
it *stronger*. The bootstrap numerator contributes nothing detectable; the flag
operates purely through "the analytic interval is narrow". The honest claim is:
narrow Wald intervals under-cover in this DGP, and the ratio is a usable
**self-contained proxy** for narrowness (a user has no cell reference to judge
"narrow" against). It is NOT evidence that the bootstrap-vs-Wald comparison
detects a pathology beyond width. Thresholds **1.75-4 are not distinguished**;
apparent flatness of the separation was an artifact of min-n weighting (under
equal cell weights separation grows monotonically, -0.056 -> -0.079 over
t = 1.75 -> 4). **t = 3 is a convention**, kept as one.

**A1 (`oe_regular`, `oe_snr` < 2): the stronger diagnostic — with semantics
opposite to its documentation.** Flagged draws **over**-cover: +0.10 (22/22
qualifying cells; stable +0.099..+0.104 across binY, s-range, and n slices).
Mechanism verified: A1 fires as `OE -> 0`, where flagged Wald intervals are
2.0x-5.2x wider (per-cell median) and median **21x wider than |truth|** vs 3x
unflagged. Unlike A2, A1 **survives the width control** (+0.036, SE 0.003
residual effect), i.e. it carries information beyond width alone. So A1 flags
intervals untrustworthy in the *width/informativeness* sense, not the coverage
sense — "flagged group degraded in coverage" was the wrong criterion for it.

**Scope:** all 24 cells share one DGP family (the same one as the coverage
grids). Nothing here supports transferring either numeric threshold to other
DGPs.

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

Two more failure modes, both found the hard way on 2026-07-16 (job `4277033`:
192/192 FAILED, exit 127, 0s elapsed — see step 2 and the note below step 3):

- **The installed package predates the weak-ID fields.** `~/Rlib/4.4-gauge` holds
  probmed 0.2.0.9000, which has none of them; the runner's guard aborts every task.
- **`module` is undefined in SLURM's batch shell.** Fixed in the scripts, but the
  same trap applies to any new sbatch file here — see the note below.

```bash
# 1. deploy (from a checkout of this package, on your workstation)
ssh hopper 'mkdir -p ~/weakid_val/{parts,logs}'
scp inst/sim/hopper/run_weakid_validation.R \
    inst/sim/hopper/collate_weakid_validation.R \
    inst/sim/hopper/submit_weakid_validation.sh \
    inst/sim/hopper/submit_weakid_collate.sh   hopper:~/weakid_val/

# 2. install probmed >= 0.3.0 into ~/Rlib/4.4-weakid (the weak-ID fields must
#    exist; both the runner and the collator abort loudly if they do not).
#    Install into 4.4-weakid, NOT 4.4-gauge -- the latter still holds the
#    0.2.0.9000 that produced the coverage grid and is kept intact for its
#    reproducibility. Both submit scripts already put 4.4-weakid first.
R CMD build --no-build-vignettes --no-manual .        # do NOT document() first:
                                                      # roxygen 8.0.0 corrupts the
                                                      # S7 docs (DESCRIPTION pins 7.3.3)
scp probmed_*.tar.gz hopper:/tmp/
ssh hopper 'bash -lc "mkdir -p \$HOME/Rlib/4.4-weakid && module load r/4.4.0-ytj2 \
  && R CMD INSTALL --library=\$HOME/Rlib/4.4-weakid /tmp/probmed_*.tar.gz"'
# verify the fields actually landed. The version string alone does NOT prove it --
# check the properties themselves. Write the probe to a file rather than fighting
# nested ssh/R quoting:
cat > /tmp/probe.R <<'EOF'
need <- c("weak_id", "weak_id_ratio", "oe_snr", "oe_regular",
          "W_ci_wald", "weak_id_ratio_threshold", "oe_snr_threshold")
cat("MISSING:", paste(setdiff(need, names(probmed::GaugePmedResult@properties)),
                      collapse = ", "), "\n")
EOF
scp /tmp/probe.R hopper:/tmp/probe.R
ssh hopper 'bash -lc "module load r/4.4.0-ytj2 \
  && export R_LIBS=\$HOME/Rlib/4.4-weakid:\$HOME/Rlib/4.4-gauge:\$HOME/Rlib/4.4-a15 \
  && Rscript /tmp/probe.R"'        # must print an empty MISSING: list

# 3. submit
ssh hopper 'cd ~/weakid_val && sbatch submit_weakid_validation.sh'   # 192-task array
ssh hopper 'cd ~/weakid_val && sbatch submit_weakid_collate.sh'      # AFTER it completes
```

**Note on `module`.** SLURM runs the batch script in a *non-login* shell, where
`module` does not exist. The scripts here therefore `source /etc/profile.d/modules.sh`
before `module load`, and hard-fail if `Rscript` is still off `PATH`. The earlier
form (`module load ... 2>/dev/null`) worked only because `sbatch` exports the
submitter's environment, so submitting from a login shell that had *already*
module-loaded R masked the bug; submitting the same script over plain `ssh` killed
all 192 tasks in 0s. Keep the `source` line in any new sbatch file here.

**Always pilot before submitting the full array — in BOTH tiers.** Standing rule.
The two tiers catch disjoint failure classes, and tier 1 alone is what let job
4277033 die 192/192: it passed cleanly on the login node while the real job could
not even start.

*Tier 1 — login-node probe (logic).* Catches estimand/output/arithmetic errors.
`REPS_PER` is 250, so running a task verbatim is a *full chunk* (hours), not a
smoke test — use a reduced-rep copy writing to a throwaway dir. Check **both**
branches: task 1 (`binY=FALSE`) and task 97 (`binY=TRUE`, first binary cell).
Task 1 never touches the binary quadrature, so it alone proves nothing there.

```bash
ssh hopper 'cd ~/weakid_val && mkdir -p /tmp/probe_parts \
  && sed -e "s/REPS_PER <- 250L/REPS_PER <- 2L/" \
         -e "s|~/weakid_val/parts|/tmp/probe_parts|g" \
         run_weakid_validation.R > /tmp/probe_run.R'
# then, with R_LIBS set as in step 2, for TASK in 1 and 97:
#   SLURM_ARRAY_TASK_ID=$TASK Rscript /tmp/probe_run.R
```

*Tier 2 — `sbatch` pilot (environment).* **Does not substitute for tier 1, and is
not substituted by it.** A login shell has `module` defined and R on `PATH`; a
SLURM batch shell has neither. Only a real submission tests module init, `PATH`,
`R_LIBS` resolution, node filesystem, and memory. Submit the actual script with a
one-task array and confirm it reaches `RUNNING` with **non-zero elapsed** — an
instant `FAILED` at `00:00:00` is the signature of an environment bug, not a code
bug (exit 127 = command not found).

```bash
ssh hopper 'cd ~/weakid_val \
  && sed "s/--array=1-192/--array=1-1/" submit_weakid_validation.sh > probe_submit.sh \
  && sbatch probe_submit.sh'
ssh hopper 'sacct -j <id> --format=JobID%14,State%12,ExitCode,Elapsed -X'
# RUNNING at 00:00:10 => environment OK. Then scancel and submit the full array.
```

When reading logs for failures, do not grep only for `error|abort|cannot` —
`Rscript: command not found` matches none of those and an all-failed run looks
clean.

Known-answer check on the output: at `s = 0.05` the exact truth is
`trW = 0.01287554` (continuous) and `0.0100586` (binary), and `trW` must be a
**single** value per cell — chunk-to-chunk variation would mean the MC truth is
back (see the git history of `run_weakid_validation.R`).

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

**The two share a data-generating process** — an earlier version of this section
claimed they were "not the same design ... complementary", which is wrong and
undersells them. `../gauge_coverage.R`'s `.pars` is
`g0=-0.2, gC=0.8, aM=0.6, aC=0.4, bA=0.5, bM=0.7, bC=0.3`, and
`run_gauge_boot_grid.R` hardcodes exactly those (`expit(-0.2+0.8*C)`,
`M <- 0.6*A + 0.4*C`, `lin <- 0.5*A + 0.7*M + tint*A*M + 0.3*C`). **`tau` and
`tint` are the same parameter** (the A x M coefficient); both grids reduce to
`OE = 0.92 + 0.6*tau`, so `W = 0.6*tau / (0.92 + 0.6*tau)` in both. They differ
only in `n`, `nsim`, truth method, and which arms were kept.

| | Local (`analytic_coverage_nsim1000.csv`) | Hopper (`gauge_boot_coverage_nsim2000.csv`) |
|---|---|---|
| Cells | 24 (4 n x 3 tau x 2 se_method) | 8 (2 n x 2 tint x 2 binY) |
| `n` | 500, 1000, 2000, 4000 | 800, 3000 (**no overlap with local**) |
| nsim | 1000 (MCSE ~0.010) | 2000 (MCSE ~0.007) |
| Truth | closed-form | Monte-Carlo (N=2e6) |
| Binary Y | no | **yes** |
| Arms reported | analytic only (committed half) | analytic **and** percentile, per row |

**They cross-validate.** Because the DGP is shared, the continuous-Y analytic arm
is the *same quantity* measured twice, by independent truth methods, at
interleaved `n`. At `tau = 0` the local grid gives 0.887 / 0.866 / 0.855 / 0.857
(n = 500 / 1000 / 2000 / 4000) and hopper gives 0.867 (n=800) and 0.859 (n=3000) —
on one curve. Closed-form and MC truth agreeing is real evidence, and it is the
strongest check either file contains. (Consistent with hopper's stored `trW` being
off by only ~0.002: exact `W` at `tint=0.8` is `0.48/1.4 = 0.342857`, stored as
0.341.)

**Building the manuscript table: never aggregate.** Rows are `tau x n`; each cell
stands alone. Pooling these 20 continuous cells yields ~0.88 and **hides that
coverage degrades as n grows** (0.887 -> 0.857 at `tau = 0`) — the marginal
reverses the conditional, the same Simpson's-paradox trap that makes the `weak_id`
flag look inverted when read across cells instead of within them. Keep binary-Y in
a **separate panel**: `trW` there is -0.035 at `tint=0` (not 0), so those are a
different estimand value, not comparable cells. And local's percentile column is
*absent*, not zero — only its analytic half was committed.

**The finding both grids carry:** analytic/Wald coverage for `W` is sub-nominal in
**every** continuous cell (0.855-0.923, n=20 cells), never reaching 0.95, and it
gets *worse* with n. Percentile over-covers (~1.000). This is metric-independent
and does not depend on any threshold.

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
`module load r/4.4.0-ytj2`; the working directory and library path default to
`GAUGE_BOOT_DIR=$HOME/gauge_boot` and `GAUGE_R_LIBS=$HOME/Rlib/4.4-gauge:$HOME/Rlib/4.4-a15`
(both overridable from the submitting shell — `sbatch` exports its environment);
adjust for another account/cluster.

## Post-fix rerun (2026-08-22, issue #38)

Every number above was produced with the corner-EIF weight bug (PR #34: one
propensity per fold, not per row). `~/gauge_boot/` and `~/Rlib/4.4-gauge` /
`4.4-weakid` stay frozen so the committed CSV remains reproducible; the rerun
uses the same 64-task design from a separate directory and library:

```bash
# dev tarball (R CMD build --no-build-vignettes --no-manual) installed with
#   R_LIBS=$HOME/Rlib/4.4-postfix:$HOME/Rlib/4.4-a15 R CMD INSTALL --library=$HOME/Rlib/4.4-postfix probmed_*.tar.gz
# verified by a known-answer probe: 4.4-postfix reproduces the local dev numbers
# to 1e-10 (W = 0.0219254952 at the probe seed); 4.4-gauge / 4.4-weakid give
# W = 0.2096424849 (the pre-fix weights) -- positive control.
D=$HOME/gauge_boot_postfix; mkdir -p $D/logs $D/parts   # scripts copied into $D
GAUGE_BOOT_DIR=$D GAUGE_R_LIBS=$HOME/Rlib/4.4-postfix:$HOME/Rlib/4.4-a15 \
  sbatch --account=2016507 --output=$D/logs/gb_%A_%a.out --array=1-64 $D/submit_gauge_boot.sh
GAUGE_BOOT_DIR=$D GAUGE_R_LIBS=$HOME/Rlib/4.4-postfix:$HOME/Rlib/4.4-a15 \
  sbatch --account=2016507 --output=$D/logs/collate_%j.out $D/submit_collate.sh   # after the array
```

Pilot first (`--array=1,33` covers a continuous and a binary cell), per the
hopper memory: an instant `FAILED` at `00:00:00` is an environment bug, not a
code bug. Expected from the package-side measurement (n = 800, 250 datasets):
Wald ~0.94-0.97 in regular cells, percentile ~0.90-0.93, bias unchanged.
Collated output lands in `$D/gauge_boot_grid.csv`; commit it beside
`../results/gauge_boot_coverage_nsim2000.csv` as `..._postfix.csv`, do not
overwrite.
