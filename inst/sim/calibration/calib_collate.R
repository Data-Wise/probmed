# Collate calibration parts (task T17).
#
# Usage (from the package root):
#   Rscript inst/sim/calibration/calib_collate.R <partsdir> <outdir>
# Writes <outdir>/calibration.csv (the pmed summarizer's input schema, every
# attempted replication) and <outdir>/calibration_summary.csv.
#
# Decision rule (repair spec, frozen): per cell, an exact 95% binomial CI for
# coverage over ALL attempted replications, counting incomplete/invalid
# intervals as non-coverage. PASS if the CI lies inside [0.93, 0.97], FAIL if
# wholly outside, otherwise INCONCLUSIVE; direction under/over-coverage on
# FAIL. Conditional coverage (valid intervals only) is reported but never
# decides a cell.

args <- commandArgs(trailingOnly = TRUE)
partsdir <- args[1]
outdir <- args[2]
source("inst/sim/calibration/calib_design.R")
cells <- calib_cells()
# T22 reruns only the bootstrap cells: CALIB_FAMILY=boot restricts collation.
fam <- Sys.getenv("CALIB_FAMILY", "")
if (nzchar(fam)) cells <- cells[cells$family == fam, ]

parts <- list.files(partsdir, "^main_task.*\\.csv$", full.names = TRUE)
d <- do.call(rbind, lapply(parts, utils::read.csv, stringsAsFactors = FALSE))
d <- d[d$cell_id %in% cells$cell_id, ]
d <- d[order(match(d$cell_id, cells$cell_id), d$replication), ]

# Completeness: every cell, replications 1..2000 exactly once.
for (cid in cells$cell_id) {
  z <- d$replication[d$cell_id == cid]
  if (length(z) != CALIB_REPS || !setequal(z, seq_len(CALIB_REPS))) {
    stop("cell ", cid, " has ", length(z), " replications; expected 1..",
         CALIB_REPS, " exactly once")
  }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
schema <- c("cell_id", "replication", "truth", "lower", "upper", "converged",
            "original_fit_status", "interval_status", "lower_status", "upper_status")
utils::write.csv(d[, c(schema, setdiff(names(d), schema))],
                 file.path(outdir, "calibration.csv"), row.names = FALSE)

summ <- do.call(rbind, lapply(cells$cell_id, function(cid) {
  z <- d[d$cell_id == cid, ]
  ep <- if (startsWith(cid, "mbco-")) "root" else "percentile"
  valid <- z$converged & z$original_fit_status == "ok" & z$interval_status == "ok" &
    z$lower_status == ep & z$upper_status == ep & is.finite(z$lower) &
    is.finite(z$upper) & z$lower >= 0 & z$upper <= 1 & z$lower <= z$upper
  covered <- valid & z$lower <= z$truth & z$upper >= z$truth
  ci <- stats::binom.test(sum(covered), nrow(z))$conf.int
  decision <- if (ci[1] >= 0.93 && ci[2] <= 0.97) "PASS" else
    if (ci[2] < 0.93 || ci[1] > 0.97) "FAIL" else "INCONCLUSIVE"
  direction <- if (ci[2] < 0.93) "undercoverage" else
    if (ci[1] > 0.97) "overcoverage" else "none"
  est <- z$original_fit_status == "ok"
  data.frame(
    cell_id = cid, attempted = nrow(z), valid = sum(valid),
    original_separated = sum(z$original_fit_status == "separated"),
    original_failed = sum(z$original_fit_status == "failed"),
    bootstrap_incomplete = sum(z$interval_status == "bootstrap_incomplete"),
    disconnected = sum(z$interval_status == "disconnected"),
    truncated = sum(z$lower_status == "truncated" | z$upper_status == "truncated"),
    coverage = mean(covered), ci_lower = ci[1], ci_upper = ci[2],
    estimable_coverage = if (any(est)) mean(covered[est]) else NA_real_,
    conditional_coverage = if (any(valid)) mean(covered[valid]) else NA_real_,
    mean_width = if (any(valid)) mean(z$upper[valid] - z$lower[valid]) else NA_real_,
    decision = decision, direction = direction,
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(summ, file.path(outdir, "calibration_summary.csv"), row.names = FALSE)
print(summ[, c("cell_id", "valid", "coverage", "ci_lower", "ci_upper",
               "conditional_coverage", "decision", "direction")],
      digits = 3, row.names = FALSE)
cat("\nDecisions:", paste(names(table(summ$decision)), table(summ$decision),
                          collapse = ", "), "\n")
