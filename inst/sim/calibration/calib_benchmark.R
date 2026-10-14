# Cost benchmark (task T15): 20 replications per cell at the planned controls,
# then projected CPU for the full matrices and the stability subset.
#
# Usage (from the package root):
#   Rscript inst/sim/calibration/calib_benchmark.R <outdir> [cores]
# Writes <outdir>/bench_task*.csv and inst/sim/results/calib_benchmark.csv.

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 2L)
source("inst/sim/calibration/calib_design.R")
cells <- calib_cells()

worker <- function(task) {
  system2("Rscript", c("-e", shQuote(sprintf(
    'assign("commandArgs", function(trailingOnly = TRUE) c("bench", "%d", "%s"), envir = globalenv()); source("inst/sim/calibration/calib_worker.R")',
    task, outdir))), stdout = TRUE, stderr = TRUE, env = "PROBMED_LOAD=source")
}
# Longest cells first so the pool stays busy.
order_tasks <- order(-(cells$family == "boot") * cells$N, -cells$N)
invisible(parallel::mclapply(order_tasks, worker, mc.cores = cores,
                             mc.preschedule = FALSE))

parts <- list.files(outdir, "^bench_task.*\\.csv$", full.names = TRUE)
d <- do.call(rbind, lapply(parts, utils::read.csv, stringsAsFactors = FALSE))
agg <- do.call(rbind, lapply(split(d, d$cell_id), function(z) data.frame(
  cell_id = z$cell_id[1], reps = nrow(z),
  secs_mean = mean(z$seconds), secs_max = max(z$seconds),
  original_ok = sum(z$original_fit_status == "ok"),
  original_separated = sum(z$original_fit_status == "separated"),
  original_failed = sum(z$original_fit_status == "failed"),
  interval_ok = sum(z$interval_status == "ok"),
  bootstrap_incomplete = sum(z$interval_status == "bootstrap_incomplete"),
  covered = sum(z$interval_status == "ok" & z$lower <= z$truth & z$upper >= z$truth),
  stringsAsFactors = FALSE)))
agg <- agg[match(cells$cell_id, agg$cell_id), ]
agg$projected_core_hours_main <- agg$secs_mean * CALIB_REPS / 3600
boot <- startsWith(agg$cell_id, "boot-")
# Stability subset: 100 reps per bootstrap cell at 4x n_boot and 4x n_sim.
agg$projected_core_hours_stability <- ifelse(boot, agg$secs_mean * 100 * 16 / 3600, 0)
utils::write.csv(agg, "inst/sim/results/calib_benchmark.csv", row.names = FALSE)
print(agg[, c("cell_id", "secs_mean", "original_separated", "bootstrap_incomplete",
              "projected_core_hours_main")], digits = 3, row.names = FALSE)
cat(sprintf("\nProjected core-hours: MBCO %.1f, bootstrap %.1f, stability %.1f, total %.1f\n",
            sum(agg$projected_core_hours_main[!boot]),
            sum(agg$projected_core_hours_main[boot]),
            sum(agg$projected_core_hours_stability),
            sum(agg$projected_core_hours_main) + sum(agg$projected_core_hours_stability)))
