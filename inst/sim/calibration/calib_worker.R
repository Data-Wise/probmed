# Calibration worker (tasks T16-T17). One call runs one chunk of one cell.
#
# Usage (from the package root, with probmed installed or loadable):
#   Rscript inst/sim/calibration/calib_worker.R <mode> <task_id> <outdir> [chunks_per_cell]
#     mode = main       36 cells x chunks_per_cell tasks (default 20 -> 100 reps each)
#     mode = stability  15 bootstrap cells x 1 task, reps 1-100, n_boot 4000, n_sim 40000
#     mode = bench      36 tasks, reps 1-20 per cell at the planned controls
#
# Output: <outdir>/<mode>_task<id>.csv with the columns the pmed summarizer
# requires (cell_id, replication, truth, lower, upper, converged,
# original_fit_status, interval_status, lower_status, upper_status) plus
# estimate, timing and failure diagnostics.
#
# Set PROBMED_LOAD=source to use pkgload::load_all(".") instead of
# library(probmed).

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
task <- as.integer(args[2])
outdir <- args[3]
chunks <- if (length(args) >= 4) as.integer(args[4]) else 20L

if (identical(Sys.getenv("PROBMED_LOAD"), "source")) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressMessages(library(probmed))
}
source("inst/sim/calibration/calib_design.R")
fixtures <- load_calib_fixtures("inst/sim/calibration/calib_fixtures.dcf")
cells <- calib_cells()

plan <- switch(mode,
  main = {
    ci <- (task - 1L) %/% chunks + 1L
    per <- CALIB_REPS %/% chunks
    ch <- (task - 1L) %% chunks
    list(ci = ci, reps = ch * per + seq_len(per), n_boot = BOOT_N_BOOT, n_sim = BOOT_N_SIM)
  },
  stability = {
    boot_ci <- which(cells$family == "boot")
    list(ci = boot_ci[task], reps = 1:100, n_boot = 4L * BOOT_N_BOOT, n_sim = 4L * BOOT_N_SIM)
  },
  bench = list(ci = task, reps = 1:20, n_boot = BOOT_N_BOOT, n_sim = BOOT_N_SIM),
  stop("unknown mode: ", mode)
)
cell <- cells[plan$ci, ]
g <- (plan$ci - 1L) * CALIB_REPS + plan$reps
streams <- calib_streams(min(g), max(g))

map_mbco <- function(r) {
  dg <- r@diagnostics$pmed
  st <- dg$status
  ep <- if (!is.null(dg$endpoint_status)) unname(dg$endpoint_status) else c("failed", "failed")
  interval <- switch(st, ok = "ok", disconnected = "disconnected", "failed")
  list(interval = interval, lower_status = ep[1], upper_status = ep[2])
}

one_rep <- function(rep, stream) {
  assign(".Random.seed", stream, envir = globalenv())
  t0 <- proc.time()[["elapsed"]]
  out <- list(cell_id = cell$cell_id, replication = rep, truth = NA_real_,
              estimate = NA_real_, lower = NA_real_, upper = NA_real_,
              converged = FALSE, original_fit_status = "failed",
              interval_status = "failed", lower_status = "failed",
              upper_status = "failed", n_valid = NA_integer_,
              fail_separated = NA_integer_, fail_other = NA_integer_,
              error = "")
  if (cell$family == "mbco") {
    d <- gen_mbco(cell$k, cell$N)
    out$truth <- truth_mbco(cell$k)
    r <- tryCatch(pmed(Y ~ X + M, formula_m = M ~ X, data = d, treatment = "X",
                       mediator = "M", method = "mbco"),
                  error = function(e) e)
    if (inherits(r, "error")) {
      out$error <- conditionMessage(r)
    } else {
      mp <- map_mbco(r)
      out$original_fit_status <- "ok"
      out$estimate <- r@estimate
      out$lower <- r@ci_lower
      out$upper <- r@ci_upper
      out$converged <- isTRUE(r@converged)
      out$interval_status <- mp$interval
      out$lower_status <- mp$lower_status
      out$upper_status <- mp$upper_status
    }
  } else {
    cfg <- BOOT_CONFIGS[[cell$k]]
    d <- gen_boot(cell$k, cell$N, fixtures)
    out$truth <- truth_boot(cell$k, d, fixtures)
    fy <- if (cfg$z) Y ~ X + M + C + Z else Y ~ X + M + C
    fam <- function(f) if (f == "gaussian") stats::gaussian() else stats::binomial()
    r <- tryCatch(suppressWarnings(pmed(
      fy, formula_m = M ~ X + C, data = d, treatment = "X", mediator = "M",
      family_m = fam(cfg$m), family_y = fam(cfg$y), method = "model_bootstrap",
      n_boot = plan$n_boot, n_sim = plan$n_sim, seed = NULL
    )), error = function(e) e)
    if (inherits(r, "error")) {
      msg <- conditionMessage(r)
      out$error <- msg
      if (grepl("separation", msg)) out$original_fit_status <- "separated"
    } else {
      dg <- r@diagnostics
      out$original_fit_status <- "ok"
      out$estimate <- r@estimate
      out$lower <- r@ci_lower
      out$upper <- r@ci_upper
      out$converged <- isTRUE(r@converged)
      out$interval_status <- dg$interval_status
      out$lower_status <- unname(dg$endpoint_status[1])
      out$upper_status <- unname(dg$endpoint_status[2])
      out$n_valid <- dg$n_valid
      out$fail_separated <- dg$failures[["separated"]]
      out$fail_other <- sum(dg$failures) - dg$failures[["separated"]]
    }
  }
  out$seconds <- proc.time()[["elapsed"]] - t0
  as.data.frame(out, stringsAsFactors = FALSE)
}

res <- do.call(rbind, Map(one_rep, plan$reps, streams))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, file.path(outdir, sprintf("%s_task%04d.csv", mode, task)),
                 row.names = FALSE)
cat(sprintf("%s task %d: %s reps %d-%d, %.1f s total\n", mode, task, cell$cell_id,
            min(plan$reps), max(plan$reps), sum(res$seconds)))
