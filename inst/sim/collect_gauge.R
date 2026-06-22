#!/usr/bin/env Rscript
# Collect the per-cell rds written by the gauge coverage SLURM array into one
# summary table. Run in the sim working dir (where out/ lives):
#   Rscript collect_gauge.R           # prints the table
#   Rscript collect_gauge.R out.csv   # also writes a CSV
# Reports how many of the 24 cells have landed so a partial run is obvious.
args <- commandArgs(trailingOnly = TRUE)
outdir <- Sys.getenv("GAUGE_SIM_OUT", "out")
f <- list.files(outdir, pattern = "^cell_.*\\.rds$", full.names = TRUE)
if (!length(f)) stop("no cell_*.rds in ", outdir, " yet")
res <- do.call(rbind, lapply(f, readRDS))
res <- res[order(res$se_method, res$tau, res$n), ]
cat(sprintf("%d / 24 cells present\n", nrow(res)))
print(res[, c("se_method","n","tau","true_W","cov_W","mcse_W","cov_P","mcse_P")],
      row.names = FALSE)
if (length(args)) { write.csv(res, args[1], row.names = FALSE); cat("wrote", args[1], "\n") }
