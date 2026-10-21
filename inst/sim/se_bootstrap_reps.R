#!/usr/bin/env Rscript
# Fallback SE candidate: nonparametric bootstrap OF the reps=R point estimator
# (resample rows, refit the whole cross-fit R times per resample). Cell 1.
# Cost B*R cross-fits per dataset, so reduced size: nrep=60, B=50, R=10.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
src <- readLines("inst/sim/phi_decomposition.R")
i1 <- grep("^gen <- function", src); i2 <- grep("^## ---- cells", src)
# eval(parse()) sources function DEFINITIONS from a committed repo script, not
# external input -- the standard way to reuse a CLI script's functions.
eval(parse(text = src[i1:(i2 - 1)]))

args <- commandArgs(trailingOnly = TRUE)
nrep <- if (length(args) >= 1) as.integer(args[1]) else 60L
B    <- if (length(args) >= 2) as.integer(args[2]) else 50L
n <- 800L; K <- 5L; R <- 10L; s <- 1; tau <- 0; binY <- FALSE
tr <- truth(s, tau, binY)
Wof <- function(phi) { t <- colMeans(phi); OE <- t["11"] - t["00"]
  unname((OE - (t["10"] - t["00"]) - (t["01"] - t["00"])) / OE) }
W_reps_avg <- function(d, seed) {              # the reps=R point estimator, as shipped
  phis <- lapply(seq_len(R), function(p) { set.seed(seed + p); .corner_fit(d, K, binY, "C")$phi })
  Wof(Reduce(`+`, phis) / R) }

rows <- vector("list", nrep); t0 <- proc.time()[["elapsed"]]
for (r in seq_len(nrep)) {
  seed <- 1000000L + r; set.seed(seed); d <- gen(n, s, tau, binY)
  W <- W_reps_avg(d, seed)
  set.seed(seed * 7L)
  Wb <- vapply(seq_len(B), function(b) {
    db <- d[sample.int(n, n, replace = TRUE), , drop = FALSE]
    W_reps_avg(db, seed * 1000L + b) }, numeric(1))
  so <- gauge_stats(phi_oracle(d, s, tau, binY))
  rows[[r]] <- c(W = W, se_boot = sd(Wb), lo_pct = unname(quantile(Wb, .025)), hi_pct = unname(quantile(Wb, .975)),
                 W_or = unname(so["W"]), se_or = unname(so["seW"]))
  if (r %% 10 == 0) message(sprintf("  %d/%d datasets, %.0fs", r, nrep, proc.time()[["elapsed"]] - t0))
}
x <- as.data.frame(do.call(rbind, rows)); emp <- sd(x$W); z <- qnorm(0.975)
cat(sprintf("bootstrap of reps=%d estimator: cell 1, nrep=%d, B=%d, %.0fs | empSD(W_avg)=%.4f sd(W_or)=%.4f\n",
            R, nrep, B, proc.time()[["elapsed"]] - t0, emp, sd(x$W_or)))
res <- rbind(
  boot_sd_wald = c(med_ratio = median(x$se_boot)/emp, mean_ratio = mean(x$se_boot)/emp, cv = sd(x$se_boot)/mean(x$se_boot),
                   cov = mean(abs(x$W - tr$W) <= z * x$se_boot), q05 = quantile(x$se_boot,.05)/emp, q95 = quantile(x$se_boot,.95)/emp),
  boot_percentile = c(med_ratio = NA, mean_ratio = NA, cv = NA, cov = mean(tr$W >= x$lo_pct & tr$W <= x$hi_pct),
                      q05 = NA, q95 = NA),
  oracle = c(med_ratio = median(x$se_or)/emp, mean_ratio = mean(x$se_or)/emp, cv = sd(x$se_or)/mean(x$se_or),
             cov = mean(abs(x$W_or - tr$W) <= z * x$se_or), q05 = quantile(x$se_or,.05)/emp, q95 = quantile(x$se_or,.95)/emp))
print(round(res, 3))
saveRDS(list(rows = x, res = res), file.path(Sys.getenv("GAUGE_SIM_OUT", "gauge_sim_out"), "se_bootstrap_reps_cell1.rds"))
