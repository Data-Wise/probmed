# T20 positive control: is model_bootstrap's large-N overcoverage Monte Carlo
# noise? Reruns boot-1-N1000 (Gaussian M and Y) on the main run's own streams,
# so every dataset is identical, with the exact P_med inside each bootstrap draw
# in place of the n_sim simulation. If coverage drops from the main run's level
# toward 0.95 and the interval narrows, the inflation is MC noise.
#
# Usage (from the package root):
#   Rscript inst/sim/calibration/calib_t20_exact.R <main_calibration.csv.gz> <out.csv> [reps] [cores]

args <- commandArgs(trailingOnly = TRUE)
main_csv <- args[1]
out_csv <- args[2]
reps <- if (length(args) >= 3) seq_len(as.integer(args[3])) else 1:200
cores <- if (length(args) >= 4) as.integer(args[4]) else max(1L, parallel::detectCores() - 2L)

suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("inst/sim/calibration/calib_design.R")
fixtures <- load_calib_fixtures("inst/sim/calibration/calib_fixtures.dcf")
cells <- calib_cells()
ci <- which(cells$cell_id == "boot-1-N1000")
cell <- cells[ci, ]

# Exact P_med for Gaussian M and Y: Y_t - Y_c ~ N(a b (x - x*), 2 (b^2 sm^2 + sy^2)).
# Covariates enter additively and cancel, so standardization is irrelevant here.
exact_core <- function(a, b, c_prime = 0, x_ref = 0, x_value = 1, sigma_y = 1,
                       sigma_m = 1, n_sim = 10000, family_y = NULL,
                       family_m = NULL, ...) {
  stopifnot(.is_gaussian(family_m), .is_gaussian(family_y))
  stats::pnorm(a * b * (x_value - x_ref) / sqrt(2 * (b^2 * sigma_m^2 + sigma_y^2)))
}

# Known-answer check of the replacement against the simulation it replaces.
set.seed(1)
g <- stats::gaussian()
mc <- .pmed_core_simple(a = 0.5, b = 0.4, c_prime = 0.3, x_ref = 0, x_value = 1, sigma_y = 1, sigma_m = 1, n_sim = 4e6,
                        family_y = g, family_m = g)
ex <- exact_core(a = 0.5, b = 0.4, family_y = g, family_m = g)
cat(sprintf("known-answer: MC %.4f vs exact %.4f\n", mc, ex))
stopifnot(abs(mc - ex) < 0.002)

ns <- asNamespace("probmed")
unlockBinding(".pmed_core_simple", ns)
assign(".pmed_core_simple", exact_core, envir = ns)

gidx <- (ci - 1L) * CALIB_REPS + reps
streams <- calib_streams(min(gidx), max(gidx))
one <- function(k) {
  assign(".Random.seed", streams[[k]], envir = globalenv())
  d <- gen_boot(cell$k, cell$N, fixtures)
  truth <- truth_boot(cell$k, d, fixtures)
  r <- tryCatch(suppressWarnings(pmed(
    Y ~ X + M + C, formula_m = M ~ X + C, data = d, treatment = "X",
    mediator = "M", family_m = stats::gaussian(), family_y = stats::gaussian(),
    method = "model_bootstrap", n_boot = BOOT_N_BOOT, n_sim = BOOT_N_SIM, seed = NULL
  )), error = function(e) NULL)
  if (is.null(r)) return(data.frame(replication = reps[k], truth = truth,
                                    estimate = NA, lower = NA, upper = NA))
  data.frame(replication = reps[k], truth = truth, estimate = r@estimate,
             lower = r@ci_lower, upper = r@ci_upper)
}
res <- do.call(rbind, parallel::mclapply(seq_along(reps), one, mc.cores = cores))
utils::write.csv(res, out_csv, row.names = FALSE)

m <- utils::read.csv(main_csv)
m <- m[m$cell_id == "boot-1-N1000" & m$replication %in% reps, ]
k <- merge(m, res, by = "replication", suffixes = c(".mc", ".exact"))
stopifnot(isTRUE(all.equal(k$truth.mc, k$truth.exact)))
cov <- function(l, u, t) !is.na(l) & l <= t & u >= t
cm <- cov(k$lower.mc, k$upper.mc, k$truth.mc)
ce <- cov(k$lower.exact, k$upper.exact, k$truth.exact)
bt <- function(x) sprintf("%.3f [%.3f, %.3f]", mean(x),
                          stats::binom.test(sum(x), length(x))$conf.int[1],
                          stats::binom.test(sum(x), length(x))$conf.int[2])
cat(sprintf("reps %d | coverage MC %s | exact %s\n", nrow(k), bt(cm), bt(ce)))
cat(sprintf("median width MC %.4f | exact %.4f | ratio %.3f\n",
            stats::median(k$upper.mc - k$lower.mc, na.rm = TRUE),
            stats::median(k$upper.exact - k$lower.exact, na.rm = TRUE),
            stats::median(k$upper.exact - k$lower.exact, na.rm = TRUE) /
              stats::median(k$upper.mc - k$lower.mc, na.rm = TRUE)))
cat(sprintf("point-estimate sd(MC - exact) %.4f | sd(exact estimate) %.4f\n",
            stats::sd(k$estimate.mc - k$estimate.exact, na.rm = TRUE),
            stats::sd(k$estimate.exact, na.rm = TRUE)))
cat(sprintf("discordant coverage: MC-only %d, exact-only %d\n", sum(cm & !ce), sum(ce & !cm)))
