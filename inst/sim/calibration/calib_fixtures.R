# Freeze the calibration fixtures (task T15 prerequisite).
#
# Calibrates the Bernoulli outcome intercepts of bootstrap configurations 2, 3
# and 5 by deterministic quadrature + uniroot to their target marginal
# prevalences, records convergence evidence at two quadrature tolerances
# (agreement required within 1e-8), and writes the exact MBCO truths.
#
# Usage (from the package root):
#   Rscript inst/sim/calibration/calib_fixtures.R
# Output: inst/sim/calibration/calib_fixtures.dcf (hash it with
#   shasum -a 256 and record the digest in calib_manifest.dcf).

source("inst/sim/calibration/calib_design.R")

out <- list(base_seed = CALIB_BASE_SEED, reps_per_cell = CALIB_REPS,
            rng_kind = "L'Ecuyer-CMRG / Inversion / Rejection",
            n_boot = BOOT_N_BOOT, n_sim = BOOT_N_SIM)

for (cfg in BOOT_CONFIGS) {
  if (cfg$y != "binomial") next
  k <- cfg$id
  root <- stats::uniroot(function(iy) boot_prevalence(k, iy) - cfg$prevalence,
                         c(-10, 5), tol = 1e-14)
  iy <- root$root
  prev_hi <- boot_prevalence(k, iy, rel.tol = 1e-12)
  prev_lo <- boot_prevalence(k, iy, rel.tol = 1e-9)
  # Convergence evidence: the prevalence at the root, at two quadrature
  # tolerances, and the truth functional of a fixed probe dataset.
  set.seed(1)
  probe <- data.frame(C = stats::rbinom(200, 1, 0.5), Z = stats::rnorm(200))
  fx <- list(intercept_y = stats::setNames(list(iy), paste0("boot-", k)))
  t_hi <- truth_boot(k, probe, fx, rel.tol = 1e-12)
  t_lo <- truth_boot(k, probe, fx, rel.tol = 1e-9)
  stopifnot(abs(prev_hi - cfg$prevalence) < 1e-10,
            abs(prev_hi - prev_lo) < 1e-8,
            abs(t_hi - t_lo) < 1e-8)
  out[[sprintf("intercept_y_boot_%d", k)]] <- sprintf("%.15f", iy)
  out[[sprintf("prevalence_boot_%d", k)]] <- sprintf("%.15f", prev_hi)
  out[[sprintf("uniroot_iter_boot_%d", k)]] <- root$iter
  out[[sprintf("uniroot_estim_prec_boot_%d", k)]] <- signif(root$estim.prec, 3)
  out[[sprintf("quadrature_agreement_prev_boot_%d", k)]] <- signif(abs(prev_hi - prev_lo), 3)
  out[[sprintf("quadrature_agreement_truth_boot_%d", k)]] <- signif(abs(t_hi - t_lo), 3)
}
for (k in seq_along(MBCO_PAIRS)) {
  out[[sprintf("truth_mbco_%d", k)]] <- sprintf("%.15f", truth_mbco(k))
}
out$truth_boot_1 <- sprintf("%.15f", truth_boot(1L, NULL, NULL))

path <- "inst/sim/calibration/calib_fixtures.dcf"
write.dcf(as.data.frame(out, check.names = FALSE), path)
cat(readLines(path), sep = "\n")
