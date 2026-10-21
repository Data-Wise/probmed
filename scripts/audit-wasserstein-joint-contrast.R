# Independently recompute the saved frozen diagnostics, without loading study code.
root <- "docs/reports/data/wasserstein-joint-2026-10-20"
d <- utils::read.csv(file.path(root, "replicates.csv"), colClasses = c(error = "character"))
s <- utils::read.csv(file.path(root, "summary.csv"))
stopifnot(nrow(d) == 34000L, length(unique(d$seed)) == 14000L, all(d$error == ""))
for (previous in c("wasserstein-observed-2026-10-20", "wasserstein-boundary-2026-10-20")) {
  old <- utils::read.csv(file.path("docs/reports/data", previous, "replicates.csv"))
  stopifnot(length(intersect(d$seed, old$seed)) == 0L)
}
for (i in seq_len(nrow(s))) {
  z <- s[i, ]
  x <- d[d$cell == z$cell & d$method == z$method, ]
  stopifnot(nrow(x) == 2000L, length(unique(x$seed)) == 2000L)
  bias <- mean(x$estimate) - unique(x$truth)
  sd <- stats::sd(x$estimate)
  coverage <- sum(x$truth >= x$lower & x$truth <= x$upper) / nrow(x)
  ratio <- mean(x$se) / sd
  expected <- c(
    bias, sd / sqrt(nrow(x)), sd, abs(bias) / sd, ratio,
    coverage, sqrt(coverage * (1 - coverage) / nrow(x)), mean(x$upper - x$lower)
  )
  actual <- unlist(z[c(
    "bias", "bias_mcse", "empirical_sd", "standardized_bias",
    "se_ratio", "coverage", "coverage_mcse", "mean_width"
  )], use.names = FALSE)
  stopifnot(
    max(abs(expected - actual)) < 1e-11,
    z$pass == (abs(bias) / sd < 0.10 && ratio >= 0.90 && ratio <= 1.10 &&
                 coverage >= 0.92 && coverage <= 0.98)
  )
  if (z$method == "raw_projected") {
    empirical <- stats::cov(x[c("d1", "d2")])
    estimated <- matrix(c(mean(x$v11), mean(x$v12), mean(x$v12), mean(x$v22)), 2)
    eigenvalues <- eigen(solve(empirical, estimated), only.values = TRUE)$values
    delta <- cbind(x$d1 - x$true_d1, x$d2 - x$true_d2)
    quadratic <- (x$v22 * delta[, 1]^2 - 2 * x$v12 * delta[, 1] * delta[, 2] +
                    x$v11 * delta[, 2]^2) / (x$v11 * x$v22 - x$v12^2)
    covered <- quadratic <= stats::qchisq(0.95, 2)
    stopifnot(
      all(covered == x$joint_covered),
      all(!covered | (x$lower <= x$truth & x$truth <= x$upper))
    )
    cat(
      "Cell", z$cell, "covariance eigenvalues", eigenvalues,
      "width MCSE", stats::sd(x$upper - x$lower) / sqrt(nrow(x)), "\n"
    )
  }
}
cat("Audited all 17 summary rows, 14000 disjoint fresh seeds and joint containment.\n")
