# Independent saved-data audit; no study helper is sourced.
root <- "docs/reports/data/wasserstein-selection-2026-10-20"
d <- utils::read.csv(file.path(root, "replicates.csv"), colClasses = c(error = "character"))
s <- utils::read.csv(file.path(root, "summary.csv"))
stopifnot(
  nrow(d) == 44000L, length(unique(d$seed)) == 14000L,
  all(d$error == ""), all(d$seed == 20561020L + d$cell * 10000L + d$replicate)
)
for (previous in c("observed", "boundary", "joint")) {
  old_root <- file.path("docs/reports/data", paste0("wasserstein-", previous, "-2026-10-20"))
  old <- utils::read.csv(file.path(old_root, "replicates.csv"))
  stopifnot(length(intersect(d$seed, old$seed)) == 0L)
}
for (i in seq_len(nrow(s))) {
  z <- s[i, ]
  x <- d[d$cell == z$cell & d$method == z$method, ]
  stopifnot(nrow(x) == 2000L, length(unique(x$seed)) == 2000L)
  bias <- mean(x$estimate) - unique(x$truth)
  sd <- stats::sd(x$estimate)
  coverage <- mean(x$truth >= x$lower & x$truth <= x$upper)
  ratio <- mean(x$se) / sd
  expected <- c(
    bias, sd / sqrt(nrow(x)), sd, abs(bias) / sd, ratio,
    coverage, sqrt(coverage * (1 - coverage) / nrow(x)), mean(x$upper - x$lower),
    mean(x$selected != 0), mean(x$ambiguous), mean(x$lower <= x$estimate & x$estimate <= x$upper)
  )
  actual <- unlist(z[c(
    "bias", "bias_mcse", "empirical_sd", "standardized_bias", "se_ratio",
    "coverage", "coverage_mcse", "mean_width", "selected_rate", "ambiguous_rate", "point_in_interval"
  )], use.names = FALSE)
  stopifnot(
    max(abs(expected - actual)) < 1e-11, z$failed == 0L,
    z$pass == (sd > 0 && abs(bias) / sd < .10 && ratio >= .90 && ratio <= 1.10 && coverage >= .92 && coverage <= .98)
  )
  if (grepl("fieller", z$method)) {
    # Independent linear-null test at truth: no geometric projection code.
    p <- x$truth
    a <- (1 - p)^2
    b <- -p^2
    positive <- (a * x$d1 + b * x$d2)^2 / (a^2 * x$v11 + 2 * a * b * x$v12 + b^2 * x$v22)
    negative <- (-a * x$d1 + b * x$d2)^2 / (a^2 * x$v11 - 2 * a * b * x$v12 + b^2 * x$v22)
    accepted <- pmin(positive, negative) <= stats::qchisq(.95, 1)
    stopifnot(all(accepted == (x$lower <= p & p <= x$upper)))
  }
  if (grepl("select", z$method)) {
    lambda <- stats::qnorm(if (z$method == "select95_fieller") .975 else .9975)
    weak1 <- abs(x$d1) / sqrt(x$v11) <= lambda
    weak2 <- abs(x$d2) / sqrt(x$v22) <= lambda
    selected <- ifelse(weak1 & !weak2, 1L, ifelse(!weak1 & weak2, 2L, 0L))
    stopifnot(
      all(selected == x$selected), all((weak1 & weak2) == x$ambiguous),
      all(x$estimate[selected == 1L] == 0), all(x$estimate[selected == 2L] == 1)
    )
    j <- selected == 1L
    stopifnot(
      all(abs(x$center1[j]) < 1e-12),
      all(abs(x$center2[j] - (x$d2[j] - x$v12[j] * x$d1[j] / x$v11[j])) < 1e-12)
    )
    j <- selected == 2L
    stopifnot(
      all(abs(x$center2[j]) < 1e-12),
      all(abs(x$center1[j] - (x$d1[j] - x$v12[j] * x$d2[j] / x$v22[j])) < 1e-12)
    )
  }
}
cat("Audited 22 summaries, 44000 rows, 14000 disjoint seeds, Fieller tests and selection centers.\n")
