# Score candidate MBCO cutoff rules offline from the composite-null harness
# (Phase 8, T4). Usage: Rscript inst/sim/calibration/calib_mbco_null_score.R <partsdir> <out.csv>
# or source() it and call score_rules(d).
#
# Rules (prespecified in tasks/plan-mbco-null.md):
#   chi2      : c = qchisq(1 - alpha, 1) everywhere (status quo)
#   switch_crit / switch_logN : c = qchisq(1 - sqrt(alpha), 1) when both branch
#               statistics are below lambda (lambda = chi2 cutoff, or log N),
#               otherwise the chi2 cutoff
#   smooth    : interpolate linearly from the double-null cutoff at max(LR_a0, LR_b0) = 0
#               to the chi2 cutoff at max(...) >= chi2 cutoff
# Acceptance (prespecified): min coverage over the lattice >= 0.94 at each N,
# and coverage at (0, 0) <= 0.97.

alpha <- 0.05
cut_rule <- function(rule, lr_a0, lr_b0, N) {
  c1 <- stats::qchisq(1 - alpha, 1)
  c0 <- stats::qchisq(1 - sqrt(alpha), 1)
  m <- pmax(lr_a0, lr_b0)
  switch(rule,
    chi2 = rep(c1, length(m)),
    switch_crit = ifelse(m < c1, c0, c1),
    switch_logN = ifelse(m < log(N), c0, c1),
    smooth = c0 + (c1 - c0) * pmin(1, m / c1),
    stop("unknown rule ", rule))
}
RULES <- c("chi2", "switch_crit", "switch_logN", "smooth")

score_rules <- function(d) {
  pts <- split(d, list(d$h_a, d$h_b, d$N, d$cov), drop = TRUE)
  do.call(rbind, lapply(pts, function(z) {
    do.call(rbind, lapply(RULES, function(rule) {
      acc <- z$lr_true <= cut_rule(rule, z$lr_a0, z$lr_b0, z$N)
      ci <- stats::binom.test(sum(acc, na.rm = TRUE), sum(!is.na(acc)))$conf.int
      data.frame(rule, h_a = z$h_a[1], h_b = z$h_b[1], N = z$N[1], cov = z$cov[1],
                 reps = sum(!is.na(acc)), coverage = mean(acc, na.rm = TRUE),
                 ci_lower = ci[1], ci_upper = ci[2], stringsAsFactors = FALSE)
    }))
  }))
}

summarize_rules <- function(sc) {
  do.call(rbind, lapply(split(sc, list(sc$rule, sc$N), drop = TRUE), function(z) {
    o <- z[z$h_a == 0 & z$h_b == 0, ]
    data.frame(rule = z$rule[1], N = z$N[1], min_coverage = min(z$coverage),
               at_min = paste0("(", z$h_a[which.min(z$coverage)], ",", z$h_b[which.min(z$coverage)],
                               if (z$cov[which.min(z$coverage)]) ",cov" else "", ")"),
               origin_coverage = max(o$coverage),
               pass = min(z$coverage) >= 0.94 && max(o$coverage) <= 0.97)
  }))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  d <- do.call(rbind, lapply(list.files(args[1], "^null_task.*\\.csv$", full.names = TRUE),
                             utils::read.csv))
  sc <- score_rules(d)
  utils::write.csv(sc, args[2], row.names = FALSE)
  print(summarize_rules(sc), row.names = FALSE, digits = 4)
}
