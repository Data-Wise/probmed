source("inst/sim/i6_width_investigation.R", local = TRUE)

fixture <- data.frame(
  scenario = rep("s6", 4), N = rep(1000L, 4), rep = 1:4,
  boot_attempted = rep(2000L, 4), boot_failed_d = rep(0L, 4),
  ci_lower_d = c(0.4, 0.4, 0.4, 0.4),
  ci_upper_d = c(0.6, 0.6, 0.6, 0.6),
  ci_lower_e = c(0.45, 0.45, 0.45, 0.45),
  ci_upper_e = c(0.55, 0.55, 0.55, 0.55),
  exact_median = rep(0.5, 4)
)
i6_validate_parts(fixture, expected_reps = 4L, scenarios = "s6", sample_sizes = 1000L)
set.seed(1)
ans <- i6_analyze_cell(fixture, 0.5, n_boot = 50L)
stopifnot(abs(ans$observed_ratio - 2) < 1e-12,
          abs(ans$ratio_ci_lo - 2) < 1e-12,
          abs(ans$predicted_ratio - sqrt(1 + 0.25 /
            (200000 * (0.1 / (2 * stats::qnorm(0.975)))^2))) < 1e-12,
          ans$d_wider == 4L, ans$d_only_covered == 0L,
          abs(ans$mean_lower_shift + 0.05) < 1e-12,
          abs(ans$mean_upper_shift - 0.05) < 1e-12,
          abs(ans$mean_midpoint_shift) < 1e-12)

fails <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
stopifnot(fails(i6_validate_parts(rbind(fixture, fixture[1, ]),
                                 expected_reps = 4L, scenarios = "s6",
                                 sample_sizes = 1000L)))
bad <- fixture; bad$boot_failed_d[1] <- 1L
stopifnot(fails(i6_validate_parts(bad, expected_reps = 4L,
                                 scenarios = "s6", sample_sizes = 1000L)))
shifted <- fixture; shifted$ci_lower_d <- shifted$ci_lower_d + 0.01
set.seed(1)
stopifnot(i6_analyze_cell(shifted, 0.5, n_boot = 50L)$mean_midpoint_shift > 0)
cat("I6 width diagnostic fixtures: PASS\n")
