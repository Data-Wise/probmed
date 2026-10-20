# Diagnostic only: this script does not alter the frozen I6 gate.
# Usage: Rscript inst/sim/i6_width_investigation.R <parts_dir> <output.csv>

i6_gaussian_truth <- c(
  s1 = stats::pnorm(0.4 * 0.5 / sqrt(2 * (0.16 + 1))),
  s6 = 0.5,
  s7 = 0.5,
  s8 = stats::pnorm(1.5 * 2.5 / sqrt(2 * (1.5^2 + 2)))
)

i6_validate_parts <- function(x, expected_reps = 2000L,
                              scenarios = c("s1", "s3", "s4", "s6", "s7", "s8"),
                              sample_sizes = c(50L, 100L, 250L, 500L, 1000L)) {
  required <- c("scenario", "N", "rep", "boot_attempted", "boot_failed_d",
                "ci_lower_d", "ci_upper_d", "ci_lower_e", "ci_upper_e",
                "exact_median")
  if (!all(required %in% names(x))) stop("missing required I6 columns")
  if (anyNA(x[c("scenario", "N", "rep")]) ||
      anyDuplicated(x[c("scenario", "N", "rep")])) stop("missing or duplicate I6 key")
  if (!all(x$boot_attempted == 2000L)) stop("unexpected bootstrap count")
  grid <- expand.grid(scenario = scenarios, N = sample_sizes,
                      stringsAsFactors = FALSE)
  if (nrow(x) != nrow(grid) * expected_reps) stop("incomplete I6 grid")
  for (j in seq_len(nrow(grid))) {
    ids <- x$rep[x$scenario == grid$scenario[j] & x$N == grid$N[j]]
    if (!identical(sort(as.integer(ids)), seq_len(expected_reps)))
      stop("incomplete I6 replication IDs")
  }
  g <- x[x$scenario %in% names(i6_gaussian_truth), , drop = FALSE]
  if (anyNA(g$boot_failed_d) || any(g$boot_failed_d != 0L))
    stop("failed Gaussian d resamples")
  v <- g[c("ci_lower_d", "ci_upper_d", "ci_lower_e", "ci_upper_e", "exact_median")]
  if (!all(vapply(v, function(z) all(is.finite(z) & z >= 0 & z <= 1), logical(1))))
    stop("nonfinite or out-of-range Gaussian value")
  if (any(g$ci_lower_d > g$ci_upper_d) || any(g$ci_lower_e > g$ci_upper_e))
    stop("reversed Gaussian interval")
  invisible(g)
}

i6_analyze_cell <- function(cell, truth, k = 200000L, n_boot = 5000L) {
  n <- nrow(cell)
  wd <- cell$ci_upper_d - cell$ci_lower_d
  we <- cell$ci_upper_e - cell$ci_lower_e
  mean_e <- mean(we)
  if (!is.finite(mean_e) || mean_e <= 0) stop("zero or invalid exact width")
  ratio <- mean(wd) / mean_e
  paired <- replicate(n_boot, {
    ix <- sample.int(n, n, replace = TRUE)
    mean(wd[ix]) / mean(we[ix])
  })
  p_var <- mean(cell$exact_median * (1 - cell$exact_median))
  exact_sd_proxy <- mean_e / (2 * stats::qnorm(0.975))
  predicted <- sqrt(1 + p_var / (k * exact_sd_proxy^2))
  cover_d <- cell$ci_lower_d <= truth & truth <= cell$ci_upper_d
  cover_e <- cell$ci_lower_e <= truth & truth <= cell$ci_upper_e
  data.frame(
    n = n, exact_mean_width = mean_e, d_mean_width = mean(wd),
    observed_ratio = ratio,
    ratio_ci_lo = unname(stats::quantile(paired, 0.025)),
    ratio_ci_hi = unname(stats::quantile(paired, 0.975)),
    predicted_ratio = predicted, observed_minus_predicted = ratio - predicted,
    mean_lower_shift = mean(cell$ci_lower_d - cell$ci_lower_e),
    mean_upper_shift = mean(cell$ci_upper_d - cell$ci_upper_e),
    mean_midpoint_shift = mean((cell$ci_lower_d + cell$ci_upper_d -
                                  cell$ci_lower_e - cell$ci_upper_e) / 2),
    d_wider = sum(wd > we), d_only_covered = sum(cover_d & !cover_e),
    e_only_covered = sum(cover_e & !cover_d),
    coverage_d = mean(cover_d), coverage_e = mean(cover_e)
  )
}

i6_width_investigation <- function(parts_dir, output_csv, n_boot = 5000L,
                                   seed = 20261020L) {
  files <- sort(list.files(parts_dir, pattern = "^i6_.*\\.csv$", full.names = TRUE))
  if (length(files) != 300L) stop("expected 300 I6 part files")
  x <- do.call(rbind, lapply(files, utils::read.csv))
  g <- i6_validate_parts(x)
  set.seed(seed)
  groups <- split(g, list(g$scenario, g$N), drop = TRUE)
  out <- do.call(rbind, lapply(groups, function(cell) {
    ans <- i6_analyze_cell(cell, i6_gaussian_truth[[cell$scenario[1]]],
                           n_boot = n_boot)
    cbind(scenario = cell$scenario[1], N = cell$N[1], ans)
  }))
  rownames(out) <- NULL
  out <- out[order(out$scenario, out$N), , drop = FALSE]
  utils::write.csv(out, output_csv, row.names = FALSE)
  out
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("usage: Rscript inst/sim/i6_width_investigation.R <parts_dir> <output.csv>")
  print(i6_width_investigation(args[[1]], args[[2]]), digits = 6, row.names = FALSE)
}
