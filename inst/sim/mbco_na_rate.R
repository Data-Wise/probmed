# MBCO interval status rate and runtime (repair task T3).
#
# The verified inversion returns a missing interval when any profile point on
# its scan fails, in either sign basin. This measures how often that happens,
# and what each fit costs, before the repair is merged.
#
# Cells: 4 effect patterns x n in {100, 500} x {no covariate, one covariate};
# 200 seeded datasets per cell. Gaussian M and Y, continuous X, delta = 1.
#
# Usage (from the package root):
#   Rscript inst/sim/mbco_na_rate.R [reps] [cores]
# Output: inst/sim/results/mbco_na_rate.csv (one row per cell)

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 200L
cores <- if (length(args) >= 2) as.integer(args[2]) else
  max(1L, parallel::detectCores() - 2L)

suppressMessages(pkgload::load_all(".", quiet = TRUE))

effects <- list(
  null_a0 = c(a = 0, b = 0.4),
  null_b0 = c(a = 0.5, b = 0),
  weak = c(a = 0.15, b = 0.15),
  strong = c(a = 0.5, b = 0.4)
)
grid <- expand.grid(
  effect = names(effects), n = c(100L, 500L), cov = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

one_fit <- function(effect, n, cov, seed) {
  set.seed(seed)
  ab <- effects[[effect]]
  X <- stats::rnorm(n)
  C1 <- stats::rnorm(n)
  gam <- if (cov) 0.3 else 0
  M <- ab[["a"]] * X + gam * C1 + stats::rnorm(n)
  Y <- 0.2 * X + ab[["b"]] * M + gam * C1 + stats::rnorm(n)
  d <- data.frame(X = X, M = M, Y = Y, C1 = C1)
  fy <- if (cov) Y ~ X + M + C1 else Y ~ X + M
  fm <- if (cov) M ~ X + C1 else M ~ X
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(
    pmed(fy, formula_m = fm, data = d, treatment = "X", mediator = "M",
      method = "mbco"),
    error = function(e) NULL
  )
  el <- proc.time()[["elapsed"]] - t0
  if (is.null(r)) {
    return(c(pmed = "error", ie = "error", secs = el))
  }
  c(
    pmed = r@diagnostics$pmed$status %||% "none",
    ie = r@diagnostics$ie$status %||% "none",
    secs = el
  )
}
`%||%` <- function(x, y) if (is.null(x)) y else x

rows <- lapply(seq_len(nrow(grid)), function(g) {
  cell <- grid[g, ]
  seeds <- 1e5 * g + seq_len(reps)
  out <- parallel::mclapply(seeds, function(s)
    one_fit(cell$effect, cell$n, cell$cov, s), mc.cores = cores)
  out <- do.call(rbind, out)
  share <- function(col, lev) mean(out[, col] == lev)
  secs <- as.numeric(out[, "secs"])
  data.frame(
    effect = cell$effect, n = cell$n, cov = cell$cov, reps = reps,
    pmed_ok = share("pmed", "ok"), pmed_truncated = share("pmed", "truncated"),
    pmed_disconnected = share("pmed", "disconnected"),
    pmed_failed = share("pmed", "failed"), pmed_error = share("pmed", "error"),
    ie_ok = share("ie", "ok"), ie_truncated = share("ie", "truncated"),
    ie_disconnected = share("ie", "disconnected"),
    ie_failed = share("ie", "failed"),
    secs_median = stats::median(secs),
    secs_p95 = stats::quantile(secs, 0.95, names = FALSE)
  )
})
res <- do.call(rbind, rows)
dir.create("inst/sim/results", showWarnings = FALSE)
utils::write.csv(res, "inst/sim/results/mbco_na_rate.csv", row.names = FALSE)
print(res, digits = 3)
