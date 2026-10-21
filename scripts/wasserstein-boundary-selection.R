# Research only: Gaussian-experiment uncertainty is not a uniform guarantee.
source("scripts/wasserstein-joint-contrast.R")
ws_point <- function(x, v, lambda) {
  weak <- sweep(abs(x), 2, sqrt(diag(v)), "/") <= lambda
  p <- sqrt(abs(x[, 1])) / rowSums(sqrt(abs(x)))
  den <- sweep(x^2, 2, 4 * diag(v), "+")
  r <- sweep(1 / den, 2, diag(v), "*")
  cross <- v[1, 2] * x[, 1] * x[, 2] / (den[, 1] * den[, 2])
  shift <- (-(1 + 2 * p) * r[, 1] + (3 - 2 * p) * r[, 2]) / 8 -
    (1 - 2 * p) * cross / 4
  result <- stats::plogis(stats::qlogis(p) - shift)
  result[weak[, 1] & !weak[, 2]] <- 0
  result[!weak[, 1] & weak[, 2]] <- 1
  result
}
ws_interval <- function(d, v) {
  wj_interval(d, v, level = stats::pchisq(stats::qchisq(.95, 1), 2))
}
ws_fit <- function(d, v, lambda, z) {
  weak <- abs(d) / sqrt(diag(v)) <= lambda
  center <- d
  if (sum(weak) == 1L) {
    j <- which(weak)
    center <- d - v[, j] * d[j] / v[j, j]
  }
  draws <- sweep(z %*% chol(v), 2, center, "+")
  list(
    estimate = ws_point(matrix(d, 1), v, lambda),
    se = stats::sd(ws_point(draws, v, lambda)), center = center,
    selected = if (sum(weak) == 1L) which(weak) else 0L,
    ambiguous = all(weak)
  )
}
