# Research only: binary joint-contrast ellipse projection, not a package API.
source("scripts/wasserstein-boundary-candidates.R")
wj_interval <- function(d, v, level = 0.95) {
  if (length(d) != 2L || any(!is.finite(d)) ||
        !is.matrix(v) || !identical(dim(v), c(2L, 2L)) || any(!is.finite(v))) {
    stop("Two finite contrasts and a finite 2 by 2 covariance required")
  }
  if (!isTRUE(all.equal(v, t(v), tolerance = 1e-12))) stop("Covariance must be symmetric")
  if (min(eigen(v, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("Covariance must be positive definite")
  }
  stopifnot(length(level) == 1L, is.finite(level), level > 0, level < 1)
  m <- solve(v)
  critical <- stats::qchisq(level, 2)
  origin <- drop(crossprod(d, m %*% d))
  if (origin <= critical) {
    return(c(0, 1))
  }
  # Opposite rays form a line; squaring the inner product includes both signs.
  distance <- function(p) {
    rays <- cbind(c(p^2, (1 - p)^2), c(p^2, -(1 - p)^2))
    inner <- drop(crossprod(rays, m %*% d))
    norms <- colSums(rays * (m %*% rays))
    max(0, origin - max(inner^2 / norms)) - critical
  }
  center <- wb_ratio(d)
  lower <- if (d[1]^2 / v[1, 1] <= critical) {
    0
  } else {
    stats::uniroot(distance, c(0, center), tol = 1e-11)$root
  }
  upper <- if (d[2]^2 / v[2, 2] <= critical) {
    1
  } else {
    stats::uniroot(distance, c(center, 1), tol = 1e-11)$root
  }
  c(lower, upper)
}
