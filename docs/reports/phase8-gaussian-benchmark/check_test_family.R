# Deterministic numerical checks of the ideal Gaussian product test family.
# Run with Rscript --vanilla. No random draws or saved replications are used.

chi1_tail <- function(cutoff) stats::pchisq(cutoff, 1, lower.tail = FALSE)
ray_tail <- function(cutoff) {
  q <- chi1_tail(cutoff)
  (q^2 + 2 * q + exp(-cutoff / 2)) / 4
}

distance_cutoff <- function(eta) {
  stats::uniroot(function(cutoff) ray_tail(cutoff) - (0.05 - eta),
                 c(0, 20), tol = 1e-13)$root
}

chi1 <- stats::qchisq(0.95, 1)
chi2 <- stats::qchisq(0.95, 2)
eta <- c(0.01, 0.04)
cutoffs <- c(chi1, vapply(eta, distance_cutoff, numeric(1)), chi2)
names(cutoffs) <- c("chi1", "eta_01", "eta_04", "chi2")

# Reflect the second coordinate before calling this function for a negative
# target. Each branch contributes a single interval to the vertical slice.
tube_slice <- function(x, radius, target) {
  if (target == 0) {
    if (abs(x) <= radius) return(matrix(c(-Inf, Inf), nrow = 1))
    return(matrix(c(-radius, radius), nrow = 1))
  }

  # With w = (s - x) / radius, the squared stationary-point equation
  # is w^2 * ((x + radius * w)^4 + target^2) = target^2.
  polynomial <- c(-target^2, 0, x^4 + target^2, 4 * x^3 * radius,
                  6 * x^2 * radius^2, 4 * x * radius^3, radius^4)
  roots <- polyroot(polynomial)
  w <- c(-1, 1, Re(roots[abs(Im(roots)) < 1e-9 & abs(Re(roots)) < 1]))
  s <- x + radius * w
  dy <- radius * sqrt(pmax(0, 1 - w^2))

  intervals <- matrix(numeric(), nrow = 0, ncol = 2)
  for (branch in c(-1, 1)) {
    keep <- branch * s > 0
    if (!any(keep)) next
    lo <- min(target / s[keep] - dy[keep])
    hi <- max(target / s[keep] + dy[keep])
    if (abs(x) <= radius && branch < 0) lo <- -Inf
    if (abs(x) <= radius && branch > 0) hi <- Inf
    intervals <- rbind(intervals, c(lo, hi))
  }
  if (nrow(intervals) == 2) {
    intervals <- intervals[order(intervals[, 1]), , drop = FALSE]
    if (intervals[1, 2] >= intervals[2, 1]) {
      intervals <- matrix(c(intervals[1, 1], max(intervals[, 2])), nrow = 1)
    }
  }
  intervals
}

slice_probability <- function(x, radius, target, v, inner_sq) {
  intervals <- tube_slice(x, radius, target)
  lo <- intervals[, 1]
  hi <- intervals[, 2]
  prob <- sum(stats::pnorm(hi - v) - stats::pnorm(lo - v))
  if (x^2 < inner_sq) {
    halfwidth <- sqrt(inner_sq - x^2)
    clipped_lo <- pmax(lo, -halfwidth)
    clipped_hi <- pmin(hi, halfwidth)
    prob <- prob - sum(pmax(0, stats::pnorm(clipped_hi - v) -
                              stats::pnorm(clipped_lo - v)))
  }
  prob
}

acceptance <- function(target, u, v, cutoff, inner_sq, tolerance, mesh) {
  if (target < 0) v <- -v
  radius <- sqrt(cutoff)
  inner <- sqrt(inner_sq)
  bounds <- c(u - 9, u + 9)
  knots <- c(-radius, -inner, 0, inner, radius, seq(bounds[1], bounds[2], by = mesh))
  knots <- sort(unique(c(bounds, knots[knots > bounds[1] & knots < bounds[2]])))
  density <- function(x) {
    stats::dnorm(x - u) * vapply(x, slice_probability, numeric(1), radius = radius,
                                target = abs(target), v = v, inner_sq = inner_sq)
  }
  parts <- lapply(seq_len(length(knots) - 1L), function(i) {
    stats::integrate(density, knots[i], knots[i + 1L], subdivisions = 1500L,
                     rel.tol = tolerance, abs.tol = tolerance)$value
  })
  sum(unlist(parts))
}

compare <- function(target, u, v, tolerance, mesh) {
  inner_sq <- c(0, stats::qchisq(eta, 2, ncp = 2 * abs(target)), 0)
  vapply(seq_along(cutoffs), function(i) {
    acceptance(target, u, v, cutoffs[i], inner_sq[i], tolerance, mesh)
  }, numeric(1))
}

cases <- data.frame(
  target = c(0, 0, 0, 0.25, 0.25, 1, -1, 1, 1, 1, 4, 1),
  u = c(0, 3, 12, 0.5, 2, 1, 1, 2, 0.5, 12, 2, 0),
  v = c(0, 0, 0, 0.5, 0.125, 1, -1, 0.5, 2, 1 / 12, 2, 0)
)
checked <- t(mapply(compare, cases$target, cases$u, cases$v,
                    MoreArgs = list(tolerance = 1e-11, mesh = 0.25)))
recheck <- t(mapply(compare, cases$target, cases$u, cases$v,
                    MoreArgs = list(tolerance = 1e-10, mesh = 0.41)))
colnames(checked) <- names(cutoffs)

axis_formula <- function(u) {
  radius <- sqrt(cutoffs)
  tail <- stats::pnorm(-radius - u) + stats::pnorm(u - radius)
  inner_sq <- c(0, stats::qchisq(eta, 2), 0)
  stopifnot(all(inner_sq < cutoffs))
  1 - tail * chi1_tail(cutoffs) - stats::pchisq(inner_sq, 2, ncp = u^2)
}
axis <- cases$target == 0
axis_expected <- t(vapply(cases$u[axis], axis_formula, numeric(4)))
partition_delta <- max(abs(checked - recheck))
axis_delta <- max(abs(checked[axis, ] - axis_expected))
on_null <- cases$u * cases$v == cases$target
stopifnot(partition_delta < 2e-7, axis_delta < 1e-10,
          max(abs(checked[6, ] - checked[7, ])) < 1e-10,
          max(abs(checked[8, ] - checked[9, ])) < 2e-8,
          all(checked[on_null, ] >= 0.95 - 1e-10),
          max(abs(ray_tail(cutoffs[2:3]) - (0.05 - eta))) < 1e-12)

print(cbind(cases, checked), digits = 9, row.names = FALSE)
cat("Maximum partition change:", format(partition_delta, digits = 4), "\n")
cat("Maximum axis formula discrepancy:", format(axis_delta, digits = 4), "\n")

ceiling_cutoff <- stats::uniroot(function(cutoff) {
  0.95 + ray_tail(cutoff) - chi1_tail(cutoff)^2 - 0.97
}, c(4, 10), tol = 1e-13)$root
ceiling_eta <- 0.05 - ray_tail(ceiling_cutoff)
pilot_cutoffs <- c(cutoffs[2:3], ceiling_cutoff)
pilot_eta <- c(eta, ceiling_eta)
print(data.frame(eta = pilot_eta, cutoff = pilot_cutoffs,
                 origin_coverage = 1 - pilot_eta - chi1_tail(pilot_cutoffs)^2,
                 strong_path_coverage = 1 - chi1_tail(pilot_cutoffs),
                 strong_path_width_ratio = sqrt(pilot_cutoffs / chi1)),
      digits = 9, row.names = FALSE)
