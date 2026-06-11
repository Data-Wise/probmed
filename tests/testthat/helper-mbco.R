# Independent MBCO P_med interval: full-dimensional constrained fit (no profiling).
# Mirrors the standalone prototype's Approach-A and supports one covariate C1.
# Model: M ~ 1 + X (+ C1); Y ~ 1 + X + M (+ C1). delta = 1.

.oracle_ll_free <- function(data, covs) {
  fm <- stats::lm(stats::reformulate(c("X", covs), "M"), data)
  fy <- stats::lm(stats::reformulate(c("X", "M", covs), "Y"), data)
  as.numeric(stats::logLik(fm)) + as.numeric(stats::logLik(fy))
}

.oracle_ll_constrained <- function(data, qstar, a_sign, covs) {
  n <- nrow(data)
  Xc <- if (length(covs)) as.matrix(data[, covs, drop = FALSE]) else NULL
  ll_joint <- function(i_m, a, gm, Vm, i_y, cp, b, gy, Vy) {
    mu_m <- i_m + a * data$X + (if (is.null(Xc)) 0 else Xc %*% gm)
    mu_y <- i_y + cp * data$X + b * data$M + (if (is.null(Xc)) 0 else Xc %*% gy)
    sum(stats::dnorm(data$M, mu_m, sqrt(Vm), log = TRUE)) +
      sum(stats::dnorm(data$Y, mu_y, sqrt(Vy), log = TRUE))
  }
  if (abs(qstar) < 1e-8) {
    fm <- stats::lm(stats::reformulate(c("X", covs), "M"), data)
    fy <- stats::lm(stats::reformulate(c("X", covs), "Y"), data)
    return(as.numeric(stats::logLik(fm)) + as.numeric(stats::logLik(fy)))
  }
  fm <- stats::lm(stats::reformulate(c("X", covs), "M"), data)
  fy <- stats::lm(stats::reformulate(c("X", "M", covs), "Y"), data)
  k <- length(covs)
  sb <- sign(qstar) * a_sign
  st <- c(
    stats::coef(fm)[1], stats::coef(fm)["X"],
    if (k) stats::coef(fm)[covs] else numeric(0), log(stats::sigma(fm)^2),
    stats::coef(fy)[1], stats::coef(fy)["X"],
    log(abs(stats::coef(fy)["M"]) + 1e-3),
    if (k) stats::coef(fy)[covs] else numeric(0)
  )
  nll <- function(p) {
    i <- 1
    i_m <- p[i]
    i <- i + 1
    a <- p[i]
    i <- i + 1
    gm <- if (k) p[i:(i + k - 1)] else numeric(0)
    i <- i + k
    Vm <- exp(p[i])
    i <- i + 1
    i_y <- p[i]
    i <- i + 1
    cp <- p[i]
    i <- i + 1
    b <- sb * exp(p[i])
    i <- i + 1
    gy <- if (k) p[i:(i + k - 1)] else numeric(0)
    Vy <- a^2 * b^2 / (2 * qstar^2) - b^2 * Vm # delta fixed at 1
    if (!is.finite(Vy) || Vy <= 1e-8) {
      return(1e10)
    }
    -ll_joint(i_m, a, gm, Vm, i_y, cp, b, gy, Vy)
  }
  o <- try(
    stats::optim(unname(st), nll,
      method = "Nelder-Mead",
      control = list(maxit = 6000, reltol = 1e-11)
    ),
    silent = TRUE
  )
  if (inherits(o, "try-error") || o$value >= 1e9) {
    return(NA_real_)
  }
  -o$value
}

oracle_mbco_ci <- function(data, level = 0.95, covs = character(0)) {
  fm <- stats::lm(stats::reformulate(c("X", covs), "M"), data)
  fy <- stats::lm(stats::reformulate(c("X", "M", covs), "Y"), data)
  a <- stats::coef(fm)["X"]
  b <- stats::coef(fy)["M"]
  Vm <- stats::sigma(fm)^2
  Vy <- stats::sigma(fy)^2
  est <- as.numeric(stats::pnorm(a * b / sqrt(2 * (b^2 * Vm + Vy))))
  a_sign <- sign(as.numeric(a))
  if (a_sign == 0) a_sign <- 1
  ll1 <- .oracle_ll_free(data, covs)
  crit <- stats::qchisq(level, df = 1)
  excess <- function(ps) {
    ps <- min(max(ps, 1e-5), 1 - 1e-5)
    ll0 <- .oracle_ll_constrained(data, stats::qnorm(ps), a_sign, covs)
    if (is.na(ll0)) {
      return(crit + 1e3)
    }
    -2 * (ll0 - ll1) - crit
  }
  endpoint <- function(dir) {
    inner <- est
    for (k in seq_len(200)) {
      outer <- min(max(est + dir * 0.01 * k, 1e-4), 1 - 1e-4)
      if (excess(outer) > 0) {
        return(tryCatch(stats::uniroot(excess, sort(c(inner, outer)), tol = 1e-4)$root,
          error = function(e) NA_real_
        ))
      }
      inner <- outer
      if (outer <= 1e-4 || outer >= 1 - 1e-4) {
        return(outer)
      }
    }
    outer
  }
  c(estimate = est, lower = endpoint(-1), upper = endpoint(+1))
}
