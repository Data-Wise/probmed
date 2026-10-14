# Independent full likelihood and unsquared probability constraint. No package
# likelihood, profiling, optimization, or constraint helpers are called.
repair_full_oracle <- function(d, target, delta = 1, scales = c(0.5, 1, 2), estimand = "pmed") {
  if ((estimand == "pmed" && target == 0.5) || (estimand == "ie" && target == 0)) {
    a0 <- as.numeric(stats::logLik(stats::lm(M ~ C, d))) +
      as.numeric(stats::logLik(stats::lm(Y ~ X + M + C, d)))
    b0 <- as.numeric(stats::logLik(stats::lm(M ~ X + C, d))) +
      as.numeric(stats::logLik(stats::lm(Y ~ X + C, d)))
    return(max(a0, b0))
  }
  dm <- stats::model.matrix(~ X + C, d)
  dy <- stats::model.matrix(~ X + M + C, d)
  km <- ncol(dm); ky <- ncol(dy)
  im <- seq_len(km); iy <- km + seq_len(ky)
  ia <- which(colnames(dm) == "X")
  ib <- km + which(colnames(dy) == "M")
  ivm <- km + ky + 1; ivy <- ivm + 1
  fm <- stats::lm(M ~ X + C, d); fy <- stats::lm(Y ~ X + M + C, d)
  initial <- c(stats::coef(fm), stats::coef(fy),
    log(mean(stats::residuals(fm)^2)), log(mean(stats::residuals(fy)^2)))
  objective <- function(p) {
    vm <- exp(p[ivm]); vy <- exp(p[ivy])
    em <- d$M - as.numeric(dm %*% p[im])
    ey <- d$Y - as.numeric(dy %*% p[iy])
    list(objective = -sum(stats::dnorm(d$M, dm %*% p[im], sqrt(vm), log = TRUE)) -
           sum(stats::dnorm(d$Y, dy %*% p[iy], sqrt(vy), log = TRUE)),
         gradient = c(-as.numeric(crossprod(dm, em)) / vm,
                      -as.numeric(crossprod(dy, ey)) / vy,
                      nrow(d) / 2 - sum(em^2) / (2 * vm),
                      nrow(d) / 2 - sum(ey^2) / (2 * vy)))
  }
  constraint <- function(p) {
    a <- p[ia]; b <- p[ib]; vm <- exp(p[ivm]); vy <- exp(p[ivy])
    if (estimand == "ie") {
      gradient <- numeric(length(p))
      gradient[ia] <- b; gradient[ib] <- a
      return(list(constraints = a * b - target, jacobian = gradient))
    }
    v <- 2 * (b^2 * vm + vy)
    q <- a * delta * b / sqrt(v)
    gradient <- numeric(length(p))
    gradient[ia] <- delta * b / sqrt(v)
    gradient[ib] <- a * delta / sqrt(v) - 2 * a * delta * b^2 * vm / v^1.5
    gradient[ivm] <- -a * delta * b^3 * vm / v^1.5
    gradient[ivy] <- -a * delta * b * vy / v^1.5
    list(constraints = stats::pnorm(q) - target,
         jacobian = gradient * stats::dnorm(q))
  }
  fits <- list()
  for (sa in c(-1, 1)) for (scale in scales) {
    start <- initial
    start[ivm] <- start[ivm] + log(scale)
    start[ivy] <- start[ivy] + log(scale)
    if (estimand == "ie") {
      start[ib] <- sa * sign(target) * max(abs(initial[ib]), 0.1)
      start[ia] <- target / start[ib]
    } else {
      start[ib] <- sa * sign(stats::qnorm(target) * delta) * max(abs(initial[ib]), 0.1)
      start[ia] <- sa * abs(stats::qnorm(target)) *
        sqrt(2 * (start[ib]^2 * exp(start[ivm]) + exp(start[ivy]))) / abs(delta * start[ib])
    }
    fit <- nloptr::nloptr(start, eval_f = objective, eval_g_eq = constraint,
      lb = c(rep(-Inf, km + ky), -20, -20), ub = c(rep(Inf, km + ky), 20, 20),
      opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-10,
                  ftol_rel = 1e-11, maxeval = 4000))
    if (fit$status %in% 1:4 &&
        abs(constraint(fit$solution)$constraints) <= 1e-8) {
      fits[[length(fits) + 1L]] <- -objective(fit$solution)$objective
    }
  }
  if (!length(fits)) stop("Independent oracle did not converge")
  max(unlist(fits))
}
