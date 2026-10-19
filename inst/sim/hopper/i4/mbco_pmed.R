# mbco_pmed.R -- MBCO confidence interval for P_med (continuous single-mediator
# model), implemented via the analytic constraint so the constrained fit is an
# ordinary (well-converging) optimization rather than a nonlinear OpenMx mxConstraint.
#
# Model:  M = iM + a*X + eM,        Var(eM) = Vm
#         Y = iY + cp*X + b*M + eY, Var(eY) = Vy
# P_med = Phi(d),  d = a*b / sqrt(2*(b^2*Vm + Vy)).
#
# Constraint P_med = p*  <=>  d = q* = qnorm(p*).  Solving for Vy:
#         Vy = a^2*b^2 / (2*q*^2) - b^2*Vm.
# We fix Vy by this formula (q* != 0) and maximize the joint likelihood over the
# remaining parameters; this is the MBCO profile log-likelihood. The interval is
# the set of p* the LR test does not reject (Wilks: -2*(L0 - L1) <= chi^2_{1,1-a}).
# This is the same procedure as RMediation::mbco (Tofighi & Kelley, 2020).
#
# Edge cases handled: q*=0 -> null indirect (b=0); Vy<=0 -> infeasible (rejected);
# sign(ab) pinned to sign(q*) so the fit targets p*, not 1-p*.

.ll_joint <- function(iM, a, Vm, iY, cp, b, Vy, data) {
  sum(dnorm(data$M, iM + a * data$X, sqrt(Vm), log = TRUE)) +
  sum(dnorm(data$Y, iY + cp * data$X + b * data$M, sqrt(Vy), log = TRUE))
}

.ll_free <- function(data) {
  as.numeric(logLik(lm(M ~ X, data))) + as.numeric(logLik(lm(Y ~ X + M, data)))
}

.ll_constrained <- function(data, qstar, a_sign) {
  if (abs(qstar) < 1e-8) {
    # p* = 0.5: ab = 0 is a COMPOSITE null, {a = 0} U {b = 0}; fit both branches
    # and keep the better one. Imposing only b = 0 misfits data generated under
    # a = 0, b != 0 and rejects the null every time (0% coverage).
    ll_b0 <- as.numeric(logLik(lm(M ~ X, data))) + as.numeric(logLik(lm(Y ~ X, data)))
    ll_a0 <- as.numeric(logLik(lm(M ~ 1, data))) + as.numeric(logLik(lm(Y ~ X + M, data)))
    return(max(ll_b0, ll_a0))
  }
  fm <- lm(M ~ X, data); fy <- lm(Y ~ X + M, data)
  Vm_f <- sigma(fm)^2; Vy_f <- sigma(fy)^2; a_f <- coef(fm)["X"]; b_f <- coef(fy)["M"]

  # The constraint P_med = p* (i.e. a*b / sqrt(2(b^2 Vm + Vy)) = q*) can be solved
  # for Vy OR for a; the two parametrizations are well-conditioned in complementary
  # regimes, so we optimize both and keep the better fit. This makes the LR test
  # correct across the whole neighborhood of p* = 0.5, where ab = 0 is composite.

  # Fit helper shared by both branches: accept an optim() result only if it
  # actually converged (finding #15 -- unchecked convergence let a premature
  # Nelder-Mead stop masquerade as an optimum and silently corrupt the LR
  # statistic). A non-converged fit gets one retry from its own endpoint
  # (cheap, often clears a stall that wasn't a real optimum); if that also
  # fails to converge, the start is discarded rather than trusted.
  fit_one <- function(start, nll) {
    o <- try(optim(start, nll, method = "Nelder-Mead",
                   control = list(maxit = 3000, reltol = 1e-10)), silent = TRUE)
    if (inherits(o, "try-error")) return(NULL)
    if (o$convergence != 0) {
      o2 <- try(optim(o$par, nll, method = "Nelder-Mead",
                      control = list(maxit = 3000, reltol = 1e-10)), silent = TRUE)
      if (inherits(o2, "try-error") || o2$convergence != 0) return(NULL)
      o <- o2
    }
    o
  }

  # --- Branch 1: solve for Vy (well-conditioned when b is small -> the b = 0 null).
  # Near p* = 0.5, b -> 0 / Vy -> inf; the analytic warm start b ~ |q*| sqrt(2 Vy_f)/|a_f|
  # keeps Vy ~ Vy_f so Nelder-Mead does not under-fit and spuriously inflate the LR stat.
  # `a` is sign-pinned to a_sign (log-magnitude parametrized), mirroring how `b` is
  # already pinned via sb -- Vy depends only on a^2, so this excludes no genuine
  # solution, only the wrong-sign-product ones that would target 1-p* instead of p*
  # (finding #17: the header's claimed "sign(ab) pinned to sign(q*)" previously held
  # only via b's sign, leaving a free to flip and silently break the pin).
  sb <- sign(qstar) * a_sign
  nll1 <- function(p) {
    iM <- p[1]; a <- a_sign * exp(p[2]); Vm <- exp(p[3]); iY <- p[4]; cp <- p[5]; b <- sb * exp(p[6])
    Vy <- a^2 * b^2 / (2 * qstar^2) - b^2 * Vm
    if (!is.finite(Vy) || Vy <= 1e-8) return(1e10)
    -.ll_joint(iM, a, Vm, iY, cp, b, Vy, data)
  }
  b0 <- abs(qstar) * sqrt(2 * Vy_f) / max(abs(a_f), 1e-3)
  best <- 1e10
  for (st in list(
        c(coef(fm)[1], log(abs(a_f) + 1e-3), log(Vm_f), coef(fy)[1], coef(fy)[2], log(abs(b_f) + 1e-3)),
        c(coef(fm)[1], log(abs(a_f) + 1e-3), log(Vm_f), coef(fy)[1], coef(fy)[2], log(b0 + 1e-6)))) {
    o <- fit_one(st, nll1)
    if (!is.null(o) && o$value < best) best <- o$value
  }

  # --- Branch 2: solve for a (well-conditioned when b is bounded from 0 -> the a = 0
  # null). a = q* sqrt(2(b^2 Vm + Vy)) / b; b, Vm, Vy free. Without this branch, data
  # generated under a = 0, b != 0 are misfit and the null is over-rejected.
  # Two starts, mirroring branch 1's asymmetry fix (finding #16): the data-fitted b_f,
  # and its sign-flipped counterpart -b_f. A single from-data start can strand
  # Nelder-Mead on the wrong side of b = 0 under small-N/noisy data, since nothing
  # here (unlike branch 1's sb) pins b's sign to hedge against that.
  nll2 <- function(p) {
    iM <- p[1]; Vm <- exp(p[2]); iY <- p[3]; cp <- p[4]; b <- p[5]; Vy <- exp(p[6])
    if (abs(b) < 1e-6) return(1e10)
    a <- qstar * sqrt(2 * (b^2 * Vm + Vy)) / b
    -.ll_joint(iM, a, Vm, iY, cp, b, Vy, data)
  }
  for (b_start in c(as.numeric(b_f), -as.numeric(b_f))) {
    o2 <- fit_one(c(coef(fm)[1], log(Vm_f), coef(fy)[1], coef(fy)[2], b_start, log(Vy_f)), nll2)
    if (!is.null(o2) && o2$value < best) best <- o2$value
  }

  if (best >= 1e9) return(NA_real_)
  -best
}

#' P_med point estimate and MBCO confidence interval (continuous outcome).
#' Endpoints are found by root-finding on the inverted likelihood-ratio test, so
#' the interval is resolution-independent (no grid artifact at large n).
#' @param data data.frame with columns X, M, Y
#' @param level confidence level (default 0.95)
#' @return list(estimate, ci_lower, ci_upper)
pmed_mbco_ci <- function(data, level = 0.95) {
  fm <- lm(M ~ X, data); fy <- lm(Y ~ X + M, data)
  a <- coef(fm)["X"]; b <- coef(fy)["M"]; Vm <- sigma(fm)^2; Vy <- sigma(fy)^2
  dhat <- a * b / sqrt(2 * (b^2 * Vm + Vy))
  est  <- as.numeric(pnorm(dhat))
  a_sign <- sign(as.numeric(a)); if (a_sign == 0) a_sign <- 1
  ll1 <- .ll_free(data); crit <- qchisq(level, df = 1)
  # excess(ps) = LR statistic at P_med = ps, minus the chi-square cutoff
  excess <- function(ps) {
    ps  <- min(max(ps, 1e-5), 1 - 1e-5)
    ll0 <- .ll_constrained(data, qnorm(ps), a_sign)
    if (is.na(ll0)) return(crit + 1e3)             # infeasible -> rejected
    -2 * (ll0 - ll1) - crit
  }
  # Step outward from the estimate until the statistic exceeds the cutoff, then
  # locate the crossing exactly with uniroot. The step grows geometrically
  # (0.01 * 1.15^k) rather than linearly: a linear step (0.01*k, k<=80, max
  # reach 0.8) can silently under-search a genuinely wide interval when the
  # point estimate itself lands near a boundary from small-N sampling noise --
  # the search would exhaust its budget short of both a real bracket and the
  # [1e-4, 1-1e-4] edge, and previously returned that unverified midpoint as
  # if it were the confidence bound (adversarial-audit finding #14). The
  # geometric schedule keeps fine resolution near the estimate (small early
  # steps, so uniroot still gets a tight bracket for ordinary interior cases)
  # while guaranteeing the full [0,1] range is covered well within the
  # iteration budget: cumulative reach after k steps is 0.01*(1.15^k-1)/0.15,
  # which already exceeds 1 by k ~ 20 even starting from the midpoint.
  endpoint <- function(dir) {
    inner <- est
    outer <- est
    step <- 0.01
    max_iter <- 60
    for (k in 1:max_iter) {
      outer <- min(max(outer + dir * step, 1e-4), 1 - 1e-4)
      if (excess(outer) > 0)
        return(tryCatch(uniroot(excess, sort(c(inner, outer)), tol = 1e-4)$root,
                        error = function(e) NA_real_))
      inner <- outer
      if (outer <= 1e-4 || outer >= 1 - 1e-4) return(outer)   # hit the boundary
      step <- step * 1.15
    }
    # Genuine non-convergence: neither a bracket nor the [0,1] edge was found
    # within max_iter geometric steps, despite the schedule covering the full
    # range several times over. Report honestly rather than guess.
    warning("pmed_mbco_ci(): endpoint search (dir = ", dir, ") found neither a ",
            "bracket nor the [0,1] boundary within ", max_iter, " geometric steps; ",
            "returning NA. This indicates a non-monotonic or unusually shaped ",
            "likelihood-ratio profile -- inspect the data before trusting other outputs.")
    NA_real_
  }
  list(estimate = est, ci_lower = endpoint(-1), ci_upper = endpoint(+1))
}
