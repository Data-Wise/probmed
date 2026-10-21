# internal: cross-fit one-step corner pseudo-outcomes phi[i, (a, a')].
#
# Shared engine for the corner means theta(a, a') = E[Y(a, M(a'))] underlying the
# whole P_med-for-modern-estimands program (gauge / incremental / Sobol). The
# conditional mean of phi[, "aa'"] given C is the nested corner mean nu(a, a', C);
# columns are the four corners "11", "10", "01", "00". Each entry is the
# triply-robust efficient influence function of Tchetgen Tchetgen & Shpitser
# (2012), built from a propensity model pi(C) = P(A = 1 | C), a mediator-density
# proxy q(C) = P(A = 1 | M, C), and an outcome model mu(a, M, C).
#
# Returns the influence matrix `phi` AND the fitted propensity `g(C) = P(A = 1 | C)`;
# `g` is used by the incremental estimator to build the tilt weights and is ignored
# by the gauge / Sobol estimators (which use `$phi` only). Centralizing this here
# removes the corner-EIF duplication that previously lived in `.gauge_fit()` and
# `.incr_fit()`.
#
# @param d data frame with columns A, M, Y and the covariates in `covars`.
# @param K integer number of cross-fitting folds.
# @param binY logical; if TRUE the outcome model is binomial (logistic), else
#   gaussian. Predictions use type = "response" so the corner pseudo-outcomes are
#   on the outcome scale in both cases.
# @param covars character vector of covariate column names.
# @return list(phi = <n x 4 matrix>, g = <length-n numeric>).
.corner_fit <- function(d, K, binY, covars) {
  n <- nrow(d); folds <- sample(rep(1:K, length.out = n))
  phi <- matrix(0, n, 4, dimnames = list(NULL, c("11", "10", "01", "00")))
  gC <- numeric(n)
  for (k in 1:K) {
    te_i <- which(folds == k)
    fk <- .corner_phi(d[folds != k, , drop = FALSE], d[te_i, , drop = FALSE], binY, covars)
    phi[te_i, ] <- fk$phi; gC[te_i] <- fk$g
  }
  list(phi = phi, g = gC)
}

# internal: the same corner influence matrix from FULL-SAMPLE nuisance fits
# (train = test = d). Not a substitute for `.corner_fit()` as a point estimator --
# the in-sample fit forfeits cross-fitting's bias protection. It consumes no RNG,
# which makes it the deterministic probe for the corner EIF used by the tests
# (row-order invariance) and by the se simulations in inst/sim/se_candidates*.R.
# A full-sample-nuisance se for `ward_residual(reps > 1)` was built on it and
# dropped (2026-08-22): once the weight bug below was fixed, the cross-fit IF se
# was already calibrated.
.corner_fit_full <- function(d, binY, covars) .corner_phi(d, d, binY, covars)

# internal: fit the three nuisance models on `tr`, evaluate the corner influence
# matrix on `te`. The eta(a, a', C) projections are fit on the A == a' subset of
# `tr` and predicted on `te`.
.corner_phi <- function(tr, te, binY, covars) {
  cf <- paste(covars, collapse = " + ")
  f_pi <- stats::as.formula(paste("A ~", cf))
  f_q  <- stats::as.formula(paste("A ~ M +", cf))
  f_om <- stats::as.formula(paste("Y ~ A * M +", cf))
  nm <- c("11", "10", "01", "00"); cor <- list(c(1, 1), c(1, 0), c(0, 1), c(0, 0))
  phi <- matrix(0, nrow(te), 4, dimnames = list(NULL, nm))
  pim <- stats::glm(f_pi, data = tr, family = stats::binomial())
  qm  <- stats::glm(f_q,  data = tr, family = stats::binomial())
  om  <- stats::glm(f_om, data = tr,
                    family = if (binY) stats::binomial() else stats::gaussian())
  ## `if`, not `ifelse`: ifelse() returns a result shaped like its TEST, and `z`
  ## is a scalar, so ifelse(z == 1, p1, 1 - p1) silently returned p1[1] -- the
  ## first test row's propensity -- for every observation in the fold. That bug
  ## shipped in 10b41aa (PR #8) and was fixed here (2026-08-22); see NEWS.
  p1 <- stats::predict(pim, newdata = te, type = "response")
  pa <- function(z) if (z == 1) p1 else 1 - p1
  q1 <- stats::predict(qm, newdata = te, type = "response")
  qa <- function(z) if (z == 1) q1 else 1 - q1
  mu <- function(z) stats::predict(om, newdata = transform(te, A = z), type = "response")
  for (j in 1:4) {
    a <- cor[[j]][1]; ap <- cor[[j]][2]
    muAM <- mu(a)
    ratio <- (qa(ap) / qa(a)) * (pa(a) / pa(ap))
    muAM_tr <- stats::predict(om, newdata = transform(tr, A = a), type = "response")
    sub <- tr$A == ap
    etam <- stats::lm(stats::reformulate(covars, "yy"),
                      data = cbind(data.frame(yy = muAM_tr[sub]),
                                   tr[sub, covars, drop = FALSE]))
    eta <- stats::predict(etam, newdata = te)
    phi[, j] <- (te$A == a) / pa(a) * ratio * (te$Y - muAM) +
                (te$A == ap) / pa(ap) * (muAM - eta) + eta
  }
  list(phi = phi, g = p1)
}
