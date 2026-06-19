#' Incremental Mediated Elasticity Result
#'
#' @description
#' S7 class for the **incremental mediated elasticity** `P_med^delta(delta)`
#' (Paper 2 of the P_med-for-modern-estimands program). Rather than a single
#' number, the estimand is a *curve*: the share of the marginal (derivative-scale)
#' total effect that flows through the mediator when the treatment-assignment odds
#' are tilted by a factor `delta`.
#'
#' Writing `Theta(delta_1, delta_2)` for the mean outcome under a `delta_1`-tilt of
#' the direct pathway and a `delta_2`-tilt of the mediated pathway, the elasticities
#' are `theta'_dir = d/d delta_1 Theta`, `theta'_med = d/d delta_2 Theta`, and
#' `theta'_tot = d/d delta Theta(delta, delta)`. The multivariate chain rule gives
#' the **exact decomposition** `theta'_tot = theta'_dir + theta'_med` (remainder
#' `R(delta) = 0` for every `delta`), and `P_med^delta(delta) = theta'_med / theta'_tot`.
#'
#' Under a linear-Gaussian model with no `A * M` interaction the curve is flat and
#' equals the classical proportion mediated; it bends with `delta` only when an
#' interaction is present.
#'
#' @param curve Data frame: one row per `delta`, columns `delta`, `dir`, `med`,
#'   `tot`, `Pmed`, `se`, `lo`, `hi`.
#' @param deltas Numeric: the tilt factors at which the curve was evaluated.
#' @param method Character: estimation method.
#' @param n Integer: sample size.
#' @param ci_level Numeric: confidence level.
#' @param call Call: original call.
#'
#' @export
IncrPmedResult <- S7::new_class(
  "IncrPmedResult", package = "probmed",
  properties = list(
    curve = S7::class_data.frame,
    deltas = S7::class_numeric,
    method = S7::class_character,
    n = S7::class_integer,
    ci_level = S7::class_numeric,
    call = S7::new_property(class = S7::class_any, default = NULL)
  ),
  validator = function(self) {
    if (nrow(self@curve) && !all(c("delta", "dir", "med", "tot", "Pmed") %in% names(self@curve)))
      "curve must contain columns delta, dir, med, tot, Pmed"
  }
)

# internal: cross-fit one-step corner pseudo-outcomes phi[i, (a, a')].
# Identical machinery to .gauge_fit(): the conditional mean of phi[, "aa'"] given C
# is the nested corner mean nu(a, a', C). Returns the influence matrix and the
# fitted propensity g(C) = P(A = 1 | C) used to build the tilt weights.
.incr_fit <- function(d, K, binY, covars) {
  cf <- paste(covars, collapse = " + ")
  f_pi <- stats::as.formula(paste("A ~", cf))
  f_q  <- stats::as.formula(paste("A ~ M +", cf))
  f_om <- stats::as.formula(paste("Y ~ A * M +", cf))
  n <- nrow(d); folds <- sample(rep(1:K, length.out = n))
  nm <- c("11", "10", "01", "00"); cor <- list(c(1, 1), c(1, 0), c(0, 1), c(0, 0))
  phi <- matrix(0, n, 4, dimnames = list(NULL, nm))
  gC <- numeric(n)
  for (k in 1:K) {
    tr <- d[folds != k, , drop = FALSE]
    te_i <- which(folds == k); te <- d[te_i, , drop = FALSE]
    pim <- stats::glm(f_pi, data = tr, family = stats::binomial())
    qm  <- stats::glm(f_q,  data = tr, family = stats::binomial())
    om  <- stats::glm(f_om, data = tr,
                      family = if (binY) stats::binomial() else stats::gaussian())
    p1 <- stats::predict(pim, newdata = te, type = "response"); gC[te_i] <- p1
    pa <- function(z) ifelse(z == 1, p1, 1 - p1)
    q1 <- stats::predict(qm, newdata = te, type = "response")
    qa <- function(z) ifelse(z == 1, q1, 1 - q1)
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
      phi[te_i, j] <- (te$A == a) / pa(a) * ratio * (te$Y - muAM) +
                      (te$A == ap) / pa(ap) * (muAM - eta) + eta
    }
  }
  list(phi = phi, g = gC)
}

#' Incremental Mediated Elasticity Curve
#'
#' @description
#' Estimate the incremental mediated elasticity `P_med^delta(delta)` -- the
#' derivative-scale proportion mediated as a function of the treatment-tilt factor
#' `delta`. The estimator is cross-fitted and built on the same triply-robust
#' corner pseudo-outcomes as [ward_residual()]; standard errors use the
#' ratio-identity influence function for `med / tot` at each `delta`.
#'
#' @details
#' At each `delta` the tilted assignment probability is
#' `q = delta * g / (delta * g + 1 - g)` with derivative
#' `q' = g (1 - g) / (delta * g + 1 - g)^2`, where `g = P(A = 1 | C)`. The
#' direct and mediated elasticities are weighted contrasts of the four corner
#' pseudo-outcomes by `q` and `q'`; by Theorem 1 (multivariate chain rule) they
#' sum to the total elasticity exactly, so `P_med^delta = med / (dir + med)`.
#'
#' **Standard-error caveat.** The reported SE uses the ratio-identity influence
#' function for the corner-mean pseudo-outcomes but treats the estimated
#' propensity `g(C)` entering the tilt weights `q`, `q'` as fixed. The interval is
#' valid when `g` is correctly specified (as here, a logistic working model that
#' is also the data-generating mechanism in the package tests) and is mildly
#' conservative otherwise. The full efficient influence function carries an
#' additional g-score term; adding it is tracked for a future release.
#'
#' @param object A data frame with columns `A` (binary treatment), `M`
#'   (mediator), `Y` (outcome), and the covariates named in `covars`.
#' @param deltas Numeric: tilt factors at which to evaluate the curve.
#' @param covars Character: covariate column names (default `"C"`).
#' @param K Integer: number of cross-fitting folds (default `5`).
#' @param ci_level Numeric: confidence level (default `0.95`).
#' @param seed Integer: random seed for fold assignment (default `1`).
#' @param ... Unused.
#'
#' @return An [IncrPmedResult] object whose `@curve` slot holds one row per
#'   `delta` with the direct, mediated, and total elasticities, the proportion
#'   mediated `Pmed`, and its standard error and Wald interval.
#'
#' @examples
#' set.seed(1)
#' n <- 1500
#' C <- rnorm(n)
#' A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
#' M <- 0.6 * A + 0.4 * C + rnorm(n)
#' Y <- 0.5 * A + 0.7 * M + 0.3 * C + rnorm(n)
#' d <- data.frame(A, M, Y, C)
#' incr_pmed(d, deltas = c(0.5, 1, 2))
#'
#' @export
incr_pmed <- S7::new_generic("incr_pmed", dispatch_args = "object")

#' @export
S7::method(incr_pmed, S7::class_data.frame) <-
  function(object, deltas = c(1 / 3, 1 / 2, 1, 2, 3), covars = "C",
           K = 5L, ci_level = 0.95, seed = 1L, ...) {
    stopifnot(all(c("A", "M", "Y") %in% names(object)), all(covars %in% names(object)))
    set.seed(seed)
    n <- nrow(object)
    binY <- all(object$Y %in% 0:1)
    fit <- .incr_fit(object, K, binY, covars)
    phi <- fit$phi; g <- fit$g
    zc <- stats::qnorm(1 - (1 - ci_level) / 2)
    rows <- lapply(deltas, function(del) {
      q  <- del * g / (del * g + 1 - g)
      qp <- g * (1 - g) / (del * g + 1 - g)^2
      a_dir <- qp * (q * (phi[, "11"] - phi[, "01"]) + (1 - q) * (phi[, "10"] - phi[, "00"]))
      a_med <- qp * (q * (phi[, "11"] - phi[, "10"]) + (1 - q) * (phi[, "01"] - phi[, "00"]))
      dir <- mean(a_dir); med <- mean(a_med); tot <- dir + med; Pmed <- med / tot
      a_tot <- a_dir + a_med
      psi <- ((a_med - med) - Pmed * (a_tot - tot)) / tot   # ratio-identity IF
      se <- stats::sd(psi) / sqrt(n)
      data.frame(delta = del, dir = dir, med = med, tot = tot, Pmed = Pmed,
                 se = se, lo = Pmed - zc * se, hi = Pmed + zc * se)
    })
    curve <- do.call(rbind, rows)
    rownames(curve) <- NULL
    IncrPmedResult(
      curve = curve, deltas = as.numeric(deltas), method = "onestep-crossfit",
      n = as.integer(n), ci_level = ci_level, call = match.call()
    )
  }

#' @export
S7::method(print, IncrPmedResult) <- function(x, ...) {
  cat("Incremental mediated elasticity P_med^delta(delta) (",
      x@method, ", n=", x@n, ")\n", sep = "")
  cn <- x@curve
  for (i in seq_len(nrow(cn))) {
    cat(sprintf("  delta=%-5.3g  P_med=%.3f  [%.3f, %.3f]  (dir=%.3f med=%.3f tot=%.3f)\n",
                cn$delta[i], cn$Pmed[i], cn$lo[i], cn$hi[i],
                cn$dir[i], cn$med[i], cn$tot[i]))
  }
  if (diff(range(cn$Pmed)) > 0.05)
    cat("  curve bends with delta => treatment-by-mediator interaction present.\n")
  invisible(x)
}
