#' Wasserstein / Transport-Scale Proportion Mediated Result
#'
#' @description
#' S7 class for the **Wasserstein proportion mediated**
#' `P_med^W = NIE^W / (NIE^W + NDE^W)`, where `NIE^W = W_2(nu_11, nu_10)` and
#' `NDE^W = W_2(nu_10, nu_00)` are the 2-Wasserstein distances between the
#' corner laws `nu(a, a') = Law(Y(a, M(a')))`. Unlike the additive split
#' `TE = NDE + NIE`, the Wasserstein decomposition handles treatment-by-mediator
#' interaction via the synergy term `G^W = NIE^W + NDE^W - TE^W >= 0`
#' (triangle gap).
#'
#' @param pmedW Numeric: transport-scale proportion mediated.
#' @param pmedW_ci Numeric length-2: confidence interval (Wald or bootstrap).
#' @param pmedW_se Numeric: standard error.
#' @param NIE_W,NDE_W,TE_W Numeric: Wasserstein NIE, NDE, and total effect.
#' @param synergy Numeric: triangle gap `NIE^W + NDE^W - TE^W >= 0`.
#' @param boundary Logical: `TRUE` when the estimator is near the NIE^W = 0 boundary.
#' @param method Character: estimation method (`"eif"`, `"bootstrap"`, or `"dr"`).
#' @param n Integer: sample size.
#' @param ci_level Numeric: confidence level.
#' @param call Call: original call.
#'
#' @export
WassersteinPmedResult <- S7::new_class(
  "WassersteinPmedResult", package = "probmed",
  properties = list(
    pmedW    = S7::class_numeric,
    pmedW_ci = S7::class_numeric,
    pmedW_se = S7::class_numeric,
    NIE_W    = S7::class_numeric,
    NDE_W    = S7::class_numeric,
    TE_W     = S7::class_numeric,
    synergy  = S7::class_numeric,
    boundary = S7::new_property(class = S7::class_logical, default = FALSE),
    method   = S7::class_character,
    n        = S7::class_integer,
    ci_level = S7::class_numeric,
    call     = S7::new_property(class = S7::class_any, default = NULL)
  ),
  validator = function(self) {
    if (length(self@NIE_W) && self@NIE_W < 1e-4)
      warning("NIE^W near 0: P_med^W is near the non-regular boundary; ",
              "SE via EIF degenerates. Bootstrap SE is reported instead.")
    NULL
  }
)

## ===================== transport + EIF engine (VERIFIED, d=1) =====================

# 1-D quadratic Wasserstein via quantile transport (Brenier monotone rearrangement).
.w2_1d <- function(y1, y2, ngrid = 4096L) {
  u <- (seq_len(ngrid) - 0.5) / ngrid
  sqrt(mean((stats::quantile(y1, u, names = FALSE, type = 7) -
             stats::quantile(y2, u, names = FALSE, type = 7))^2))
}

# Kantorovich-potential influence function of W_2^2(nu_a, nu_b) in its FIRST
# argument evaluated at corner-a draws `ya` (transport toward corner-b draws `yb`).
# phi(x) = x^2 - 2 int_0^x T(t) dt, T = Q_b o F_a; IF = phi(x) - E_a[phi].
# Verified FACT5 in pmed-modern/04-wasserstein-pmed/sims/wasserstein_verify.R.
.potIF <- function(ya, yb, ngrid = 4096L) {
  u  <- (seq_len(ngrid) - 0.5) / ngrid
  Qa <- stats::quantile(ya, u, names = FALSE, type = 7)
  Qb <- stats::quantile(yb, u, names = FALSE, type = 7)
  Tmap <- stats::approxfun(Qa, Qb, rule = 2)
  grid <- seq(min(ya), max(ya), length.out = 4096L)
  Tg   <- Tmap(grid); dt <- diff(grid)
  cumI <- c(0, cumsum((utils::head(Tg, -1) + utils::tail(Tg, -1)) / 2 * dt))
  Ifun <- stats::approxfun(grid, cumI, rule = 2)
  p <- ya^2 - 2 * (Ifun(ya) - Ifun(0))
  p - mean(p)
}

# Point estimate + EIF-based SE from corner pseudo-samples (3 independent draws, d=1).
# EIF of P_med^W: delta method through phi_{W2}=potIF/(2 W2) and phi_P=((1-P)phi_NIE - P phi_NDE)/D.
.pmedW_engine <- function(nu11, nu10, nu00, eif = TRUE) {
  NIE <- .w2_1d(nu11, nu10); NDE <- .w2_1d(nu10, nu00); TEW <- .w2_1d(nu11, nu00)
  D   <- NIE + NDE; P <- NIE / D
  out <- list(NIE_W = NIE, NDE_W = NDE, TE_W = TEW, synergy = NIE + NDE - TEW,
              pmedW = P, se = NA_real_)
  if (eif && NIE > 1e-6 && NDE > 1e-6) {
    n   <- min(length(nu11), length(nu10), length(nu00))
    c11 <- ((1 - P) * .potIF(nu11, nu10) / (2 * NIE)) / D
    c10 <- ((1 - P) * .potIF(nu10, nu11) / (2 * NIE) -
             P       * .potIF(nu10, nu00) / (2 * NDE)) / D
    c00 <- (-P * .potIF(nu00, nu10) / (2 * NDE)) / D
    out$se <- sqrt((stats::var(c11) + stats::var(c10) + stats::var(c00)) / n)
  }
  out
}

## ===================== corner-law estimator: MC g-computation ====================

# Parametric g-computation of the three corner laws via Monte-Carlo.
# Draws Mstar ~ N(muM(a'), sM) and Ystar ~ N(muY(a, Mstar, C), sY) or Bernoulli.
.corner_laws_gcomp <- function(d, covars, n_mc = 20000L, binY = FALSE) {
  cf  <- paste(covars, collapse = " + ")
  m_M <- stats::lm(stats::as.formula(paste("M ~ A +", cf)), data = d)
  f_Y <- stats::as.formula(paste("Y ~ A * M +", cf))
  m_Y <- stats::glm(f_Y, data = d,
                    family = if (binY) stats::binomial() else stats::gaussian())
  sM  <- stats::sigma(m_M)
  sY  <- if (binY) NA_real_ else stats::sigma(m_Y)
  Cdraw <- d[sample.int(nrow(d), n_mc, replace = TRUE), covars, drop = FALSE]
  corner <- function(a, ap) {
    nd_M  <- data.frame(A = ap, Cdraw)
    muM   <- stats::predict(m_M, newdata = nd_M)
    Mstar <- muM + stats::rnorm(n_mc, 0, sM)
    nd_Y  <- data.frame(A = a, M = Mstar, Cdraw)
    muY   <- stats::predict(m_Y, newdata = nd_Y, type = "response")
    if (binY) stats::rbinom(n_mc, 1, muY) else muY + stats::rnorm(n_mc, 0, sY)
  }
  list(nu11 = corner(1, 1), nu10 = corner(1, 0), nu00 = corner(0, 0))
}

## ====================== cross-fit DR corner-CDF estimator =======================

# Quantile grid for the DR corner-CDF inversion.
.U_dr <- (seq_len(2048L) - 0.5) / 2048L

# Cross-fit triply-robust corner CDF F_{a,a'}(y) estimator.
# Condition: EITHER outcome OR propensity model may be misspecified (mediator correct).
# Returns the quantile function of the estimated corner law on .U_dr.
.corner_cdf_dr <- function(d, a, ap, covars, ygrid, K = 5L, misspec = "none") {
  n <- nrow(d); folds <- sample(rep(1:K, length.out = n))
  cf <- paste(covars, collapse = " + ")
  G <- length(ygrid); Fmat <- matrix(NA_real_, n, G)
  for (k in 1:K) {
    tr  <- d[folds != k, ]; idx <- which(folds == k); te <- d[idx, ]
    f_pi <- if (misspec == "propensity") stats::as.formula("A ~ 1") else
              stats::as.formula(paste("A ~", cf))
    f_oY <- if (misspec == "outcome")   stats::as.formula(paste("Y ~ A +", cf)) else
              stats::as.formula(paste("Y ~ A * M +", cf))
    pim  <- stats::glm(f_pi, data = tr, family = stats::binomial())
    mMm  <- stats::lm(stats::as.formula(paste("M ~ A +", cf)), data = tr)
    oYm  <- stats::lm(f_oY, data = tr)
    sM   <- stats::sigma(mMm); sY <- stats::sigma(oYm)
    p1   <- stats::predict(pim, te, type = "response")
    pa   <- function(z) ifelse(z == 1, p1, 1 - p1)
    mMhat  <- function(z) stats::predict(mMm, transform(te, A = z))
    muY    <- function(z, mv) stats::predict(oYm,
                data.frame(A = z, M = mv, te[, covars, drop = FALSE]))
    r      <- stats::dnorm(te$M, mMhat(ap), sM) / stats::dnorm(te$M, mMhat(a), sM)
    alpha  <- muY(a, rep(0, nrow(te))); beta_eff <- muY(a, rep(1, nrow(te))) - alpha
    muY_obs <- alpha + beta_eff * te$M
    mM_ap   <- mMhat(ap); sd_eta <- sqrt(sY^2 + beta_eff^2 * sM^2)
    ind     <- (te$A == a)  / pa(a)  * r
    wap     <- (te$A == ap) / pa(ap)
    for (g in seq_len(G)) {
      y <- ygrid[g]
      Fmat[idx, g] <- ind  * ((te$Y <= y) - stats::pnorm(y, muY_obs, sY)) +
                      wap  * (stats::pnorm(y, muY_obs, sY) -
                               stats::pnorm(y, alpha + beta_eff * mM_ap, sd_eta)) +
                      stats::pnorm(y, alpha + beta_eff * mM_ap, sd_eta)
    }
  }
  Fhat <- cummax(pmin(pmax(colMeans(Fmat), 0), 1))
  stats::approx(Fhat, ygrid, xout = .U_dr, rule = 2, ties = "ordered")$y
}

#' Doubly-Robust Wasserstein Proportion Mediated
#'
#' @description
#' Standalone cross-fit DR estimator of `P_med^W`. Consistent when EITHER the
#' outcome model OR the propensity model is misspecified, provided the mediator
#' model is correct (multiply robust; verified in simulation).
#'
#' @param d A `data.frame` with columns `A`, `M`, `Y`, and the covariates in `covars`.
#' @param covars Character vector of covariate column names. Default `"C"`.
#' @param K Integer: cross-fitting folds (default 5).
#' @param misspec Character: `"none"` (both models correct), `"outcome"`, or
#'   `"propensity"` — for robustness testing only.
#' @param G Integer: number of `y`-grid points for the CDF estimator (default 300).
#' @param seed Integer: RNG seed for fold assignment (default 1).
#'
#' @return A named list with `NIE_W`, `NDE_W`, and `pmedW`.
#'
#' @export
pmedW_dr <- function(d, covars = "C", K = 5L, misspec = "none", G = 300L, seed = 1L) {
  set.seed(seed)
  yg  <- seq(min(d$Y) - 1, max(d$Y) + 1, length.out = G)
  q11 <- .corner_cdf_dr(d, 1, 1, covars, yg, K, misspec)
  q10 <- .corner_cdf_dr(d, 1, 0, covars, yg, K, misspec)
  q00 <- .corner_cdf_dr(d, 0, 0, covars, yg, K, misspec)
  NIE <- sqrt(mean((q11 - q10)^2))
  NDE <- sqrt(mean((q10 - q00)^2))
  list(NIE_W = NIE, NDE_W = NDE, pmedW = NIE / (NIE + NDE))
}

## ===================== entropic OT (Sinkhorn) for d > 1 =========================

# Log-domain Sinkhorn iteration (Genevay et al.). Returns transport cost.
.sinkhorn_cost <- function(X, Y, eps, iters = 400L, tol = 1e-9) {
  X <- as.matrix(X); Y <- as.matrix(Y); n <- nrow(X); m <- nrow(Y)
  C <- outer(rowSums(X^2), rowSums(Y^2), "+") - 2 * X %*% t(Y); C[C < 0] <- 0
  la <- log(rep(1 / n, n)); lb <- log(rep(1 / m, m)); f <- numeric(n); g <- numeric(m)
  lse <- function(v) { mx <- max(v); mx + log(sum(exp(v - mx))) }
  for (it in seq_len(iters)) {
    fo <- f
    for (i in seq_len(n)) f[i] <- -eps * lse(lb + (g - C[i, ]) / eps)
    for (j in seq_len(m)) g[j] <- -eps * lse(la + (f - C[, j]) / eps)
    if (max(abs(f - fo)) < tol) break
  }
  sum(exp((outer(f, g, "+") - C) / eps + outer(la, lb, "+")) * C)
}

#' Debiased Sinkhorn Divergence
#'
#' @description
#' Debiased Sinkhorn divergence `S_eps(X, Y)` (Genevay et al.), which approximates
#' `W_2^2(X, Y)` as `eps -> 0`. Use for multivariate (`d > 1`) Wasserstein P_med^W
#' via [pmedW_md()].
#'
#' @param X,Y Numeric matrices (rows = samples, columns = dimensions).
#' @param eps Numeric: regularization parameter (default 0.1).
#' @param iters Integer: maximum Sinkhorn iterations (default 400).
#'
#' @return Non-negative scalar.
#'
#' @export
sinkhorn_div <- function(X, Y, eps = 0.1, iters = 400L) {
  max(.sinkhorn_cost(X, Y, eps, iters) -
      0.5 * .sinkhorn_cost(X, X, eps, iters) -
      0.5 * .sinkhorn_cost(Y, Y, eps, iters), 0)
}

#' Multivariate Wasserstein Proportion Mediated
#'
#' @description
#' Computes `P_med^W` for multivariate outcomes (`d >= 1`) from corner pseudo-samples
#' via debiased Sinkhorn divergence. Intended for use with pre-computed corner-law draws
#' (e.g. from g-computation or DR). SE not yet implemented.
#'
#' @param NU11,NU10,NU00 Numeric matrices (rows = draws, columns = outcome dimensions).
#' @param eps Numeric: Sinkhorn regularization parameter (default 0.1).
#'
#' @return A named list with `NIE_W`, `NDE_W`, and `pmedW`.
#'
#' @export
pmedW_md <- function(NU11, NU10, NU00, eps = 0.1) {
  NIE <- sqrt(sinkhorn_div(NU11, NU10, eps))
  NDE <- sqrt(sinkhorn_div(NU10, NU00, eps))
  list(NIE_W = NIE, NDE_W = NDE, pmedW = NIE / (NIE + NDE))
}

## ========================= S7 generic + method ==================================

#' Wasserstein / Transport-Scale Proportion Mediated
#'
#' @description
#' Estimate `P_med^W = NIE^W / (NIE^W + NDE^W)`, the transport-scale proportion
#' mediated, where `NIE^W = W_2(nu_11, nu_10)` and `NDE^W = W_2(nu_10, nu_00)` are
#' the 2-Wasserstein distances between the corner laws `nu(a, a')`.
#'
#' Three estimation methods are available:
#' \itemize{
#'   \item `"eif"` (default): parametric g-computation corner samples + EIF-based SE.
#'     Valid when the g-comp working models are correct. Fast; d = 1 only.
#'   \item `"bootstrap"`: parametric g-computation + bootstrap CI. Slower but does not
#'     require the EIF to be non-degenerate; used automatically near the NIE^W = 0
#'     boundary.
#'   \item `"dr"`: cross-fit DR corner-CDF estimator (multiply robust: consistent when
#'     either the outcome or propensity model is wrong, mediator model correct) with
#'     bootstrap SE. Slowest; d = 1 only.
#' }
#'
#' @param object A `data.frame` with columns `A` (binary 0/1), `M` (mediator),
#'   `Y` (outcome), and the covariates named in `covars`.
#' @param covars Character vector of covariate column names. Default `"C"`.
#' @param n_mc Integer: Monte-Carlo draws for g-comp corner laws (default 20000).
#'   Ignored when `method = "dr"`.
#' @param binY Logical: if `TRUE`, treats `Y` as binary (uses logistic g-comp).
#' @param method Character: `"eif"` (default), `"bootstrap"`, or `"dr"`.
#' @param n_boot Integer: bootstrap replications (default 200). Used by
#'   `"bootstrap"` and `"dr"` methods, and automatically when near the boundary.
#' @param K Integer: cross-fitting folds for `method = "dr"` (default 5).
#' @param G Integer: y-grid points for `method = "dr"` (default 300).
#' @param ci_level Numeric: confidence level (default 0.95).
#' @param seed Integer: RNG seed (default 1).
#' @param ... Unused.
#'
#' @return A [WassersteinPmedResult] object.
#'
#' @examples
#' set.seed(1)
#' n <- 800; C <- rnorm(n)
#' A <- rbinom(n, 1, plogis(0.3 * C))
#' M <- 0.7 * A + 0.5 * C + rnorm(n)
#' Y <- 0.6 * A + 0.8 * M + 0.3 * C + rnorm(n)
#' wasserstein_pmed(data.frame(A, M, Y, C), covars = "C")
#'
#' @export
wasserstein_pmed <- S7::new_generic(
  "wasserstein_pmed", dispatch_args = "object",
  fun = function(object, covars = "C", n_mc = 20000L, binY = FALSE,
                 method = c("eif", "bootstrap", "dr"),
                 n_boot = 200L, K = 5L, G = 300L,
                 ci_level = 0.95, seed = 1L, ...) {
    S7::S7_dispatch()
  })

#' @export
S7::method(wasserstein_pmed, S7::class_data.frame) <-
  function(object, covars = "C", n_mc = 20000L, binY = FALSE,
           method = c("eif", "bootstrap", "dr"),
           n_boot = 200L, K = 5L, G = 300L,
           ci_level = 0.95, seed = 1L, ...) {
    stopifnot(all(c("A", "M", "Y") %in% names(object)), all(covars %in% names(object)))
    method <- match.arg(method); set.seed(seed)

    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    on_boundary <- FALSE

    if (method == "dr") {
      ## DR path: cross-fit corner CDF -> W2 quantile distances
      yg  <- seq(min(object$Y) - 1, max(object$Y) + 1, length.out = G)
      q11 <- .corner_cdf_dr(object, 1, 1, covars, yg, K, "none")
      q10 <- .corner_cdf_dr(object, 1, 0, covars, yg, K, "none")
      q00 <- .corner_cdf_dr(object, 0, 0, covars, yg, K, "none")
      NIE <- sqrt(mean((q11 - q10)^2))
      NDE <- sqrt(mean((q10 - q00)^2))
      TEW <- sqrt(mean((q11 - q00)^2))
      D   <- NIE + NDE; P <- NIE / D
      on_boundary <- NIE < 0.05 * D
      ## bootstrap SE on the DR path
      bs <- vapply(seq_len(n_boot), function(b) {
        db  <- object[sample.int(nrow(object), replace = TRUE), , drop = FALSE]
        pmedW_dr(db, covars = covars, K = K, G = G, seed = b)$pmedW
      }, numeric(1))
      se <- stats::sd(bs, na.rm = TRUE)
      ci <- stats::quantile(bs, c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2),
                            names = FALSE, na.rm = TRUE)
      synergy <- NIE + NDE - TEW
      method_out <- "dr"
    } else {
      ## g-comp path: parametric Monte-Carlo corner laws
      cl  <- .corner_laws_gcomp(object, covars, n_mc, binY)
      eng <- .pmedW_engine(cl$nu11, cl$nu10, cl$nu00, eif = (method == "eif"))
      NIE <- eng$NIE_W; NDE <- eng$NDE_W; TEW <- eng$TE_W; D <- NIE + NDE; P <- eng$pmedW
      synergy <- eng$synergy
      on_boundary <- NIE < 0.05 * D

      if (method == "eif" && !on_boundary) {
        se <- eng$se
        ci <- pmin(pmax(c(P - z * se, P + z * se), 0), 1)
        method_out <- "eif"
      } else {
        ## boundary or bootstrap requested: fall back to bootstrap CI
        if (method == "eif" && on_boundary)
          message("wasserstein_pmed: NIE^W near 0 (boundary); switching to bootstrap CI.")
        bs <- vapply(seq_len(n_boot), function(b) {
          db  <- object[sample.int(nrow(object), replace = TRUE), , drop = FALSE]
          clb <- .corner_laws_gcomp(db, covars, n_mc, binY)
          .pmedW_engine(clb$nu11, clb$nu10, clb$nu00, eif = FALSE)$pmedW
        }, numeric(1))
        se <- stats::sd(bs, na.rm = TRUE)
        ci <- stats::quantile(bs, c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2),
                              names = FALSE, na.rm = TRUE)
        method_out <- if (on_boundary) "bootstrap-boundary" else "bootstrap"
      }
    }

    WassersteinPmedResult(
      pmedW    = P,
      pmedW_ci = pmin(pmax(ci, 0), 1),
      pmedW_se = se,
      NIE_W    = NIE, NDE_W = NDE, TE_W = TEW, synergy = synergy,
      boundary = on_boundary,
      method   = method_out,
      n        = as.integer(nrow(object)),
      ci_level = ci_level,
      call     = match.call()
    )
  }

#' @export
S7::method(print, WassersteinPmedResult) <- function(x, ...) {
  cat(sprintf("Wasserstein P_med^W (%s, n=%d)\n", x@method, x@n))
  cat(sprintf("  P_med^W = %.3f  [%.3f, %.3f]  (se=%.3f)\n",
              x@pmedW, x@pmedW_ci[1], x@pmedW_ci[2], x@pmedW_se))
  cat(sprintf("  NIE^W=%.3f  NDE^W=%.3f  TE^W=%.3f  G^W(synergy)=%.3f\n",
              x@NIE_W, x@NDE_W, x@TE_W, x@synergy))
  if (x@boundary)
    cat("  ! Near NIE^W=0 boundary: bootstrap CI reported; inference is conservative.\n")
  if (x@synergy > 0.1 * (x@NIE_W + x@NDE_W))
    cat("  ! G^W > 10% of (NIE^W + NDE^W): triangle gap is non-negligible.\n")
  invisible(x)
}
