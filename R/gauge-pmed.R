#' Gauge-Calibrated Proportion Mediated Result
#'
#' @description
#' S7 class for the gauge-calibrated proportion mediated. Reports the
#' interventional proportion mediated `P_med = IIE / OE` alongside the
#' **gauge residual** `W = R / OE`, where `R = OE - IDE - IIE` is the
#' non-decomposability (treatment-by-mediator interaction) term. A `W`
#' significantly different from zero signals that the additive split of the
#' overall effect fails and the single-number `P_med` is unreliable.
#'
#' @param p_med Numeric: interventional proportion mediated, `IIE / OE`.
#' @param p_med_ci Numeric length-2: confidence interval for `p_med`.
#' @param p_med_fieller Numeric: Fieller confidence set for the ratio `p_med`
#'   (may be unbounded; empty if not requested).
#' @param fieller_type Character: type of Fieller set returned (e.g. bounded,
#'   unbounded, or empty); `NA` if not computed.
#' @param W Numeric: gauge residual `R / OE`.
#' @param W_ci Numeric length-2: confidence interval for `W`.
#' @param W_se Numeric: standard error of `W`.
#' @param W_p Numeric: two-sided p-value for `H0: W = 0`.
#' @param OE,IDE,IIE,R Numeric: overall, interventional direct, interventional
#'   indirect effects and the remainder.
#' @param theta Numeric length-4: corner means `theta(a, a')`.
#' @param method Character: estimation method.
#' @param n Integer: sample size.
#' @param ci_level Numeric: confidence level.
#' @param se_method Character: `"analytic"` (influence-function se) or
#'   `"bootstrap"` (nonparametric, valid but mildly conservative).
#' @param reps Integer: number of repeated cross-fitting fold draws averaged for
#'   the point estimate.
#' @param call Call: original call.
#'
#' @export
GaugePmedResult <- S7::new_class(
  "GaugePmedResult", package = "probmed",
  properties = list(
    p_med = S7::class_numeric, p_med_ci = S7::class_numeric,
    p_med_fieller = S7::new_property(class = S7::class_numeric, default = numeric(0)),
    fieller_type = S7::new_property(class = S7::class_character, default = NA_character_),
    W = S7::class_numeric, W_ci = S7::class_numeric,
    W_se = S7::class_numeric, W_p = S7::class_numeric,
    OE = S7::class_numeric, IDE = S7::class_numeric,
    IIE = S7::class_numeric, R = S7::class_numeric,
    theta = S7::class_numeric, method = S7::class_character,
    n = S7::class_integer, ci_level = S7::class_numeric,
    se_method = S7::new_property(class = S7::class_character, default = "analytic"),
    reps = S7::new_property(class = S7::class_integer, default = 1L),
    call = S7::new_property(class = S7::class_any, default = NULL)
  ),
  validator = function(self) {
    if (length(self@OE) && abs(self@OE) < 1e-8)
      warning("OE near 0: W = R/OE is unstable; report unnormalized R.")
  }
)

#' Gauge-Calibrated Proportion Mediated
#'
#' @description
#' Estimate the interventional proportion mediated together with the gauge
#' residual `W = R / OE` that flags non-decomposability (treatment-by-mediator
#' interaction). Uses a cross-fitted one-step estimator built on the
#' triply-robust efficient influence functions of the four corner means
#' `theta(a, a') = E[Y(a, M(a'))]`.
#'
#' @param object A `data.frame` with columns `A` (binary treatment), `M`
#'   (mediator), `Y` (outcome), and the covariates named in `covars`.
#' @param covars Character vector of covariate column names. Default `"C"`.
#' @param K Integer: number of cross-fitting folds (default 5).
#' @param ci_level Numeric: confidence level (default 0.95).
#' @param seed Integer: RNG seed for fold assignment.
#' @param fieller Logical: also compute the Fieller confidence set for the ratio
#'   `P_med = IIE/OE`, which is unbounded when the total effect `OE` is not
#'   significant (default `TRUE`).
#' @param reps Integer: number of repeated cross-fitting fold draws averaged for
#'   the point estimate (default `1`). `reps > 1` averages the corner influence
#'   matrix over independent fold assignments, removing the fold-split component
#'   of the variance; the analytic se then adds the residual fold Monte-Carlo
#'   variance of the averaged point.
#' @param se_method Character: `"analytic"` (default, influence-function se,
#'   symmetric Wald CI) or `"bootstrap"`. The analytic se for the ratios `W` and
#'   `P_med` is right-skewed with a median below the empirical SD, so the
#'   symmetric Wald CI is mildly **anti-conservative** (sub-nominal coverage
#'   \eqn{\approx 0.85}-\eqn{0.90}). Under `"bootstrap"` the CIs for `W` and
#'   `P_med` are the tail-aware **percentile** intervals of the nonparametric
#'   bootstrap (resample rows, refit) -- the appropriate construction for a
#'   skewed ratio (widening a symmetric se mis-covers it); `W_se`/`p_med` se are
#'   still reported as bootstrap dispersion summaries. `reps > 1` and
#'   `se_method = "bootstrap"` compose. Bootstrap validity follows from the
#'   estimator being a Neyman-orthogonal cross-fit (DML-type) functional
#'   (Lin et al. 2026); it requires the ratio to be regular, i.e. `OE` bounded
#'   away from 0 -- see the Fieller diagnostic for the near-null case.
#' @param B Integer: number of bootstrap resamples when `se_method = "bootstrap"`
#'   (default `200`). Cost is `B` (x `reps`) refits.
#' @param ... Unused.
#'
#' @return A [GaugePmedResult] object.
#'
#' @examples
#' set.seed(1)
#' n <- 800; C <- rnorm(n)
#' A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
#' M <- 0.6 * A + 0.4 * C + rnorm(n)
#' Y <- 0.5 * A + 0.7 * M + 0.8 * A * M + 0.3 * C + rnorm(n)
#' ward_residual(data.frame(A, M, Y, C))
#'
#' @export
ward_residual <- S7::new_generic(
  "ward_residual", dispatch_args = "object",
  fun = function(object, covars = "C", K = 5L, ci_level = 0.95,
                 seed = 1L, fieller = TRUE, reps = 1L,
                 se_method = c("analytic", "bootstrap"), B = 200L, ...) {
    S7::S7_dispatch()
  })

#' @export
S7::method(ward_residual, S7::class_data.frame) <-
  function(object, covars = "C", K = 5L, ci_level = 0.95, seed = 1L, fieller = TRUE,
           reps = 1L, se_method = c("analytic", "bootstrap"), B = 200L, ...) {
    stopifnot(all(c("A", "M", "Y") %in% names(object)), all(covars %in% names(object)))
    se_method <- match.arg(se_method)
    reps <- max(1L, as.integer(reps))
    set.seed(seed)
    binY <- all(object$Y %in% 0:1); n <- nrow(object)

    ## ---- repeated cross-fitting (reps): average the corner influence matrix over
    ## `reps` independent fold draws, removing the fold-split variance component.
    ## reps = 1 reproduces the single-draw estimator exactly.
    .gp_WP <- function(p) {                       # (W, P_med) from a corner-influence matrix
      t <- colMeans(p); oe <- t["11"] - t["00"]
      c(W = unname((oe - (t["10"] - t["00"]) - (t["01"] - t["00"])) / oe),
        P = unname((t["01"] - t["00"]) / oe))
    }
    W_reps <- P_reps <- NULL
    if (reps == 1L) {
      phi <- .corner_fit(object, K, binY, covars)$phi
    } else {
      phis <- vector("list", reps); W_reps <- P_reps <- numeric(reps)
      for (r in seq_len(reps)) {
        set.seed(seed + r)
        phis[[r]] <- .corner_fit(object, K, binY, covars)$phi
        wp <- .gp_WP(phis[[r]]); W_reps[r] <- wp["W"]; P_reps[r] <- wp["P"]
      }
      phi <- Reduce(`+`, phis) / reps
    }
    th <- colMeans(phi)
    OE <- th["11"] - th["00"]; IDE <- th["10"] - th["00"]
    IIE <- th["01"] - th["00"]; R <- OE - IDE - IIE
    Pmed <- IIE / OE; W <- R / OE
    pOE <- phi[, "11"] - phi[, "00"]; pIDE <- phi[, "10"] - phi[, "00"]
    pIIE <- phi[, "01"] - phi[, "00"]; pR <- pOE - pIDE - pIIE
    se <- function(x) stats::sd(x) / sqrt(n); zc <- stats::qnorm(1 - (1 - ci_level) / 2)
    alpha <- 1 - ci_level
    seP <- se((pIIE - Pmed * pOE) / OE); seW <- se((pR - W * pOE) / OE)
    ## reps aggregation: add the residual fold Monte-Carlo variance of the averaged point.
    if (reps > 1L) {
      seW <- sqrt(seW^2 + stats::var(W_reps) / reps)
      seP <- sqrt(seP^2 + stats::var(P_reps) / reps)
    }
    ## analytic (symmetric Wald) intervals -- the default.
    W_ci     <- c(W - zc * seW, W + zc * seW)
    p_med_ci <- c(Pmed - zc * seP, Pmed + zc * seP)
    ## ---- bootstrap (near-null remedy): the analytic IF se for W and P_med is
    ## right-skewed and median-below the empirical SD, so the symmetric Wald CI
    ## under-covers (~0.85-0.90). W = R/OE and P_med = IIE/OE are ratios, so we use
    ## the tail-aware *percentile* bootstrap interval (resample rows, refit) rather
    ## than widening a symmetric se -- the latter mis-covers a skewed ratio. Cost is
    ## B (x reps) refits. seW/seP are still reported as bootstrap dispersion summaries.
    if (se_method == "bootstrap") {
      bsamp <- vapply(seq_len(B), function(b) {
        db <- object[sample.int(n, n, replace = TRUE), , drop = FALSE]
        pb <- if (reps == 1L) .corner_fit(db, K, binY, covars)$phi
              else Reduce(`+`, lapply(seq_len(reps),
                            function(r) .corner_fit(db, K, binY, covars)$phi)) / reps
        .gp_WP(pb)
      }, numeric(2))
      seW <- stats::sd(bsamp["W", ]); seP <- stats::sd(bsamp["P", ])
      W_ci     <- stats::quantile(bsamp["W", ], c(alpha / 2, 1 - alpha / 2), names = FALSE)
      p_med_ci <- stats::quantile(bsamp["P", ], c(alpha / 2, 1 - alpha / 2), names = FALSE)
    }
    z <- W / seW

    ## Fieller confidence set for P_med = IIE/OE. When the denominator OE is not
    ## significant the set is unbounded; the Wald interval understates this.
    fbounds <- numeric(0); ftype <- NA_character_
    if (isTRUE(fieller)) {
      VOE <- stats::var(pOE)/n; VIIE <- stats::var(pIIE)/n; COVio <- stats::cov(pIIE, pOE)/n
      fa <- OE^2 - zc^2 * VOE
      fb <- -2 * IIE * OE + 2 * zc^2 * COVio
      fc <- IIE^2 - zc^2 * VIIE
      disc <- fb^2 - 4 * fa * fc
      if (fa > 0) {
        if (disc >= 0) { fbounds <- sort((-fb + c(-1,1)*sqrt(disc))/(2*fa)); ftype <- "bounded" }
        else          { fbounds <- c(NA_real_, NA_real_); ftype <- "empty" }
      } else if (fa < 0) {
        if (disc >= 0) { fbounds <- sort((-fb + c(-1,1)*sqrt(disc))/(2*fa)); ftype <- "exclusive-unbounded" }
        else          { fbounds <- c(-Inf, Inf); ftype <- "all-real" }
      } else { # fa == 0 (linear)
        fbounds <- c(-Inf, Inf); ftype <- "all-real"
      }
    }

    GaugePmedResult(
      p_med = unname(Pmed), p_med_ci = unname(p_med_ci),
      p_med_fieller = unname(fbounds), fieller_type = ftype,
      W = unname(W), W_ci = unname(W_ci), W_se = unname(seW),
      W_p = unname(2 * stats::pnorm(-abs(z))),
      OE = unname(OE), IDE = unname(IDE), IIE = unname(IIE), R = unname(R),
      theta = th, method = "onestep-crossfit", n = as.integer(n),
      ci_level = ci_level, se_method = se_method, reps = as.integer(reps),
      call = match.call()
    )
  }

#' @export
S7::method(print, GaugePmedResult) <- function(x, ...) {
  lab <- if (identical(x@se_method, "bootstrap")) "percentile" else "Wald"
  cat("Gauge-calibrated proportion mediated (", x@method, ", n=", x@n, ")\n", sep = "")
  cat(sprintf("  P_med = %.3f  %s [%.3f, %.3f]\n", x@p_med, lab, x@p_med_ci[1], x@p_med_ci[2]))
  if (!is.na(x@fieller_type)) {
    fb <- x@p_med_fieller
    cat(switch(x@fieller_type,
      "bounded" = sprintf("    Fieller 95%% CI [%.3f, %.3f]\n", fb[1], fb[2]),
      "exclusive-unbounded" = sprintf(
        "    Fieller 95%% set (-Inf, %.3f] U [%.3f, Inf) -- OE not significant => UNBOUNDED\n",
        fb[1], fb[2]),
      "all-real" = "    Fieller 95% set = all of R (OE indistinguishable from 0)\n",
      "empty" = "    Fieller set empty (degenerate)\n", ""))
  }
  cat(sprintf("  W=R/OE = %.3f  %s [%.3f, %.3f]  (p=%.3g)\n",
              x@W, lab, x@W_ci[1], x@W_ci[2], x@W_p))
  cat(sprintf("  OE=%.3f  IDE=%.3f  IIE=%.3f  R=%.3f\n", x@OE, x@IDE, x@IIE, x@R))
  if (abs(x@W) > 0.1)
    cat("  ! |W| large: additive split unreliable; interpret P_med with care.\n")
  invisible(x)
}
