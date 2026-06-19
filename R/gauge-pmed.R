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
    call = S7::new_property(class = S7::class_any, default = NULL)
  ),
  validator = function(self) {
    if (length(self@OE) && abs(self@OE) < 1e-8)
      warning("OE near 0: W = R/OE is unstable; report unnormalized R.")
  }
)

# internal: cross-fit one-step influence matrix for the four corner means
.gauge_fit <- function(d, K, binY, covars) {
  cf <- paste(covars, collapse = " + ")
  f_pi <- stats::as.formula(paste("A ~", cf))
  f_q  <- stats::as.formula(paste("A ~ M +", cf))
  f_om <- stats::as.formula(paste("Y ~ A * M +", cf))
  n <- nrow(d); folds <- sample(rep(1:K, length.out = n))
  nm <- c("11", "10", "01", "00"); cor <- list(c(1,1), c(1,0), c(0,1), c(0,0))
  phi <- matrix(0, n, 4, dimnames = list(NULL, nm))
  for (k in 1:K) {
    tr <- d[folds != k, , drop = FALSE]
    te_i <- which(folds == k); te <- d[te_i, , drop = FALSE]
    pim <- stats::glm(f_pi, data = tr, family = stats::binomial())
    qm  <- stats::glm(f_q,  data = tr, family = stats::binomial())
    om  <- stats::glm(f_om, data = tr,
                      family = if (binY) stats::binomial() else stats::gaussian())
    p1 <- stats::predict(pim, newdata = te, type = "response")
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
  phi
}

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
ward_residual <- S7::new_generic("ward_residual", dispatch_args = "object")

#' @export
S7::method(ward_residual, S7::class_data.frame) <-
  function(object, covars = "C", K = 5L, ci_level = 0.95, seed = 1L, fieller = TRUE, ...) {
    stopifnot(all(c("A", "M", "Y") %in% names(object)), all(covars %in% names(object)))
    set.seed(seed)
    binY <- all(object$Y %in% 0:1)
    phi <- .gauge_fit(object, K, binY, covars); n <- nrow(object)
    th <- colMeans(phi)
    OE <- th["11"] - th["00"]; IDE <- th["10"] - th["00"]
    IIE <- th["01"] - th["00"]; R <- OE - IDE - IIE
    Pmed <- IIE / OE; W <- R / OE
    pOE <- phi[, "11"] - phi[, "00"]; pIDE <- phi[, "10"] - phi[, "00"]
    pIIE <- phi[, "01"] - phi[, "00"]; pR <- pOE - pIDE - pIIE
    se <- function(x) stats::sd(x) / sqrt(n); zc <- stats::qnorm(1 - (1 - ci_level) / 2)
    seP <- se((pIIE - Pmed * pOE) / OE); seW <- se((pR - W * pOE) / OE); z <- W / seW

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
      p_med = unname(Pmed), p_med_ci = unname(c(Pmed - zc * seP, Pmed + zc * seP)),
      p_med_fieller = unname(fbounds), fieller_type = ftype,
      W = unname(W), W_ci = unname(c(W - zc * seW, W + zc * seW)), W_se = unname(seW),
      W_p = unname(2 * stats::pnorm(-abs(z))),
      OE = unname(OE), IDE = unname(IDE), IIE = unname(IIE), R = unname(R),
      theta = th, method = "onestep-crossfit", n = as.integer(n),
      ci_level = ci_level, call = match.call()
    )
  }

#' @export
S7::method(print, GaugePmedResult) <- function(x, ...) {
  cat("Gauge-calibrated proportion mediated (", x@method, ", n=", x@n, ")\n", sep = "")
  cat(sprintf("  P_med = %.3f  Wald [%.3f, %.3f]\n", x@p_med, x@p_med_ci[1], x@p_med_ci[2]))
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
  cat(sprintf("  W=R/OE = %.3f  [%.3f, %.3f]  (p=%.3g)\n", x@W, x@W_ci[1], x@W_ci[2], x@W_p))
  cat(sprintf("  OE=%.3f  IDE=%.3f  IIE=%.3f  R=%.3f\n", x@OE, x@IDE, x@IIE, x@R))
  if (abs(x@W) > 0.1)
    cat("  ! |W| large: additive split unreliable; interpret P_med with care.\n")
  invisible(x)
}
