#' Sobol / Variance-Scale Proportion Mediated Result
#'
#' @description
#' S7 class for the **Sobol (functional-ANOVA) proportion mediated**
#' `P_med^{sigma2} = V_med / V_T` (Paper 3 of the P_med-for-modern-estimands
#' program). On the variance scale the proportion mediated is the first-order
#' Sobol sensitivity index of the outcome with respect to the mediated pathway:
#' the share of the total interventional outcome variance `V_T` explained by the
#' mediated component `V_med`. Unlike the additive split `TE = NDE + NIE`, the
#' variance decomposition `V_T = V_dir + V_med + V_int` always holds (the
#' interaction variance `V_int` is a separate non-negative term), so
#' `P_med^{sigma2}` is well defined even under treatment-by-mediator interaction.
#'
#' @details
#' Writing the four corner means as `theta(a, a') = E[Y(a, M(a'))]`, the contrasts
#' are the direct antecedent `Delta_d`, the mediated antecedent `Delta_m`, and the
#' interaction remainder `R = theta_11 - theta_10 - theta_01 + theta_00`. With
#' Bernoulli tilt variances `c_d = p_d(1 - p_d)` and `c_m = p_m(1 - p_m)`,
#' `V_dir = c_d Delta_d^2`, `V_med = c_m Delta_m^2`, `V_int = c_d c_m R^2`, and
#' `V_T = V_dir + V_med + V_int`. The estimator is the cross-fitted one-step
#' estimator built on the shared triply-robust corner pseudo-outcomes (the same
#' engine as [ward_residual()] and [incr_pmed()]); the delta method propagates the
#' corner influence functions through the quadratic `V`-components and the
#' ratio identity to a standard error for `P_med^{sigma2}`.
#'
#' **Boundary remedy.** Because `V_med = c_m Delta_m^2` is a squared quantity, its
#' influence function `phi_{V_med} = 2 c_m Delta_m phi_{Delta_m}` degenerates as
#' `Delta_m -> 0`, deflating the Wald standard error and collapsing coverage at
#' the null `V_med = 0`. The null is equivalent to `Delta_m = 0`, which has a
#' regular root-`n` influence function. Following Williamson et al. (2021,
#' *Biometrics*), the boundary is tested with a sample-split statistic that
#' estimates `Delta_m` on one half and its standard error on the independent other
#' half (so `z ~ N(0, 1)` even at the boundary). When the split test does not
#' reject `H0: V_med = 0` (`boundary == TRUE`), the symmetric Wald interval is
#' non-regular and a one-sided upper bound `[0, Pmed_upper]` (Procedure A) is
#' reported in `ci` instead; off-boundary, `ci` is the Wald interval.
#'
#' **Coverage caveat (Procedure A is _not_ uniformly valid).** The split *test*
#' restores level at the exact null, and Procedure A is conservative-to-nominal at
#' the exact boundary and in regular cells. But as a *pre-test* (gating) rule it is
#' only *pointwise* valid: across the near-null transition (a small but strictly
#' positive `V_med`, observed at intermediate-to-large `n`) the one-sided bound
#' `Pmed_upper = c_m (|Delta_m_hat| + z * se)^2 / V_T_hat` contracts at rate
#' `O_p(1/n)` while the truth stays fixed, so conditioning the report on
#' non-rejection routes downward-selected `Delta_m_hat` to a bound that can fall
#' below the truth -- a Leeb-Potscher post-selection pathology. In a registered
#' simulation grid the gated interval's coverage in the near-null cell was
#' conservative (~1.00) for `n <= 4000`, dropped to ~0.84 at `n = 8000`, and
#' recovered (~0.93) by `n = 16000` (a *transient* transition-zone dip, not an
#' intrinsically miscalibrated bound). Uniformly valid *interval* coverage through
#' the transition requires a sample-split CI (Procedure B) whose interval is not
#' conditioned on the test outcome; that interval is **not yet implemented here**.
#' Until it ships, treat `ci` at `boundary == TRUE` as a provisional one-sided bound,
#' prefer reporting `Pmed_upper` honestly as an upper bound (not a two-sided CI),
#' and consult `ci_wald` and the split-test fields (`vmed_split_p`,
#' `vmed_split_reject`) when characterising boundary behaviour.
#'
#' @param p_med Numeric: Sobol proportion mediated `V_med / V_T`.
#' @param se Numeric: standard error of `p_med` (delta-method, ratio identity).
#' @param ci Numeric length-2: **reported** (gated) confidence interval -- the
#'   Wald interval off-boundary, or the Procedure-A upper bound `[0, Pmed_upper]`
#'   at the boundary.
#' @param ci_wald Numeric length-2: ungated symmetric Wald interval (always the
#'   Wald interval, even at the boundary).
#' @param boundary Logical: `TRUE` when the split test does not reject
#'   `H0: V_med = 0` (the Wald CI is non-regular and `ci` is the upper bound).
#' @param Pmed_upper Numeric: one-sided Procedure-A upper bound for `p_med`.
#' @param S1_med Numeric: first-order Sobol index `V_med / V_T` (equals `p_med`).
#' @param ST_med Numeric: total Sobol index `(V_med + V_int) / V_T`.
#' @param Vd,Vm,Vdm,VT Numeric: direct, mediated, interaction, and total variance
#'   components.
#' @param Dm Numeric: mediated antecedent contrast `Delta_m`.
#' @param se_Dm Numeric: standard error of `Delta_m`.
#' @param vmed_split_p Numeric: p-value of the sample-split boundary test of
#'   `H0: V_med = 0`.
#' @param vmed_split_reject Integer: `1L` if the boundary test rejects, `0L`
#'   otherwise (`NA` when the test was not run).
#' @param theta Numeric length-4: corner means `theta(a, a')`.
#' @param method Character: estimation method.
#' @param n Integer: sample size.
#' @param ci_level Numeric: confidence level.
#' @param call Call: original call.
#'
#' @references
#' Williamson, B. D., Gilbert, P. B., Carone, M., & Simon, N. (2021).
#' Nonparametric variable importance assessment using machine learning
#' techniques. *Biometrics*, 77(1), 9--22.
#'
#' @export
SobolPmedResult <- S7::new_class(
  "SobolPmedResult", package = "probmed",
  properties = list(
    p_med = S7::class_numeric, se = S7::class_numeric,
    ci = S7::class_numeric, ci_wald = S7::class_numeric,
    boundary = S7::class_logical, Pmed_upper = S7::class_numeric,
    S1_med = S7::class_numeric, ST_med = S7::class_numeric,
    Vd = S7::class_numeric, Vm = S7::class_numeric,
    Vdm = S7::class_numeric, VT = S7::class_numeric,
    Dm = S7::class_numeric, se_Dm = S7::class_numeric,
    vmed_split_p = S7::new_property(class = S7::class_numeric, default = NA_real_),
    vmed_split_reject = S7::new_property(class = S7::class_integer, default = NA_integer_),
    theta = S7::class_numeric, method = S7::class_character,
    n = S7::class_integer, ci_level = S7::class_numeric,
    call = S7::new_property(class = S7::class_any, default = NULL)
  ),
  validator = function(self) {
    if (length(self@VT) && self@VT <= 0)
      "VT (total variance) must be strictly positive"
  }
)

#' Map corner means to the Sobol variance decomposition
#'
#' @description
#' Closed-form map from the four corner means `theta(a, a')` to the Sobol /
#' functional-ANOVA variance components and the variance-scale proportion
#' mediated. Used internally by [sobol_pmed()] and exported as a helper for
#' computing the analytic (population) truth from known corner means.
#'
#' @param theta Numeric length-4: corner means named `"11"`, `"10"`, `"01"`,
#'   `"00"` (i.e. `theta(a, a')` for `(a, a')` in the four corners).
#' @param pd,pm Numeric: tilt probabilities for the direct and mediated pathways
#'   (Bernoulli design points); defaults `0.5`.
#'
#' @return A named numeric vector with elements `Vd`, `Vm`, `Vdm`, `VT`,
#'   `Pmed_sobol` (= `Vm / VT`), `S1_med` (= `Vm / VT`),
#'   `ST_med` (= `(Vm + Vdm) / VT`), and `Rint` (the interaction remainder).
#'
#' @examples
#' th <- c("11" = 1.32, "10" = 0.5, "01" = 0.42, "00" = 0)
#' sobol_from_theta(th)
#'
#' @export
sobol_from_theta <- function(theta, pd = 0.5, pm = 0.5) {
  t11 <- theta[["11"]]; t10 <- theta[["10"]]; t01 <- theta[["01"]]; t00 <- theta[["00"]]
  Dd <- (1 - pm) * (t10 - t00) + pm * (t11 - t01)
  Dm <- (1 - pd) * (t01 - t00) + pd * (t11 - t10)
  R <- t11 - t10 - t01 + t00
  cd <- pd * (1 - pd); cm <- pm * (1 - pm)
  Vd <- cd * Dd^2; Vm <- cm * Dm^2; Vdm <- cd * cm * R^2; VT <- Vd + Vm + Vdm
  c(Vd = unname(Vd), Vm = unname(Vm), Vdm = unname(Vdm), VT = unname(VT),
    Pmed_sobol = unname(Vm / VT), S1_med = unname(Vm / VT),
    ST_med = unname((Vm + Vdm) / VT), Rint = unname(R))
}

# internal: core Sobol fit on a single data frame. Returns a plain list (NOT an
# S7 object) so the boundary sample-split can recurse into it directly without
# round-tripping through the S7 generic. RNG contract (must match the sim source
# of truth byte-for-byte): set.seed(seed) here; the first RNG-consuming call is
# .corner_fit() (fold draw); then the split's sample.int(); the two split halves
# recurse with boundary_test = "none" (no further split, no extra RNG beyond
# their own corner draws). binY is fixed FALSE: the Sobol scale is defined for a
# continuous (Gaussian) outcome.
.sobol_fit <- function(d, pd = 0.5, pm = 0.5, covars = "C", K = 5L, seed = 1L,
                       level = 0.95, warn_boundary = TRUE,
                       boundary_test = c("split", "plugin", "none")) {
  boundary_test <- match.arg(boundary_test)
  set.seed(seed); n <- nrow(d)
  phi <- .corner_fit(d, K, binY = FALSE, covars)$phi
  th <- colMeans(phi)
  cd <- pd * (1 - pd); cm <- pm * (1 - pm); cdm <- cd * cm
  Dd <- (1 - pm) * (th["10"] - th["00"]) + pm * (th["11"] - th["01"])
  Dm <- (1 - pd) * (th["01"] - th["00"]) + pd * (th["11"] - th["10"])
  Rg <- th["11"] - th["10"] - th["01"] + th["00"]
  Vd <- cd * Dd^2; Vm <- cm * Dm^2; Vdm <- cdm * Rg^2; VT <- Vd + Vm + Vdm; P <- Vm / VT
  ## per-observation EIFs (centered corner phis), delta method through V-components
  cphi <- sweep(phi, 2, th)
  pDd <- (1 - pm) * (cphi[, "10"] - cphi[, "00"]) + pm * (cphi[, "11"] - cphi[, "01"])
  pDm <- (1 - pd) * (cphi[, "01"] - cphi[, "00"]) + pd * (cphi[, "11"] - cphi[, "10"])
  pR  <- cphi[, "11"] - cphi[, "10"] - cphi[, "01"] + cphi[, "00"]
  pVd <- 2 * cd * Dd * pDd; pVm <- 2 * cm * Dm * pDm; pVdm <- 2 * cdm * Rg * pR
  pVT <- pVd + pVm + pVdm
  pP  <- (pVm - as.numeric(P) * pVT) / as.numeric(VT)        # ratio-identity IF
  se  <- stats::sd(pP) / sqrt(n); zc <- stats::qnorm(1 - (1 - level) / 2)
  ## ---- V_med = 0 boundary test (Williamson 2021 VIM, Delta_m reparametrization) ----
  ## V_med = c_m Delta_m^2 with c_m = p_m(1 - p_m) > 0, so V_med = 0 <=> Delta_m = 0.
  ## The squaring degenerates phi_{V_med} = 2 c_m Delta_m phi_{Delta_m} -> 0 at the
  ## boundary (deflating the Wald SE for P); test on the regular contrast Delta_m.
  se_Dm <- stats::sd(pDm) / sqrt(n)
  z_Dm  <- as.numeric(Dm) / se_Dm
  p_plugin <- 2 * stats::pnorm(-abs(z_Dm))                   # plug-in test (DIAGNOSTIC ONLY)
  if (boundary_test == "split" && n >= 8L) {
    ## sample-split decorrelation: estimate Delta_m on one half, its SE on the
    ## INDEPENDENT other half, so z ~ N(0, 1) exactly even at the boundary.
    h  <- sample.int(n, n %/% 2L)
    f1 <- .sobol_fit(d[h, , drop = FALSE],  pd, pm, covars, K, seed = seed + 1L, level, FALSE, "none")
    f2 <- .sobol_fit(d[-h, , drop = FALSE], pd, pm, covars, K, seed = seed + 2L, level, FALSE, "none")
    z1 <- f1$Dm / f2$se_Dm; z2 <- f2$Dm / f1$se_Dm           # decorrelated (independent halves)
    vmed_split_reject <- as.integer(abs(z1) > zc && abs(z2) > zc)  # both orderings (conservative)
    vmed_split_p <- max(2 * stats::pnorm(-abs(z1)), 2 * stats::pnorm(-abs(z2)))
  } else if (boundary_test == "plugin") {
    vmed_split_reject <- as.integer(p_plugin < (1 - level)); vmed_split_p <- p_plugin
  } else {                                                    # "none" (inside the split recursion)
    vmed_split_reject <- NA_integer_; vmed_split_p <- NA_real_
  }
  ## ---- Procedure A: gate the CI; at the boundary report a one-sided upper bound ----
  boundary <- isTRUE(vmed_split_reject == 0L)
  Pmed_upper <- as.numeric(cm * (abs(as.numeric(Dm)) + stats::qnorm(level) * se_Dm)^2 / as.numeric(VT))
  ci_wald <- unname(c(P - zc * se, P + zc * se))
  ci <- if (boundary) c(0, Pmed_upper) else ci_wald          # REPORTED interval (gated)
  if (warn_boundary && boundary)
    warning("sobol_pmed: H0 V_med=0 not rejected (split test, p=", signif(vmed_split_p, 2),
            "); the symmetric Wald CI is non-regular at the boundary. Reporting the one-sided ",
            "upper bound [0, ", signif(Pmed_upper, 3), "] (Procedure A) instead. ",
            "Note: this gated bound is NOT uniformly valid across the near-null transition ",
            "(pre-test under-coverage at intermediate-to-large n); treat it as a provisional ",
            "one-sided upper bound. A uniformly valid sample-split CI (Procedure B) is not yet ",
            "implemented -- see ?SobolPmedResult 'Coverage caveat'.", call. = FALSE)
  list(P_med_sobol = unname(P), se = unname(se), ci = unname(ci), ci_wald = ci_wald,
       boundary = boundary, Pmed_upper = Pmed_upper,
       S1_med = unname(Vm / VT), ST_med = unname((Vm + Vdm) / VT),
       Vd = unname(Vd), Vm = unname(Vm), Vdm = unname(Vdm), VT = unname(VT),
       Dm = unname(Dm), se_Dm = unname(se_Dm),
       vmed_split_p = unname(vmed_split_p), vmed_split_reject = unname(vmed_split_reject),
       p_plugin = unname(p_plugin), theta = round(th, 4), n = n)
}

#' Sobol / Variance-Scale Proportion Mediated
#'
#' @description
#' Estimate the **Sobol proportion mediated** `P_med^{sigma2} = V_med / V_T` -- the
#' share of the total interventional outcome variance carried by the mediated
#' pathway, i.e. the first-order Sobol (functional-ANOVA) sensitivity index of the
#' outcome with respect to the mediator. The estimator is cross-fitted and built
#' on the same triply-robust corner pseudo-outcomes as [ward_residual()] and
#' [incr_pmed()]; the standard error uses the delta method through the quadratic
#' variance components and the ratio identity. See [SobolPmedResult] for the
#' variance decomposition and the boundary remedy.
#'
#' @details
#' At the null `V_med = 0` (equivalently `Delta_m = 0`) the squared-variance
#' influence function degenerates and the symmetric Wald interval is non-regular.
#' A Williamson et al. (2021) sample-split test on the regular contrast `Delta_m`
#' decides the boundary; when it does not reject, the reported interval `ci`
#' switches to the one-sided Procedure-A upper bound `[0, Pmed_upper]`. The
#' boundary machinery recurses through a plain internal fitter (not this generic),
#' so the point estimate, Wald interval, and fold draw are deterministic in
#' `seed` and unaffected by the test. **Caveat:** this gating (Procedure A) is
#' pointwise- but not uniformly-valid across the near-null transition; the gated
#' interval can under-cover at intermediate-to-large `n` and a uniformly valid
#' sample-split CI (Procedure B) is not yet implemented. See the **Coverage caveat**
#' in [SobolPmedResult].
#'
#' @param object A `data.frame` with columns `A` (binary treatment), `M`
#'   (mediator), `Y` (continuous outcome), and the covariates named in `covars`.
#' @param pd,pm Numeric: tilt probabilities for the direct and mediated pathways
#'   (Bernoulli design points); defaults `0.5`.
#' @param covars Character vector of covariate column names. Default `"C"`.
#' @param K Integer: number of cross-fitting folds (default `5`).
#' @param seed Integer: RNG seed for fold assignment and the boundary split
#'   (default `1`).
#' @param ci_level Numeric: confidence level (default `0.95`).
#' @param warn_boundary Logical: warn when the boundary test does not reject and
#'   the reported interval falls back to the Procedure-A upper bound
#'   (default `TRUE`).
#' @param boundary_test Character: boundary test for `H0: V_med = 0`. `"split"`
#'   (default) is the Williamson sample-split test on `Delta_m`; `"plugin"` is the
#'   (over-rejecting) plug-in Wald test, diagnostic only.
#' @param ... Unused.
#'
#' @return A [SobolPmedResult] object.
#'
#' @references
#' Williamson, B. D., Gilbert, P. B., Carone, M., & Simon, N. (2021).
#' Nonparametric variable importance assessment using machine learning
#' techniques. *Biometrics*, 77(1), 9--22.
#'
#' @examples
#' set.seed(1)
#' n <- 1500; C <- rnorm(n)
#' A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
#' M <- 0.6 * A + 0.4 * C + rnorm(n)
#' # strong mediation (off-boundary): Wald interval reported
#' Y <- 0.5 * A + 0.7 * M + 0.3 * C + rnorm(n)
#' sobol_pmed(data.frame(A, M, Y, C))
#'
#' @export
sobol_pmed <- S7::new_generic(
  "sobol_pmed", dispatch_args = "object",
  fun = function(object, pd = 0.5, pm = 0.5, covars = "C", K = 5L, seed = 1L,
                 ci_level = 0.95, warn_boundary = TRUE,
                 boundary_test = c("split", "plugin", "none"), ...) {
    S7::S7_dispatch()
  })

#' @export
S7::method(sobol_pmed, S7::class_data.frame) <-
  function(object, pd = 0.5, pm = 0.5, covars = "C", K = 5L, seed = 1L,
           ci_level = 0.95, warn_boundary = TRUE,
           boundary_test = c("split", "plugin", "none"), ...) {
    stopifnot(all(c("A", "M", "Y") %in% names(object)), all(covars %in% names(object)))
    boundary_test <- match.arg(boundary_test)
    f <- .sobol_fit(object, pd = pd, pm = pm, covars = covars, K = K, seed = seed,
                    level = ci_level, warn_boundary = warn_boundary,
                    boundary_test = boundary_test)
    SobolPmedResult(
      p_med = f$P_med_sobol, se = f$se, ci = f$ci, ci_wald = f$ci_wald,
      boundary = f$boundary, Pmed_upper = f$Pmed_upper,
      S1_med = f$S1_med, ST_med = f$ST_med,
      Vd = f$Vd, Vm = f$Vm, Vdm = f$Vdm, VT = f$VT,
      Dm = f$Dm, se_Dm = f$se_Dm,
      vmed_split_p = f$vmed_split_p,
      vmed_split_reject = if (is.na(f$vmed_split_reject)) NA_integer_ else as.integer(f$vmed_split_reject),
      theta = f$theta, method = "onestep-crossfit", n = as.integer(f$n),
      ci_level = ci_level, call = match.call()
    )
  }

#' @export
S7::method(print, SobolPmedResult) <- function(x, ...) {
  cat("Sobol / variance-scale proportion mediated P_med^sigma2 (",
      x@method, ", n=", x@n, ")\n", sep = "")
  cat(sprintf("  P_med = %.3f  se = %.3f\n", x@p_med, x@se))
  if (isTRUE(x@boundary)) {
    cat(sprintf("  CI (Procedure A, boundary) [0, %.3f]   (Wald [%.3f, %.3f])\n",
                x@Pmed_upper, x@ci_wald[1], x@ci_wald[2]))
    cat(sprintf("  ! H0 V_med=0 not rejected (split p=%.3g): Wald CI non-regular at boundary.\n",
                x@vmed_split_p))
  } else {
    cat(sprintf("  Wald %g%% CI [%.3f, %.3f]\n", 100 * x@ci_level, x@ci[1], x@ci[2]))
    if (!is.na(x@vmed_split_p))
      cat(sprintf("  H0 V_med=0 rejected (split p=%.3g).\n", x@vmed_split_p))
  }
  cat(sprintf("  S1_med = %.3f  ST_med = %.3f\n", x@S1_med, x@ST_med))
  cat(sprintf("  Vd=%.4f  Vm=%.4f  Vint=%.4f  VT=%.4f\n", x@Vd, x@Vm, x@Vdm, x@VT))
  invisible(x)
}
