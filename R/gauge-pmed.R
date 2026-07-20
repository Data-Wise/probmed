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
#' @param W_ci_wald Numeric length-2: symmetric Wald interval for `W`, retained
#'   as the reference for the weak-identification diagnostic even under
#'   `se_method = "bootstrap"` (where `W_ci` holds the percentile interval).
#' @param weak_id Logical: weak-identification flag -- `TRUE` when the percentile
#'   CI for `W` is at least 3x wider than the symmetric Wald interval, the
#'   regime where the default Wald interval is most under-covering (Zhan 2026;
#'   see below for what the grid showed this does and does not mean). `NA` unless
#'   `se_method = "bootstrap"` (both intervals are required to compare), and
#'   also `NA` (not `FALSE`) in the degenerate case where the Wald interval has
#'   zero width. The 3x threshold is calibrated above the ~2x width gap that the
#'   anti-conservative Wald interval shows even under strong identification.
#'
#'   **An absent flag does not buy nominal coverage.** Independently of the flag,
#'   the symmetric Wald interval for `W` is anti-conservative *everywhere*:
#'   **~0.88** overall across the 16,000-rep coverage grid
#'   (`inst/sim/results/gauge_boot_coverage_nsim2000.csv`), never reaching the
#'   nominal 0.95 in any cell. So `weak_id = FALSE` cannot mean "this interval is
#'   trustworthy" -- at best it means "not the worst tail of an already
#'   sub-nominal default". If you need coverage you can rely on, use
#'   `se_method = "bootstrap"` -- but know what you are buying: the percentile
#'   interval covered **1.00** in all 8 cells of that grid, which is not "mildly"
#'   conservative, it is uncalibrated in the safe direction. Neither arm is
#'   nominal; pick the error you can live with.
#'
#'   **What the 48,000-rep validation grid established** (24 cells,
#'   `inst/sim/results/weakid_validation_cells.csv` + the two sweep files;
#'   adversarially re-verified from the raw reps): draws with `weak_id = TRUE`
#'   do have worse Wald coverage in **all 24 cells** (within-cell contrast
#'   -0.03 to -0.07 depending on aggregation). But that separation is
#'   **explained almost entirely by the Wald interval's own width** -- under a
#'   flexible width control the flag's coverage effect collapses to ~0 (-0.002,
#'   SE 0.003), and the bootstrap numerator of the ratio contributes nothing
#'   detectable. Read `weak_id` as a **self-contained proxy for "your analytic
#'   interval is narrow, and narrow analytic intervals under-cover"** -- useful
#'   precisely because a single fitted result offers no reference to judge
#'   "narrow" against -- not as evidence that the bootstrap-vs-Wald comparison
#'   detects a distinct pathology. The default threshold of 3 is a
#'   **convention**: thresholds 1.75-4 were not distinguishable in the grid, and
#'   all 24 cells share one simulation design, so the number has no claim to
#'   transfer across data-generating processes. `oe_regular` is the gate that
#'   carries information *beyond* interval width (see below); consult it and the
#'   Fieller set alongside this flag.
#' @param weak_id_ratio Numeric: ratio of the percentile `W`-interval width to the
#'   Wald `W`-interval width; `NA` if not computed.
#' @param oe_snr Numeric: signal-to-noise of the denominator, `|OE| / se(OE)`.
#'   Small values indicate `OE` near 0, where `W = R/OE` is non-regular.
#' @param oe_regular Logical: `TRUE` when `oe_snr >= 2`; `FALSE` flags a
#'   near-singular `OE` (Lin et al. 2026). Computed unconditionally from the
#'   analytic se(OE) regardless of `se_method` (only the accompanying
#'   `warning()` is gated on `se_method = "bootstrap"`). `NA` only in the
#'   degenerate case `se(OE) == 0`.
#'
#'   **What `FALSE` means -- width, not under-coverage.** In the 48,000-rep
#'   validation grid, draws flagged by this gate did **not** under-cover; they
#'   *over*-covered (+0.10, in 22/22 qualifying cells, stable across every
#'   design slice). The mechanism: as `OE` approaches 0 the ratio `W = R/OE`
#'   and its intervals blow up -- flagged Wald intervals were a median **21x
#'   wider than the true |W|** (vs 3x unflagged). Such an interval contains the
#'   truth and nearly everything else; it is untrustworthy because it is
#'   **uninformative**, not because it misses. Of the two gates this is the
#'   stronger diagnostic -- its effect survives controlling for interval width
#'   (+0.036, SE 0.003), unlike `weak_id`'s. Treat `oe_regular = FALSE` as
#'   "the point estimate and CI for `W` are numerically meaningless here", not
#'   as a coverage warning.
#' @param weak_id_ratio_threshold,oe_snr_threshold Numeric: the threshold values
#'   actually used for the `weak_id`/`oe_regular` gates on this result (see the
#'   identically-named arguments of [ward_residual()]). Stored so `print()` can
#'   report the threshold that was actually applied, not a hardcoded default.
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
    W_ci_wald = S7::new_property(class = S7::class_numeric, default = numeric(0)),
    weak_id = S7::new_property(class = S7::class_logical, default = NA),
    weak_id_ratio = S7::new_property(class = S7::class_numeric, default = NA_real_),
    oe_snr = S7::new_property(class = S7::class_numeric, default = NA_real_),
    oe_regular = S7::new_property(class = S7::class_logical, default = NA),
    weak_id_ratio_threshold = S7::new_property(class = S7::class_numeric, default = 3),
    oe_snr_threshold = S7::new_property(class = S7::class_numeric, default = 2),
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
#'   skewed ratio (widening a symmetric se fails to cover it); `W_se`/`p_med` se are
#'   still reported as bootstrap dispersion summaries. `reps > 1` and
#'   `se_method = "bootstrap"` compose. Bootstrap validity follows from the
#'   estimator being a Neyman-orthogonal cross-fit (DML-type) functional
#'   (Lin et al. 2026); it requires the ratio to be regular, i.e. `OE` bounded
#'   away from 0 -- see the Fieller diagnostic for the near-null case.
#' @param B Integer: number of bootstrap resamples when `se_method = "bootstrap"`
#'   (default `200`). Cost is `B` (x `reps`) refits.
#' @param a0,a1 Reference and comparison exposure levels (defaults `0`/`1`). `A`
#'   may use any two-level coding (factor, `{1,2}`, `{-1,1}`); it is recoded to
#'   the binary indicator `(A == a1)`. The gauge `W = R/OE` is invariant to
#'   swapping `a0`/`a1` (both `R` and `OE` flip sign). An `A` with **more than
#'   two** levels is an error, not a silent subset: restricting to `{a0,a1}`
#'   would shift the covariate-averaging population and estimate a different
#'   (sub-population) gauge. Filter to the two intended levels first.
#'   Multi-valued / continuous exposures are future work.
#' @param weak_id_ratio_threshold Numeric: the `weak_id` flag fires when the
#'   percentile `W`-CI is at least this many times wider than the Wald interval
#'   (default `3`). **A convention, not a calibrated optimum**: in the 48,000-rep
#'   validation grid (`inst/sim/hopper/run_weakid_validation.R`) no threshold in
#'   1.75-4 was distinguishable, and the grid spans a single simulation design --
#'   see the `weak_id` property documentation of [GaugePmedResult] for what the
#'   flag does and does not detect. Exposed so it can be tuned without editing
#'   the source.
#' @param oe_snr_threshold Numeric: the `oe_regular` guard flags a near-singular
#'   denominator when `|OE|/se(OE)` falls below this (default `2`). Validated in
#'   the same grid as marking **uninformatively wide** (over-covering) intervals,
#'   not under-coverage -- see the `oe_regular` property documentation of
#'   [GaugePmedResult].
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
                 se_method = c("analytic", "bootstrap"), B = 200L,
                 a0 = 0, a1 = 1, weak_id_ratio_threshold = 3,
                 oe_snr_threshold = 2, ...) {
    S7::S7_dispatch()
  })

#' @export
S7::method(ward_residual, S7::class_data.frame) <-
  function(object, covars = "C", K = 5L, ci_level = 0.95, seed = 1L, fieller = TRUE,
           reps = 1L, se_method = c("analytic", "bootstrap"), B = 200L,
           a0 = 0, a1 = 1, weak_id_ratio_threshold = 3, oe_snr_threshold = 2, ...) {
    stopifnot(all(c("A", "M", "Y") %in% names(object)), all(covars %in% names(object)))
    se_method <- match.arg(se_method)
    reps <- max(1L, as.integer(reps))
    ## ---- two-level exposure contrast (a0 reference, a1 comparison).
    ## A may use any two-level coding (factor, {1,2}, {-1,1}); it is recoded to the
    ## binary indicator (A == a1) for the corner machinery. A guard rejects >2 levels:
    ## subsetting a multi-valued A to {a0,a1} would silently shift the C-averaging
    ## population and estimate a restricted sub-population gauge, not the population
    ## gauge. Multi-valued / continuous A is future work (see manuscript sec-extensions).
    a_levels <- unique(object$A)
    if (!all(c(a0, a1) %in% a_levels))
      stop("a0 / a1 not found among the levels of A: have {",
           paste(a_levels, collapse = ", "), "}.", call. = FALSE)
    if (identical(a0, a1))
      stop("a0 and a1 must differ.", call. = FALSE)
    extra <- setdiff(a_levels, c(a0, a1))
    if (length(extra))
      stop("A has more than two levels (extra: {", paste(extra, collapse = ", "),
           "}). ward_residual contrasts exactly two levels; subsetting a multi-valued ",
           "A would change the estimand. Pre-filter to the two levels you intend, or ",
           "see the manuscript for the multi-valued generalization (not yet implemented).",
           call. = FALSE)
    object$A <- as.integer(object$A == a1)
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
    ## analytic (symmetric Wald) intervals -- the default, and the A2 reference.
    W_ci_wald <- c(W - zc * seW, W + zc * seW)  # retained for the weak-ID diagnostic
    W_ci     <- W_ci_wald
    p_med_ci <- c(Pmed - zc * seP, Pmed + zc * seP)
    ## A1 regularity precondition: OE signal-to-noise |OE| / se(OE). W = R/OE is a
    ## non-regular functional as OE -> 0, where bootstrap consistency for this
    ## Neyman-orthogonal cross-fit estimator fails (Lin et al. 2026). oe_snr below
    ## oe_snr_threshold flags OE statistically indistinguishable from 0.
    ## A gradient DGM sweep (A-effect scale s) shows weak_id_ratio decreasing in
    ## oe_snr ON AVERAGE -- ~3.4 at oe_snr~0.7 down to ~1.4 at oe_snr~12 (20 draws
    ## per point; inst/sim/pilot/weak_id_pilot.R) -- so the flag grades rather than
    ## merely detecting the near-null. Do NOT read exact crossover points off that
    ## sweep: the per-draw SD approaches the mean in the weak regime, so a
    ## single-draw sweep traces a clean curve by luck (an earlier version of this
    ## comment quoted such points as if stable; they did not reproduce).
    ## The 48k-rep validation grid (inst/sim/hopper/run_weakid_validation.R,
    ## 2026-07-16) settled the gates' division of labour, OPPOSITE to the naive
    ## reading: A2's coverage separation is ~97% explained by wid_wald alone
    ## (flag coef -0.051 -> -0.002 under a flexible width control), so A2 is a
    ## self-contained narrowness proxy, NOT evidence that CI divergence carries
    ## information beyond width -- Zhan (2026)'s premise did not survive for
    ## coverage. A1 is the gate with above-width content (+0.036 residual, SE
    ## 0.003) and flags OVER-covering, uninformatively wide intervals (median 21x
    ## wider than |truth|). t=3 is a convention (1.75-4 indistinguishable).
    seOE       <- se(pOE)
    oe_snr     <- unname(abs(OE) / seOE)
    ## NA-preserving: a degenerate se(OE) == 0 yields oe_snr = NaN, which must
    ## stay NA (undefined), not collapse to a confident FALSE via isTRUE(NA).
    oe_regular <- if (is.na(oe_snr)) NA else isTRUE(oe_snr >= oe_snr_threshold)
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

    ## A2 weak-identification flag: Wald-vs-percentile CI divergence is itself a
    ## weak-ID diagnostic (Zhan 2026). Computable only when both intervals exist,
    ## i.e. under se_method = "bootstrap" (W_ci is then the percentile interval,
    ## W_ci_wald the symmetric Wald). We flag when the percentile interval is >= 3x
    ## wider than the Wald. The symmetric Wald interval for W is *anti-conservative*
    ## by construction (a skewed ratio; ~0.85-0.90 coverage), so the percentile is
    ## routinely ~2x wider even under strong identification -- the 3x threshold is
    ## calibrated to fire only on the marginal inflation beyond that baseline, where
    ## the percentile interval starts tracking the ratio's exploding tail as OE->0.
    ## (An overlap metric is not separable here: it sits near 0.5 in both regimes.)
    weak_id <- NA; weak_id_ratio <- NA_real_
    if (se_method == "bootstrap") {
      wald_w <- W_ci_wald[2] - W_ci_wald[1]
      weak_id_ratio <- if (wald_w > 0) (W_ci[2] - W_ci[1]) / wald_w else NA_real_
      ## NA-preserving: a degenerate (zero-width) Wald interval leaves the
      ## ratio undefined -- isTRUE(NA >= x) would silently report "not weakly
      ## identified" (FALSE), which is a false confident claim.
      weak_id <- if (is.na(weak_id_ratio)) NA else isTRUE(weak_id_ratio >= weak_id_ratio_threshold)
    }
    ## Gate warnings (A1 fires only when a bootstrap CI for W is actually being
    ## reported; A2 whenever the divergence is detected). isFALSE()/isTRUE() (not
    ## bare negation) since oe_regular/weak_id may now be NA in degenerate cases.
    if (isFALSE(oe_regular) && se_method == "bootstrap")
      warning("Near-singular OE (|OE|/se = ", round(oe_snr, 2), " < ",
              oe_snr_threshold, "): W = R/OE is non-regular (Lin et al. 2026); ",
              "its point estimate and CI are numerically meaningless here -- ",
              "expect an interval many times wider than the estimand. Prefer ",
              "the Fieller set for the near-null case.", call. = FALSE)
    if (isTRUE(weak_id))
      warning("Weak-identification flag: the percentile CI for W is ",
              round(weak_id_ratio, 1), "x wider than the Wald interval (>= ",
              weak_id_ratio_threshold, "x). The default Wald CI is narrow here, ",
              "and narrow Wald intervals for W under-cover -- treat it with ",
              "suspicion (Zhan 2026).", call. = FALSE)

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
      W_ci_wald = unname(W_ci_wald), weak_id = weak_id,
      weak_id_ratio = weak_id_ratio, oe_snr = oe_snr, oe_regular = oe_regular,
      weak_id_ratio_threshold = weak_id_ratio_threshold,
      oe_snr_threshold = oe_snr_threshold,
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
  if (isTRUE(x@weak_id))
    cat(sprintf(paste0("  ! weak-ID: percentile CI for W is %.1fx wider than Wald",
                       " [%.3f, %.3f] (>= %gx); Wald CI likely too narrow (under-covers).\n"),
                x@weak_id_ratio, x@W_ci_wald[1], x@W_ci_wald[2], x@weak_id_ratio_threshold))
  if (isFALSE(x@oe_regular))
    cat(sprintf(paste0("  ! near-singular OE (|OE|/se = %.2f < %g): W = R/OE",
                       " non-regular; estimate and CI uninformative here.\n"),
                x@oe_snr, x@oe_snr_threshold))
  invisible(x)
}
