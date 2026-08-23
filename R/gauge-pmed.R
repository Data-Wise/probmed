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
#' @param weak_id Logical: `TRUE` when the percentile CI for `W` is at least
#'   `weak_id_ratio_threshold` (default 3) times wider than the symmetric Wald
#'   interval. `NA` unless `se_method = "bootstrap"` (both intervals are
#'   required), and `NA` (not `FALSE`) when the Wald interval has zero width.
#'   The ratio is a scale-only adaptation of the bootstrap-vs-asymptotic-
#'   discrepancy idea in Zhan (2026); it is not Zhan's statistic, which is a
#'   Kolmogorov-Smirnov distance for linear IV.
#'
#'   **What it measures** (2026-08-22, corrected corner weights;
#'   `docs/specs/FINDINGS-2026-08-22-postfix-bootstrap-gates.md`): a
#'   bootstrap-vs-Wald width discrepancy, a symptom of bootstrap instability
#'   near `OE = 0` -- **not a coverage diagnostic**. Against exact truth (4
#'   cells x 200 datasets, n = 800) draws flagged at *any* threshold from 1.2
#'   to 3 have the same Wald coverage as unflagged draws (0.94-0.95 vs
#'   0.94-0.97). In the regular cells the ratio is 1.02 (median), 1.22 (q95),
#'   1.49 (q99); near the null it tracks `oe_regular = FALSE` (P(flag |
#'   irregular) 0.82-0.88 at thresholds >= 1.5 vs 0.15-0.40 given regular) and
#'   at 3x adds 2 of the 98 near-null draws that pass A1. No `warning()` is
#'   issued on it; `print()` reports the ratio when it exceeds the threshold.
#'   `oe_regular` is the gate to consult, with the Fieller set.
#'
#'   **History.** The flag shipped with its 3x threshold "calibrated above the
#'   ~2x width gap of the anti-conservative Wald interval" and was validated on a
#'   48,000-rep grid as a "narrowness proxy" whose coverage separation was
#'   explained by Wald width (`inst/sim/results/weakid_validation_cells.csv`).
#'   All of that measured the estimator *before* the corner-EIF weight fix
#'   (NEWS 0.3.0.9000): the ~30% se shortfall and the 2x gap were the bug. The
#'   threshold stays at 3 -- with no signal to calibrate to, a new number would
#'   imply a calibration that does not exist -- and the field is kept because a
#'   user running the bootstrap may want to see its instability.
#' @param weak_id_ratio Numeric: ratio of the percentile `W`-interval width to the
#'   Wald `W`-interval width; `NA` if not computed.
#' @param oe_snr Numeric: signal-to-noise of the denominator, `|OE| / se(OE)`.
#'   Small values indicate `OE` near 0, where `W = R/OE` is non-regular.
#' @param oe_regular Logical: `TRUE` when `oe_snr >= 2`; `FALSE` flags a
#'   near-singular `OE` -- the denominator of the ratio `W = R/OE`, whose
#'   interval estimation is the classical Fieller (1954) problem. Computed unconditionally from the
#'   analytic se(OE) regardless of `se_method` (only the accompanying
#'   `warning()` is gated on `se_method = "bootstrap"`). `NA` only in the
#'   degenerate case `se(OE) == 0`.
#'
#'   **What `FALSE` means -- width, not under-coverage.** (Re-measured with the
#'   corrected corner weights, 2026-08-22, same reading: in the near-null cell
#'   flagged draws cover 0.99 at a median Wald width 8.7x |W|, while the
#'   near-null draws that *pass* the gate cover 0.90 at 4.1x --
#'   `docs/specs/FINDINGS-2026-08-22-postfix-bootstrap-gates.md`.) In the 48,000-rep
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
#'   variance of the averaged point. With the corrected corner weights (see
#'   `se_method`) that component is small -- n = 800, continuous: empSD 0.063 at
#'   `reps = 1` vs 0.061 at `reps = 10`, oracle 0.056 -- so the default is
#'   usually enough.
#' @param se_method Character: `"analytic"` (default) or `"bootstrap"`. Under
#'   `"analytic"` the se is the influence-function se of the cross-fit one-step
#'   estimator and the CI is symmetric Wald. Measured against exact truth with
#'   the corrected corner weights (`inst/sim/se_shipped_check.R`, 250 datasets
#'   per cell, n = 800): Wald coverage 0.936 (continuous `Y`, strong
#'   identification), 0.972 (binary `Y`), 0.912 (intermediate `OE`), with the
#'   se's dispersion across datasets matching an oracle that knows the nuisances
#'   (CV 0.23 vs 0.18). Under `"bootstrap"` the CIs for `W` and `P_med` are
#'   **percentile** intervals of a nonparametric bootstrap that resamples rows
#'   and refits the whole cross-fit (a fresh partition each resample);
#'   `W_se`/`p_med` se are then bootstrap dispersion summaries. Re-measured
#'   with the corrected weights (`inst/sim/boot_gates_check.R`, B = 200, 200
#'   datasets per cell): percentile coverage 0.905 (continuous) / 0.930
#'   (binary) / 0.920 (intermediate) at the same width as the Wald interval,
#'   and a bootstrap SD that is unstable for binary `Y` (CV 2.2 across
#'   datasets vs 0.33 analytic) -- an option for users who want a refit
#'   bootstrap, not a remedy for anything, at `B` (x `reps`) refits. The
#'   bootstrap-consistency result of Lin and Han (2026) for cross-fit DML
#'   functionals holds the nuisances fixed and does **not** cover this
#'   refit-per-resample scheme. `reps > 1` and `se_method = "bootstrap"`
#'   compose. Near the null (`OE` not bounded away from 0) the ratio is
#'   non-regular for every arm -- even an oracle's se explodes -- use the
#'   Fieller set.
#'
#'   **History.** Every coverage figure this package reported before
#'   2026-08-22 -- the 2,000-rep coverage grid
#'   (`inst/sim/results/gauge_boot_coverage_nsim2000.csv`: Wald ~0.85-0.90,
#'   percentile ~1.00), the 48,000-rep gate-validation grid, and the variance
#'   decomposition (`inst/sim/results/phi_decomp_summary.csv`;
#'   `docs/specs/FINDINGS-2026-08-22-phi-decomposition.md`) -- was produced
#'   with a corner-EIF weight bug: `ifelse()` on a scalar test returned the
#'   **first fold row's** propensity (and mediator-density proxy) for every row,
#'   so the inverse-probability weights were one arbitrary constant per fold
#'   (NEWS, 0.3.0.9000). The "fold-split noise" those grids diagnosed was this
#'   bug. They are superseded. The bootstrap arm and both gates were
#'   re-measured with the fix the same day
#'   (`docs/specs/FINDINGS-2026-08-22-postfix-bootstrap-gates.md`,
#'   `inst/sim/results/boot_gates_postfix.csv`): see above, and the `weak_id` /
#'   `oe_regular` property docs of [GaugePmedResult].
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
#'   (default `3`). **Not calibrated to coverage, and cannot be** (2026-08-22,
#'   corrected weights): flagged and unflagged draws cover alike at every
#'   threshold from 1.2 to 3 -- see the `weak_id` property of [GaugePmedResult].
#'   Left at 3 so the flag stays rare outside the near-null regime (0% of
#'   regular-cell draws, 19% near the null). Exposed so it can be tuned without
#'   editing the source.
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
    ## (With the corrected corner weights -- 2026-08-22, see .corner_phi -- this
    ## cross-fit IF se is calibrated at every reps: n = 800, 250 datasets, Wald
    ## coverage 0.936 continuous / 0.972 binary at reps = 1, 0.936 / 0.968 at
    ## reps = 10, se dispersion matching an oracle's; inst/sim/se_shipped_check.R,
    ## inst/sim/results/se_shipped_postfix.csv.
    ## A full-sample-nuisance se was built and measured on this branch and dropped:
    ## no better once the weights were right.)
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
    ## History: a gradient sweep (inst/sim/pilot/weak_id_pilot.R) and the 48k-rep
    ## grid (inst/sim/hopper/run_weakid_validation.R, 2026-07-16) characterized A2
    ## as a narrowness proxy and A1 as the gate with above-width content that flags
    ## OVER-covering intervals (median 21x wider than |truth|). Both were measured
    ## on the pre-fix corner weights. The post-fix re-measurement (2026-08-22,
    ## docs/specs/FINDINGS-2026-08-22-postfix-bootstrap-gates.md) keeps A1's
    ## reading (flagged near-null draws cover 0.99 at 8.7x |W|) and finds A2
    ## uninformative for coverage at every threshold.
    seOE       <- se(pOE)
    oe_snr     <- unname(abs(OE) / seOE)
    ## NA-preserving: a degenerate se(OE) == 0 yields oe_snr = NaN, which must
    ## stay NA (undefined), not collapse to a confident FALSE via isTRUE(NA).
    oe_regular <- if (is.na(oe_snr)) NA else isTRUE(oe_snr >= oe_snr_threshold)
    ## ---- bootstrap: W = R/OE and P_med = IIE/OE are ratios, so the alternative
    ## arm is the tail-aware *percentile* bootstrap interval (resample rows, refit)
    ## rather than a widened symmetric se. Cost is B (x reps) refits. seW/seP are
    ## then reported as bootstrap dispersion summaries. (The "Wald under-covers
    ## ~0.85-0.90" motivation this arm was added under came from the pre-fix weight
    ## bug. Post-fix, 2026-08-22: percentile coverage 0.905 / 0.930 / 0.920 in the
    ## regular cells at the same width as Wald, bootstrap SD unstable for binary Y
    ## (CV 2.2) -- an option, not a remedy; see the se_method docs.)
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

    ## A2: the percentile/Wald width ratio, a scale-only adaptation of Zhan (2026)'s
    ## bootstrap-vs-asymptotic-discrepancy idea (his statistic is a KS distance in
    ## linear IV; this is not it). Computable only under se_method = "bootstrap"
    ## (W_ci is then the percentile interval, W_ci_wald the symmetric Wald).
    ## Post-fix (2026-08-22, FINDINGS-2026-08-22-postfix-bootstrap-gates) it is a
    ## bootstrap-instability discrepancy, NOT a coverage diagnostic: flagged draws
    ## cover like unflagged ones at every threshold 1.2-3; regular-cell ratio
    ## median 1.02, q99 1.49; near-null median 1.74, tracking oe_regular = FALSE.
    ## The 3x threshold was set pre-fix "above the ~2x baseline" -- a baseline
    ## that was the weight bug's ~30% se shortfall -- and is kept unchanged; no
    ## warning is raised on it (print() shows the ratio). See the weak_id docs.
    weak_id <- NA; weak_id_ratio <- NA_real_
    if (se_method == "bootstrap") {
      wald_w <- W_ci_wald[2] - W_ci_wald[1]
      weak_id_ratio <- if (wald_w > 0) (W_ci[2] - W_ci[1]) / wald_w else NA_real_
      ## NA-preserving: a degenerate (zero-width) Wald interval leaves the
      ## ratio undefined -- isTRUE(NA >= x) would silently report "not weakly
      ## identified" (FALSE), which is a false confident claim.
      weak_id <- if (is.na(weak_id_ratio)) NA else isTRUE(weak_id_ratio >= weak_id_ratio_threshold)
    }
    ## Gate warning: A1 fires only when a bootstrap CI for W is actually being
    ## reported. A2 (weak_id) no longer warns (2026-08-22: not a coverage
    ## diagnostic); print() shows its ratio. isFALSE() (not bare negation) since
    ## oe_regular may be NA in degenerate cases.
    if (isFALSE(oe_regular) && se_method == "bootstrap")
      warning("Near-singular OE (|OE|/se = ", round(oe_snr, 2), " < ",
              oe_snr_threshold, "): W = R/OE is non-regular; ",
              "its point estimate and CI are numerically meaningless here -- ",
              "expect an interval many times wider than the estimand. Prefer ",
              "the Fieller set for the near-null case.", call. = FALSE)

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
    cat(sprintf(paste0("  ! bootstrap/Wald width discrepancy: percentile CI for W is %.1fx wider than Wald",
                       " [%.3f, %.3f] (>= %gx); a bootstrap-instability symptom, not a coverage diagnostic -- see oe_regular.\n"),
                x@weak_id_ratio, x@W_ci_wald[1], x@W_ci_wald[2], x@weak_id_ratio_threshold))
  if (isFALSE(x@oe_regular))
    cat(sprintf(paste0("  ! near-singular OE (|OE|/se = %.2f < %g): W = R/OE",
                       " non-regular; estimate and CI uninformative here.\n"),
                x@oe_snr, x@oe_snr_threshold))
  invisible(x)
}
