#' P_med Sensitivity Result (M--Y Confounding)
#'
#' @description
#' S7 result for the sensitivity of an indirect-to-total effect ratio to an
#' additive bias in its numerator. It stores the original ratio, the ratio on
#' the bias grid, and the bias needed to make the indirect effect zero or to
#' reach a specified ratio threshold.
#'
#' @param p_med Numeric: the point estimate `indirect / total` (b = 0).
#' @param indirect Numeric: the indirect (mediated) component, the numerator.
#' @param total Numeric: the total effect, the denominator.
#' @param bias_grid Numeric: the grid of numerator biases `b` evaluated.
#' @param p_med_bias Numeric: `P_med` at each `b`, i.e. `(indirect + b)/total`.
#' @param tipping_indirect Numeric: bias `b` at which `indirect + b = 0`
#'   (the indirect effect vanishes). `NA` if `total == 0`.
#' @param tipping_threshold Numeric: bias `b` at which `P_med` hits
#'   `threshold`; `NA` if not requested.
#' @param threshold Numeric: the ratio target used for `tipping_threshold`
#'   (`NA` if not requested).
#' @param evalue Numeric: optional E-value summary (length-2: estimate and CI
#'   bound) from the \pkg{EValue} package; `numeric(0)` if not computed.
#' @param scale Character: label of the ratio scale (e.g. `"Delta"`,
#'   `"incremental"`); free text supplied by the caller.
#' @param call Call: the originating call.
#'
#' @seealso [pmed_sensitivity()]
#' @export
PmedSensitivity <- S7::new_class(
  "PmedSensitivity", package = "probmed",
  properties = list(
    p_med = S7::class_numeric,
    indirect = S7::class_numeric,
    total = S7::class_numeric,
    bias_grid = S7::class_numeric,
    p_med_bias = S7::class_numeric,
    tipping_indirect = S7::class_numeric,
    tipping_threshold = S7::new_property(class = S7::class_numeric, default = NA_real_),
    threshold = S7::new_property(class = S7::class_numeric, default = NA_real_),
    evalue = S7::new_property(class = S7::class_numeric, default = numeric(0)),
    scale = S7::new_property(class = S7::class_character, default = NA_character_),
    call = S7::new_property(class = S7::class_any, default = NULL)
  ),
  validator = function(self) {
    if (length(self@total) && abs(self@total) < 1e-12) {
      return("total (denominator) is ~0: P_med is undefined.")
    }
    if (length(self@bias_grid) != length(self@p_med_bias)) {
      return("bias_grid and p_med_bias must have equal length.")
    }
    NULL
  }
)

#' Extract `(indirect, total)` from a `P_med`-family result.
#'
#' Internal generic so each paper's result class can plug in. The default
#' handles `GaugePmedResult` and any list / object exposing `indirect`+`total`
#' (or `IIE`+`OE`) components. Returns a length-2 numeric `c(indirect, total)`.
#' @noRd
.pmed_components <- function(object) {
  ## GaugePmedResult (Paper 1): numerator = IIE, denominator = OE.
  if (inherits(object, "probmed::GaugePmedResult") ||
      methods::is(object, "probmed::GaugePmedResult")) {
    return(c(indirect = object@IIE, total = object@OE))
  }
  ## Generic S4/S7 with IIE/OE or indirect/total properties.
  get_one <- function(obj, names) {
    for (nm in names) {
      val <- tryCatch(methods::slot(obj, nm), error = function(e) NULL)
      if (is.null(val)) val <- tryCatch(obj[[nm]], error = function(e) NULL)
      if (!is.null(val)) return(as.numeric(val)[1])
    }
    NULL
  }
  ind <- get_one(object, c("indirect", "IIE", "iie", "nie", "NIE"))
  tot <- get_one(object, c("total", "OE", "oe", "te", "TE"))
  if (is.null(ind) || is.null(tot)) {
    stop("Cannot extract indirect/total components. Pass `indirect=`/`total=` ",
         "directly, or supply an object exposing them.", call. = FALSE)
  }
  c(indirect = ind, total = tot)
}

#' M--Y Confounding Sensitivity for a Proportion-Mediated Estimand
#'
#' @description
#' For a mean-based proportion mediated, this function adds a specified bias
#' `b` to the indirect effect while holding the total effect fixed. It reports
#' `(indirect + b) / total` over the bias grid and the bias needed to make the
#' indirect effect zero or bring the ratio to a user-specified threshold.
#'
#' @details
#' The additive bias model applies to the mean-based gauge and incremental
#' shares from [ward_residual()] and [incr_pmed()]. In those frameworks the bias
#' represents a shift in the cross-world mean when their assumption A4 (no
#' exposure-induced mediator--outcome confounding) fails.
#'
#' The tie-adjusted stochastic-superiority probability from [pmed()] and the
#' manuscript's A2/A5 sensitivity analyses are outside the scope of this helper.
#' For a Sobol variance share, callers would need to supply a separately
#' justified perturbation of the variance component.
#'
#' The optional E-value summary (`evalue = TRUE`) is a coarse RR-oriented
#' robustness number; it is reported only as a familiar single figure and is
#' guarded behind the \pkg{EValue} package (Suggests).
#'
#' @param object A [GaugePmedResult] or another object exposing indirect and
#'   total effect components. Ignored if both `indirect` and `total`
#'   are supplied directly.
#' @param indirect,total Numeric: supply the components directly instead of (or
#'   to override) extraction from `object`.
#' @param bias_grid Numeric: numerator biases `b` to evaluate. Default spans
#'   `+/- |indirect|` on an 11-point grid.
#' @param threshold Numeric or `NA`: a ratio value whose tipping bias to
#'   report (e.g. `0` to ask "how much bias makes the share vanish"). Default
#'   `NA` (skip).
#' @param evalue Logical: also compute an \pkg{EValue} E-value summary if the
#'   package is installed (default `FALSE`).
#' @param se Numeric or `NA`: standard error of `indirect`, used only for the
#'   E-value CI bound when `evalue = TRUE`.
#' @param scale Character: label for the `P_med` scale (free text).
#' @param ... Unused.
#'
#' @return A [PmedSensitivity] object.
#'
#' @examples
#' ## From mean-based effect components:
#' pmed_sensitivity(indirect = 0.30, total = 0.50, threshold = 0)
#'
#' \donttest{
#' ## From a Paper 1 fit:
#' set.seed(1)
#' n <- 400; C <- rnorm(n)
#' A <- rbinom(n, 1, plogis(0.5 * C))
#' M <- 0.5 * A + 0.3 * C + rnorm(n)
#' Y <- 0.4 * A + 0.6 * M + 0.2 * C + rnorm(n)
#' fit <- ward_residual(data.frame(A, M, Y, C))
#' pmed_sensitivity(fit, threshold = 0)
#' }
#' @export
pmed_sensitivity <- function(object = NULL, indirect = NULL, total = NULL,
                             bias_grid = NULL, threshold = NA_real_,
                             evalue = FALSE, se = NA_real_,
                             scale = NA_character_, ...) {
  if (is.null(indirect) || is.null(total)) {
    comp <- .pmed_components(object)
    if (is.null(indirect)) indirect <- unname(comp[["indirect"]])
    if (is.null(total))    total    <- unname(comp[["total"]])
  }
  indirect <- as.numeric(indirect)[1]
  total <- as.numeric(total)[1]
  if (!is.finite(total) || abs(total) < 1e-12) {
    stop("total (denominator) is ~0 or non-finite: P_med undefined.", call. = FALSE)
  }

  p_med <- indirect / total

  if (is.null(bias_grid)) {
    span <- max(abs(indirect), .Machine$double.eps)
    bias_grid <- seq(-span, span, length.out = 11L)
  }
  bias_grid <- as.numeric(bias_grid)
  p_med_bias <- (indirect + bias_grid) / total

  ## Tipping bias: indirect + b = 0  =>  b = -indirect.
  tipping_indirect <- -indirect

  ## Tipping bias for a P_med threshold:  (indirect + b)/total = thr
  ##   => b = thr * total - indirect.
  tipping_threshold <- NA_real_
  if (!is.na(threshold)) {
    tipping_threshold <- threshold * total - indirect
  }

  ev <- numeric(0)
  if (isTRUE(evalue)) {
    if (!requireNamespace("EValue", quietly = TRUE)) {
      warning("Package 'EValue' not installed; skipping E-value summary.",
              call. = FALSE)
    } else {
      ## Treat the standardized indirect component as an RR-like estimate.
      ## Coarse, single-number robustness figure only.
      est <- exp(indirect)
      lo <- if (is.finite(se)) exp(indirect - 1.96 * se) else NA_real_
      ev_est <- tryCatch(
        EValue::evalues.RR(est = est,
                           lo = if (is.finite(lo)) lo else NA,
                           hi = NA)["E-values", "point"],
        error = function(e) NA_real_)
      ev_ci <- tryCatch(
        EValue::evalues.RR(est = est,
                           lo = if (is.finite(lo)) lo else NA,
                           hi = NA)["E-values", "lower"],
        error = function(e) NA_real_)
      ev <- c(point = unname(ev_est), lower = unname(ev_ci))
    }
  }

  PmedSensitivity(
    p_med = p_med, indirect = indirect, total = total,
    bias_grid = bias_grid, p_med_bias = p_med_bias,
    tipping_indirect = tipping_indirect,
    tipping_threshold = tipping_threshold,
    threshold = as.numeric(threshold),
    evalue = ev, scale = as.character(scale),
    call = match.call()
  )
}

#' @export
S7::method(print, PmedSensitivity) <- function(x, ...) {
  lbl <- if (!is.na(x@scale)) paste0(" [", x@scale, "]") else ""
  cat("P_med M-Y confounding sensitivity (A4)", lbl, "\n", sep = "")
  cat(sprintf("  P_med = %.3f  (indirect=%.3f / total=%.3f)\n",
              x@p_med, x@indirect, x@total))
  cat(sprintf("  Tipping bias (indirect -> 0): b = %.3f\n", x@tipping_indirect))
  if (!is.na(x@tipping_threshold)) {
    cat(sprintf("  Tipping bias (P_med -> %.3f): b = %.3f\n",
                x@threshold, x@tipping_threshold))
  }
  rng <- range(x@p_med_bias)
  cat(sprintf("  P_med over bias grid [%.3f, %.3f]: [%.3f, %.3f]\n",
              min(x@bias_grid), max(x@bias_grid), rng[1], rng[2]))
  if (length(x@evalue)) {
    cat(sprintf("  E-value (point) = %.3f\n", x@evalue[["point"]]))
  }
  invisible(x)
}
