#' M--Y Confounding Sensitivity for the Incremental Mediated Elasticity
#'
#' @description
#' Per-`delta` M--Y unmeasured-confounding sensitivity for an [IncrPmedResult].
#' For each tilt factor `delta`, reports the additive numerator bias that zeroes
#' the incremental mediated share (`tipping`) and, optionally, the bias that
#' drives `P_med^delta` to a user `threshold`. Reuses [pmed_sensitivity()];
#' the additive-offset model is exact for the Paper 2 (mean-based) share.
#'
#' @param object An [IncrPmedResult] from [incr_pmed()].
#' @param threshold Numeric or `NA`: a `P_med` value whose tipping bias to
#'   report per `delta` (default `0` = "bias that makes the share vanish").
#' @param ... Unused.
#'
#' @return A data frame: one row per `delta` with columns `delta`, `med`,
#'   `tot`, `Pmed`, `tipping`, `tipping_threshold`.
#'
#' @examples
#' set.seed(1)
#' n <- 800; C <- rnorm(n)
#' A <- rbinom(n, 1, plogis(-0.2 + 0.8 * C))
#' M <- 0.6 * A + 0.4 * C + rnorm(n)
#' Y <- 0.5 * A + 0.7 * M + 0.3 * C + rnorm(n)
#' fit <- incr_pmed(data.frame(A, M, Y, C), deltas = c(0.5, 1, 2))
#' incr_sensitivity(fit, threshold = 0)
#'
#' @export
incr_sensitivity <- function(object, threshold = 0, ...) {
  if (!inherits(object, "probmed::IncrPmedResult") &&
      !methods::is(object, "probmed::IncrPmedResult")) {
    stop("`object` must be an IncrPmedResult (from incr_pmed()).", call. = FALSE)
  }
  .incr_sensitivity_from_curve(object@curve, threshold = threshold)
}

#' Internal: build the per-delta sensitivity table from an IncrPmedResult curve.
#' @noRd
.incr_sensitivity_from_curve <- function(curve, threshold = 0) {
  rows <- lapply(seq_len(nrow(curve)), function(i) {
    s <- probmed::pmed_sensitivity(
      indirect = curve$med[i], total = curve$tot[i], threshold = threshold
    )
    data.frame(
      delta = curve$delta[i], med = curve$med[i], tot = curve$tot[i],
      Pmed = s@p_med, tipping = s@tipping_indirect,
      tipping_threshold = s@tipping_threshold
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
