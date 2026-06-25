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
