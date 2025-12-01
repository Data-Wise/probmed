#' Print PmedResult
#'
#' @name print-PmedResult
S7::method(print, PmedResult) <- function(x, digits = 3, ...) {
  cat("\n")
  cat("P_med: Probability of Mediated Shift\n")
  cat("====================================\n\n")
  
  cat("Estimate:", format(x@estimate, digits = digits), "\n")
  
  if (!is.na(x@ci_lower)) {
    cat(paste0(x@ci_level * 100, "% CI: ["),
        format(x@ci_lower, digits = digits), ", ",
        format(x@ci_upper, digits = digits), "]\n")
  }
  
  cat("\n")
  cat("Inference:", x@method, "\n")
  if (!is.na(x@n_boot)) {
    cat("Bootstrap samples:", x@n_boot, "\n")
  }
  
  cat("\nTreatment contrast: X =", x@x_value, "vs. X* =", x@x_ref, "\n")
  cat("\nInterpretation:\n")
  cat("  P(Y_{X*, M_X} > Y_{X, M_X}) =", format(x@estimate, digits = digits), "\n")
  cat("\n")
  
  invisible(x)
}

#' Summary Method
#'
#' @name summary-PmedResult
S7::method(summary, PmedResult) <- function(object, ...) {
  print(object, ...)
  
  if (length(object@boot_estimates) > 0) {
    cat("Bootstrap Distribution:\n")
    print(summary(object@boot_estimates))
    cat("\nStandard Error:", stats::sd(object@boot_estimates, na.rm = TRUE), "\n")
  }
  
  cat("\nSource:", object@source_extract@source_package, "\n")
  cat("Sample size:", object@source_extract@n_obs, "\n")
  
  invisible(object)
}

#' Plot Bootstrap Distribution
#'
#' @name plot-PmedResult
S7::method(plot, PmedResult) <- function(x, ...) {
  if (length(x@boot_estimates) == 0) {
    stop("No bootstrap distribution to plot. Use method = 'parametric_bootstrap' or 'nonparametric_bootstrap'")
  }
  
  graphics::hist(
    x@boot_estimates,
    main = "Bootstrap Distribution of P_med",
    xlab = "P_med",
    col = "skyblue",
    border = "white",
    xlim = c(0, 1),
    ...
  )
  
  graphics::abline(v = x@estimate, col = "red", lwd = 2, lty = 2)
  graphics::abline(v = c(x@ci_lower, x@ci_upper), col = "blue", lwd = 2, lty = 2)
  
  graphics::legend(
    "topright",
    legend = c("Point Estimate", paste0(x@ci_level * 100, "% CI")),
    col = c("red", "blue"),
    lty = 2,
    lwd = 2
  )
  
  invisible(x)
}