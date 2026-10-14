#' Deterministic separation check for a binary-response design
#'
#' Complete or quasi-complete separation exists exactly when some nonzero
#' `beta` satisfies `s_i * x_i' beta >= 0` for every row, with
#' `s_i = 2 y_i - 1`, and strict inequality for at least one row (Albert and
#' Anderson, 1984). This solves the linear program
#'
#'   maximize   sum_i s_i x_i' beta
#'   subject to s_i x_i' beta >= 0 for all i,  -1 <= beta_j <= 1,
#'
#' (Konis, 2007) with `lpSolve`, splitting `beta = u - v`, `u, v >= 0`. A
#' positive optimum means separation; zero means none. Fitted probabilities
#' near 0/1 or glm warnings are not used: they are not proof either way.
#'
#' @param X Numeric design matrix (including any intercept column).
#' @param y 0/1 response.
#' @param tol Optimum above `tol * max(1, sum(abs(X)))` counts as separation.
#' @return `"separated"`, `"ok"`, or `"unchecked"` when `lpSolve` is not
#'   installed or the LP fails.
#' @keywords internal
.detect_separation <- function(X, y, tol = 1e-9) {
  if (!requireNamespace("lpSolve", quietly = TRUE)) return("unchecked")
  X <- as.matrix(X)
  if (length(unique(y)) < 2L) return("separated")
  s <- 2 * as.numeric(y) - 1
  A <- X * s
  p <- ncol(X)
  obj <- colSums(A)
  con <- rbind(cbind(A, -A), diag(2 * p))
  dir <- c(rep(">=", nrow(A)), rep("<=", 2 * p))
  rhs <- c(rep(0, nrow(A)), rep(1, 2 * p))
  sol <- tryCatch(
    lpSolve::lp("max", c(obj, -obj), con, dir, rhs),
    error = function(e) NULL
  )
  if (is.null(sol) || sol$status != 0) return("unchecked")
  if (sol$objval > tol * max(1, sum(abs(X)))) "separated" else "ok"
}

#' Separation status of an extracted model's binary equations
#'
#' @return Named character vector over the binary sides only
#'   (`mediator`, `outcome`); empty when both equations are Gaussian or no
#'   data are stored.
#' @keywords internal
.extract_separation <- function(extract, design = NULL, rows = NULL) {
  data <- extract@data
  out <- character(0)
  if (is.null(data)) return(out)
  binary_m <- !.is_gaussian(extract@family_m)
  binary_y <- !.is_gaussian(extract@family_y)
  if (!binary_m && !binary_y) return(out)
  if (is.null(design)) design <- .design_matrices(extract)
  if (is.null(rows)) rows <- seq_len(nrow(data))
  if (binary_m) {
    out[["mediator"]] <- .detect_separation(design$m[rows, , drop = FALSE],
                                            data[[extract@mediator]][rows])
  }
  if (binary_y) {
    out[["outcome"]] <- .detect_separation(design$y[rows, , drop = FALSE],
                                           data[[extract@outcome]][rows])
  }
  out
}

#' Stop if the original fit of a binary equation is separated
#' @keywords internal
.check_original_separation <- function(extract) {
  status <- .extract_separation(extract)
  sep <- names(status)[status == "separated"]
  if (length(sep)) {
    stop("pmed() detected complete or quasi-complete separation in the ",
         paste(sep, collapse = " and "), " model; its maximum-likelihood ",
         "estimates do not exist, so P_med cannot be estimated from this fit.",
         call. = FALSE)
  }
  invisible(status)
}
