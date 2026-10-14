#' Observed-sample covariate standardization
#'
#' For a binary mediator or outcome, covariates do not cancel in the
#' cross-world contrast, so P_med is standardized over the empirical covariate
#' distribution of the analysis rows: each of the `n_sim` pairs draws one row
#' uniformly and uses its covariate values in both arms (Definition 1 with the
#' target held fixed at the observed rows). Gaussian additive models are not
#' standardized: their covariate shifts cancel exactly.
#'
#' medfit stores predictors as model-matrix column names (`"C"`, `"Fb"`,
#' `"log(abs(C) + 1)"`), so the design is rebuilt column by column from the
#' stored data.
#'
#' @name standardize
#' @keywords internal
NULL

#' Data variables a model-matrix column label refers to
#' @keywords internal
.label_vars <- function(label, data) {
  direct <- tryCatch(all.vars(str2lang(label)), error = function(e) character(0))
  direct <- intersect(direct, names(data))
  if (length(direct)) return(direct)
  # Factor dummy column: <variable><level>
  is_cat <- vapply(data, function(v) is.factor(v) || is.character(v), logical(1))
  cands <- names(data)[is_cat]
  cands[vapply(cands, function(v) startsWith(label, v), logical(1))]
}

#' Rebuild one model-matrix column from the stored data
#' @keywords internal
.design_column <- function(label, data) {
  n <- nrow(data)
  if (label %in% names(data) && (is.numeric(data[[label]]) || is.logical(data[[label]]))) {
    return(as.numeric(data[[label]]))
  }
  # `label` is a term of the user's own fitted formula (e.g. `log(abs(C) + 1)`),
  # which glm() already evaluated on this data; re-evaluating it here crosses no
  # new trust boundary.
  val <- tryCatch(eval(str2lang(label), data, asNamespace("stats")), error = function(e) NULL)
  if (is.numeric(val) && length(val) == n) return(as.numeric(val))
  for (v in .label_vars(label, data)) {
    x <- data[[v]]
    if (is.factor(x) || is.character(x)) {
      level <- substring(label, nchar(v) + 1L)
      if (level %in% unique(as.character(x))) return(as.numeric(as.character(x) == level))
    }
  }
  stop("pmed() cannot rebuild the model column `", label, "` from the stored ",
       "data. Use the formula interface or simpler covariate terms.",
       call. = FALSE)
}

#' Rebuild the mediator and outcome design matrices of an extracted model
#'
#' Column names match the coefficient names without their `m_` / `y_` prefix.
#' @keywords internal
.design_matrices <- function(extract) {
  data <- extract@data
  build <- function(prefix, predictors) {
    labels <- setdiff(predictors, "(Intercept)")
    cols <- lapply(labels, .design_column, data = data)
    D <- do.call(cbind, c(list(), cols))
    if (is.null(D)) D <- matrix(numeric(0), nrow(data), 0)
    colnames(D) <- labels
    if (paste0(prefix, "(Intercept)") %in% names(extract@estimates)) {
      D <- cbind("(Intercept)" = 1, D)
    }
    D
  }
  list(m = build("m_", extract@mediator_predictors),
       y = build("y_", extract@outcome_predictors))
}

#' Does this model need covariate standardization?
#' @keywords internal
.needs_standardization <- function(extract) {
  length(.covariate_terms(extract)) > 0 &&
    !(.is_gaussian(extract@family_m) && .is_gaussian(extract@family_y))
}

#' Per-row covariate linear predictors at x = 0 and m = 0
#'
#' `coef_m` / `coef_y` are named by design column (no prefix). Returns the
#' mediator and outcome linear predictors with the treatment and mediator
#' contributions removed, one value per analysis row.
#' @keywords internal
.std_linear_predictors <- function(design, coef_m, coef_y, treatment, mediator) {
  keep_m <- setdiff(colnames(design$m), treatment)
  keep_y <- setdiff(colnames(design$y), c(treatment, mediator))
  list(
    m = drop(design$m[, keep_m, drop = FALSE] %*% coef_m[keep_m]),
    y = drop(design$y[, keep_y, drop = FALSE] %*% coef_y[keep_y])
  )
}

#' Split a joint coefficient vector into mediator and outcome parts
#' @keywords internal
.split_coefs <- function(theta, design) {
  list(
    m = stats::setNames(unname(theta[paste0("m_", colnames(design$m))]), colnames(design$m)),
    y = stats::setNames(unname(theta[paste0("y_", colnames(design$y))]), colnames(design$y))
  )
}

#' Standardization inputs for one coefficient vector, or NULLs if not needed
#'
#' `theta` is a joint coefficient vector named `m_<col>` / `y_<col>`.
#' @keywords internal
.std_for <- function(extract, theta, design = NULL) {
  if (!.needs_standardization(extract)) return(list(m = NULL, y = NULL))
  if (is.null(design)) design <- .design_matrices(extract)
  co <- .split_coefs(theta, design)
  .std_linear_predictors(design, co$m, co$y, extract@treatment, extract@mediator)
}
