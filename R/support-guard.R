#' Supported-model guard for single-mediator P_med
#'
#' The core draws Gaussian responses (identity link) or Bernoulli responses
#' (logit or probit link) and nothing else. Fitting a GLM is not the same as
#' being able to generate its response law, so every other family, link and
#' non-additive structure is rejected here, before any random draw.
#'
#' @name support-guard
#' @keywords internal
NULL

#' Is a (family, link) pair one the core can simulate?
#' @keywords internal
.supported_family <- function(family) {
  if (is.null(family)) return(TRUE)
  fam <- family$family
  link <- family$link
  (identical(fam, "gaussian") && identical(link, "identity")) ||
    (identical(fam, "binomial") && link %in% c("logit", "probit"))
}

#' Describe a family for an error message
#' @keywords internal
.family_label <- function(family) {
  if (is.null(family)) return("gaussian(identity)")
  paste0(family$family, "(", family$link, ")")
}

#' Model terms that involve a variable other than as its own main effect
#'
#' Returns the term labels (e.g. `"X:M"`, `"I(M^2)"`, `"log(X)"`) that mention
#' any of `vars` without being exactly one of them.
#' @keywords internal
.offending_terms <- function(labels, vars) {
  labels <- setdiff(labels, "(Intercept)")
  bad <- vapply(labels, function(lab) {
    if (lab %in% vars) return(FALSE)
    used <- tryCatch(all.vars(str2lang(lab)), error = function(e) lab)
    any(used %in% vars)
  }, logical(1))
  labels[bad]
}

#' Reject structures a formula would fit but the core cannot represent
#'
#' Called by the formula method before any model is fitted.
#' @keywords internal
.check_supported_formulas <- function(formula_y, formula_m, treatment, mediator) {
  vars <- c(treatment, mediator)
  bad_y <- .offending_terms(attr(stats::terms(formula_y), "term.labels"), vars)
  bad_m <- .offending_terms(attr(stats::terms(formula_m), "term.labels"), vars)
  if (!is.null(attr(stats::terms(formula_y), "offset")) ||
      !is.null(attr(stats::terms(formula_m), "offset"))) {
    stop("pmed() does not support `offset()` terms: the simulation core ",
         "would ignore them.", call. = FALSE)
  }
  bad <- c(bad_m, bad_y)
  if (length(bad)) {
    stop(
      "pmed() supports additive single-mediator models only: the treatment `",
      treatment, "` and mediator `", mediator, "` must enter as main effects. ",
      "Unsupported term(s): ", paste0("`", unique(bad), "`", collapse = ", "),
      ". Interactions and transformed treatment/mediator terms are not ",
      "implemented.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate an extracted single-mediator model before computation
#'
#' Checks the response families and links, the term structure, and that every
#' coefficient was estimated (rank deficiency leaves aliased coefficients `NA`).
#' @keywords internal
.check_supported_extract <- function(extract) {
  for (side in c("mediator", "outcome")) {
    family <- if (side == "mediator") extract@family_m else extract@family_y
    if (!.supported_family(family)) {
      stop(
        "pmed() cannot simulate a ", .family_label(family), " ", side,
        ". Supported: gaussian(identity), binomial(logit), binomial(probit). ",
        "Count, ordinal and other families are not implemented.",
        call. = FALSE
      )
    }
  }
  if (.needs_standardization(extract) && is.null(extract@data)) {
    # Covariates do not cancel under a nonlinear link; standardizing over
    # them needs the analysis rows.
    stop("pmed() needs the source data to standardize over covariates (",
         paste0("`", .covariate_terms(extract), "`", collapse = ", "),
         ") with a binary mediator or outcome.", call. = FALSE)
  }
  vars <- c(extract@treatment, extract@mediator)
  bad <- c(
    .offending_terms(extract@mediator_predictors, vars),
    .offending_terms(extract@outcome_predictors, vars)
  )
  if (length(bad)) {
    stop(
      "pmed() supports additive single-mediator models only. Unsupported ",
      "term(s): ", paste0("`", unique(bad), "`", collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!extract@treatment %in% extract@mediator_predictors ||
      !extract@mediator %in% extract@outcome_predictors) {
    stop(
      "pmed() needs the treatment as a predictor of the mediator and the ",
      "mediator as a predictor of the outcome.",
      call. = FALSE
    )
  }
  .check_supported_data(extract)
  est <- extract@estimates
  model_coefs <- est[grepl("^(m|y)_", names(est))]
  if (anyNA(model_coefs)) {
    stop(
      "pmed() found unestimated (aliased) coefficients: ",
      paste0("`", names(model_coefs)[is.na(model_coefs)], "`", collapse = ", "),
      ". The design is rank deficient.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate the numerical controls shared by every method
#' @keywords internal
.check_controls <- function(x_ref, x_value, n_boot, ci_level, seed, n_sim,
                            method) {
  scalar_num <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  count <- function(x) scalar_num(x) && x >= 1 && x == round(x)
  if (!scalar_num(x_ref) || !scalar_num(x_value)) {
    stop("`x_ref` and `x_value` must be single finite numbers.", call. = FALSE)
  }
  if (method %in% c("parametric_bootstrap", "nonparametric_bootstrap",
                    "model_bootstrap") &&
      !count(n_boot)) {
    stop("`n_boot` must be a single positive whole number.", call. = FALSE)
  }
  if (!count(n_sim)) {
    stop("`n_sim` must be a single positive whole number.", call. = FALSE)
  }
  if (!scalar_num(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!is.null(seed) && !scalar_num(seed)) {
    stop("`seed` must be NULL or a single finite number.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Covariates the model uses beyond the treatment and mediator
#' @keywords internal
.covariate_terms <- function(extract) {
  vars <- c(extract@treatment, extract@mediator, "(Intercept)")
  unique(c(setdiff(extract@mediator_predictors, vars),
           setdiff(extract@outcome_predictors, vars)))
}

#' Data-level checks: missing rows, variable types, binary responses
#'
#' The nonparametric bootstrap resamples `extract@data`, so rows a model fit
#' dropped for missingness would come back; and a binomial response that is not
#' coded 0/1 (proportions, successes/failures) cannot be simulated as Bernoulli.
#' @keywords internal
.check_supported_data <- function(extract) {
  data <- extract@data
  if (is.null(data)) return(invisible(TRUE))
  labels <- c(extract@mediator_predictors, extract@outcome_predictors)
  labels <- setdiff(labels, "(Intercept)")
  used <- unique(c(extract@treatment, extract@mediator, extract@outcome,
                   unlist(lapply(labels, .label_vars, data = data))))
  used <- intersect(used, names(data))
  if (anyNA(data[, used, drop = FALSE])) {
    stop("pmed() found missing values in the model variables (",
         paste0("`", used[vapply(data[used], anyNA, logical(1))], "`",
                collapse = ", "),
         "). Use complete cases before fitting; the bootstrap would ",
         "otherwise resample rows the fit dropped.", call. = FALSE)
  }
  for (v in c(extract@treatment, extract@mediator)) {
    if (!is.numeric(data[[v]]) && !is.logical(data[[v]])) {
      stop("pmed() needs a numeric `", v, "` (code a binary variable as 0/1).",
           call. = FALSE)
    }
  }
  sides <- list(list(fam = extract@family_m, var = extract@mediator),
                list(fam = extract@family_y, var = extract@outcome))
  for (sd in sides) {
    if (!.is_gaussian(sd$fam)) {
      y <- data[[sd$var]]
      if (!is.null(y) && !all(y %in% c(0, 1))) {
        stop("pmed() simulates a binomial `", sd$var, "` as Bernoulli and ",
             "needs it coded 0/1; proportions and successes/failures are not ",
             "supported.", call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}
