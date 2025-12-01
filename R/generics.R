#' Extract Mediation Structure
#'
#' @description
#' Generic function to extract mediation structure from fitted models.
#' See method documentation for specific arguments.
#'
#' @param object Fitted model object
#' @param ... Additional arguments passed to methods
#'
#' @return MediationExtract object
#' @export
extract_mediation <- S7::new_generic(
  "extract_mediation",
  dispatch_args = "object"
)

#' Compute P_med
#'
#' @description
#' Generic function to compute P_med from various inputs.
#'
#' @param object Either a formula, MediationExtract, or fitted model
#' @param ... Additional arguments passed to methods
#'
#' @return PmedResult object
#' @export
pmed <- S7::new_generic("pmed", dispatch_args = "object")
