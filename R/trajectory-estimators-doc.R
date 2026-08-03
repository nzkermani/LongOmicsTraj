#' Estimate longitudinal molecular trajectories
#'
#' Estimate feature-level longitudinal molecular trajectories from repeated
#' measurements stored in a `LongOmicsTraj` object.
#'
#' The available methods include empirical visit-level summaries, linear
#' mixed-effects models, generalised additive models and maSigPro-style
#' polynomial regression.
#'
#' @param object A `LongOmicsTraj` object.
#' @param ... Additional method-specific arguments passed to the selected
#'   trajectory-estimation method.
#' @param assay Name of the assay or experiment to analyse.
#' @param visits Character vector giving the chronological visit order.
#' @param group_col Metadata column containing clinical or treatment groups.
#' @param visit_col Metadata column containing visit labels.
#' @param estimator Character label recorded for the fitted estimator.
#' @param overwrite Logical; overwrite existing visit-level trajectory
#'   results.
#' @param verbose Logical; print progress messages.
#' @param formula_fixed Fixed-effects model formula used by the linear
#'   mixed-effects estimator.
#' @param subject_col Metadata column containing subject identifiers.
#' @param BPPARAM A `BiocParallelParam` object controlling parallel execution.
#' @param time_col Metadata column containing numeric time values.
#' @param degree Polynomial degree used by the regression-based trajectory
#'   estimator.
#'
#' @return A `LongOmicsTraj` object with updated visit-level molecular model
#' results stored in the `vmm` component.
#'
#' @name lot_trajectory_estimators
#'
#' @aliases lot_from_empirical
#' @aliases lot_from_lmm
#' @aliases lot_from_gam
#' @aliases lot_from_masigpro
#' @aliases lot_from_empirical,LongOmicsTraj-method
#' @aliases lot_from_lmm,LongOmicsTraj-method
#' @aliases lot_from_gam,LongOmicsTraj-method
#' @aliases lot_from_masigpro,LongOmicsTraj-method
#'
#' @examples
#' data("lot_glucold")
#'
#' visits <- c(
#'   "Baseline",
#'   "6 Months",
#'   "30 Months"
#' )
#'
#' fitted <- lot_from_empirical(
#'   object = lot_glucold,
#'   assay = "transcriptomics",
#'   visits = visits,
#'   overwrite = TRUE
#' )
#'
#' utils::head(
#'   as.data.frame(
#'     lot_vmm(fitted)
#'   )
#' )
NULL