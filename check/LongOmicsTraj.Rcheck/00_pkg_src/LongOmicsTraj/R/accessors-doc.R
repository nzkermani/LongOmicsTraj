#' Access components of a LongOmicsTraj object
#'
#' Accessor functions for retrieving the principal data structures stored in a
#' `LongOmicsTraj` object. These functions provide a stable interface to the
#' object and should be used instead of direct slot access.
#'
#' \describe{
#'   \item{`lot_mae()`}{
#'     Returns the underlying
#'     \code{\link[MultiAssayExperiment]{MultiAssayExperiment}} containing the
#'     molecular assays and sample metadata.
#'   }
#'   \item{`lot_vmm()`}{
#'     Returns the visit-level molecular model matrix containing estimated
#'     molecular values for each feature, group and visit.
#'   }
#'   \item{`lot_deltas()`}{
#'     Returns estimated molecular changes between consecutive visits.
#'   }
#'   \item{`lot_thresholds()`}{
#'     Returns the feature-specific thresholds used to classify molecular
#'     changes as increasing, decreasing or stable.
#'   }
#'   \item{`lot_ots()`}{
#'     Returns the Omics Trajectory Signature assignments obtained by combining
#'     the classified transitions across visits.
#'   }
#'   \item{`lot_relationships()`}{
#'     Returns relationships among molecular objects, such as
#'     feature-to-cluster assignments.
#'   }
#'   \item{`lot_results()`}{
#'     Returns additional analysis results stored by trajectory estimators,
#'     clustering procedures and diagnostic functions.
#'   }
#' }
#'
#' @param object A `LongOmicsTraj` object.
#'
#' @return
#' The requested component of the `LongOmicsTraj` object. The exact class
#' depends on the accessor:
#'
#' \itemize{
#'   \item `lot_mae()` returns a `MultiAssayExperiment`.
#'   \item `lot_vmm()`, `lot_deltas()`, `lot_thresholds()`,
#'     `lot_ots()` and `lot_relationships()` return
#'     `S4Vectors::DataFrame` objects.
#'   \item `lot_results()` returns a named list.
#' }
#'
#' @name lot_accessors
#'
#' @aliases lot_mae
#' @aliases lot_vmm
#' @aliases lot_deltas
#' @aliases lot_thresholds
#' @aliases lot_ots
#' @aliases lot_relationships
#' @aliases lot_results
#' @aliases lot_mae,LongOmicsTraj-method
#' @aliases lot_vmm,LongOmicsTraj-method
#' @aliases lot_deltas,LongOmicsTraj-method
#' @aliases lot_thresholds,LongOmicsTraj-method
#' @aliases lot_ots,LongOmicsTraj-method
#' @aliases lot_relationships,LongOmicsTraj-method
#' @aliases lot_results,LongOmicsTraj-method
#'
#' @examples
#' data("lot_glucold_ots")
#'
#' ## Extract the underlying experiment container.
#' mae <- lot_mae(lot_glucold_ots)
#' names(MultiAssayExperiment::experiments(mae))
#'
#' ## Inspect the visit-level molecular estimates.
#' vmm <- lot_vmm(lot_glucold_ots)
#' dim(vmm)
#' utils::head(as.data.frame(vmm))
#'
#' ## Examine the estimated changes between adjacent visits.
#' deltas <- lot_deltas(lot_glucold_ots)
#' dim(deltas)
#' utils::head(
#'   as.data.frame(deltas)[
#'     ,
#'     c(
#'       "object_id",
#'       "group",
#'       "transition",
#'       "delta_hat"
#'     )
#'   ]
#' )
#'
#' ## Inspect the thresholds used to classify molecular changes.
#' thresholds <- lot_thresholds(lot_glucold_ots)
#' utils::head(
#'   as.data.frame(thresholds)[
#'     ,
#'     c(
#'       "object_id",
#'       "threshold_method",
#'       "delta_g"
#'     )
#'   ]
#' )
#'
#' ## Inspect the resulting Omics Trajectory Signatures.
#' ots <- lot_ots(lot_glucold_ots)
#' table(
#'   as.data.frame(ots)$topology_label
#' )
#'
#' ## Relationships and additional results may be empty when no
#' ## clustering analysis has been performed.
#' relationships <- lot_relationships(lot_glucold_ots)
#' nrow(relationships)
#'
#' results <- lot_results(lot_glucold_ots)
#' names(results)
NULL