# ============================================================
# LongOmicsTraj generics
# ============================================================
#' @import data.table
#' @importFrom data.table .I .N .EACHI .BY
NULL

#' Create a LongOmicsTraj object
#' @export
setGeneric("lot_create", function(...) {
  standardGeneric("lot_create")
})

#' Access the MultiAssayExperiment
#' @export
setGeneric("lot_mae", function(object) {
  standardGeneric("lot_mae")
})

#' Access the Visit Mean Matrix
#' @export
setGeneric("lot_vmm", function(object) {
  standardGeneric("lot_vmm")
})

#' Access adjacent deltas
#' @export
setGeneric("lot_deltas", function(object) {
  standardGeneric("lot_deltas")
})

#' Access thresholds
#' @export
setGeneric("lot_thresholds", function(object) {
  standardGeneric("lot_thresholds")
})

#' Access Overall Trajectory States
#' @export
setGeneric("lot_ots", function(object) {
  standardGeneric("lot_ots")
})

#' Access object relationships
#' @export
setGeneric("lot_relationships", function(object) {
  standardGeneric("lot_relationships")
})

#' Access stored results
#' @export
setGeneric("lot_results", function(object) {
  standardGeneric("lot_results")
})

#' Map trajectories to OTS labels
#' @export
setGeneric(
  "lot_map_to_ots",
  function(object, visits, overwrite = FALSE, ...)
    standardGeneric("lot_map_to_ots")
)
#' Query topology labels
#' @export
setGeneric("lot_query", function(object, ...) {
  standardGeneric("lot_query")
})

#' Compare topology assignments
#' @export
setGeneric("lot_compare", function(object, ...) {
  standardGeneric("lot_compare")
})

#' Visualise topology outputs
#' @export
setGeneric("lot_visualise", function(object, ...) {
  standardGeneric("lot_visualise")
})
#' Compute topology thresholds
#' @export
setGeneric("lot_compute_thresholds", function(object, ...) {
  standardGeneric("lot_compute_thresholds")
})

#' Query OTS topology labels
#' @export
setGeneric("lot_query", function(object, ...) {
  standardGeneric("lot_query")
})

#' Build VMM from linear mixed-effects models
#' @export
setGeneric("lot_from_lmm", function(object, ...) {
  standardGeneric("lot_from_lmm")
})

setGeneric("lot_from_gam", function(object, ...) {
  standardGeneric("lot_from_gam")
})
setGeneric("lot_from_masigpro", function(object, ...) {
  standardGeneric("lot_from_masigpro")
})
setGeneric("lot_compare_estimators", function(...) {
  standardGeneric("lot_compare_estimators")
})
setGeneric("lot_consensus", function(x, min_agreement = 2, mode = "majority") {
  standardGeneric("lot_consensus")
})
setGeneric("lot_from_clusters", function(object, ...) {
  standardGeneric("lot_from_clusters")
})
setGeneric("lot_compute_deltas", function(object, ...) {
  standardGeneric("lot_compute_deltas")
})
