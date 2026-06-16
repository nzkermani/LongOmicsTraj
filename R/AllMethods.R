# ============================================================
# LongOmicsTraj accessor methods
# ============================================================

#' @export
setMethod("lot_mae", "LongOmicsTraj", function(object) {
  object@mae
})

#' @export
setMethod("lot_vmm", "LongOmicsTraj", function(object) {
  object@vmm
})

#' @export
setMethod("lot_deltas", "LongOmicsTraj", function(object) {
  object@deltas
})

#' @export
setMethod("lot_thresholds", "LongOmicsTraj", function(object) {
  object@thresholds
})

#' @export
setMethod("lot_ots", "LongOmicsTraj", function(object) {
  object@ots
})

#' @export
setMethod("lot_relationships", "LongOmicsTraj", function(object) {
  object@relationships
})

#' @export
setMethod("lot_results", "LongOmicsTraj", function(object) {
  object@results
})

#' Build VMM from empirical visit means
#' @export
setGeneric("lot_from_empirical", function(object, ...) {
  standardGeneric("lot_from_empirical")
})
