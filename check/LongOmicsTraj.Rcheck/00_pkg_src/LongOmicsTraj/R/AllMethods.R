# ============================================================
# LongOmicsTraj accessor methods
# ============================================================

#' @export
#' @rdname lot_accessors
#' @export
setMethod("lot_mae", "LongOmicsTraj", function(object) {
  object@mae
})

#' @export
#' @rdname lot_accessors
#' @export
setMethod("lot_vmm", "LongOmicsTraj", function(object) {
  object@vmm
})

#' @export
#' @rdname lot_accessors
#' @export
setMethod("lot_deltas", "LongOmicsTraj", function(object) {
  object@deltas
})

#' @export
#' @rdname lot_accessors
#' @export
setMethod("lot_thresholds", "LongOmicsTraj", function(object) {
  object@thresholds
})

#' @export
#' @rdname lot_accessors
#' @export
setMethod("lot_ots", "LongOmicsTraj", function(object) {
  object@ots
})

#' @export
#' @rdname lot_accessors
#' @export
setMethod("lot_relationships", "LongOmicsTraj", function(object) {
  object@relationships
})

#' @export
#' @rdname lot_accessors
#' @export
setMethod("lot_results", "LongOmicsTraj", function(object) {
  object@results
})


