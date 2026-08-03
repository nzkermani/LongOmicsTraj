#' Internal OTS computation
#'
#' Converts molecular deltas into transition states using feature-specific
#' thresholds and combines those states into Omics Trajectory Signature
#' labels.
#'
#' @param deltas A data frame or `DataFrame` containing molecular-change
#'   estimates.
#' @param thresholds A data frame or `DataFrame` containing feature-specific
#'   thresholds.
#'
#' @return A `DataFrame` containing Omics Trajectory Signature assignments.
#'
#' @keywords internal
.compute_ots <- function(
    deltas,
    thresholds
) {
  
  required_deltas <- c(
    "object_id",
    "object_type",
    "assay",
    "group",
    "estimator",
    "transition_index",
    "from_visit",
    "to_visit",
    "transition",
    "delta_hat"
  )
  
  required_thresholds <- c(
    "object_id",
    "assay",
    "delta_g"
  )
  
  missing_deltas <- setdiff(
    required_deltas,
    colnames(deltas)
  )
  
  missing_thresholds <- setdiff(
    required_thresholds,
    colnames(thresholds)
  )
  
  if (length(missing_deltas) > 0L) {
    stop(
      "`deltas` is missing required columns: ",
      paste(
        missing_deltas,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  if (length(missing_thresholds) > 0L) {
    stop(
      "`thresholds` is missing required columns: ",
      paste(
        missing_thresholds,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  d <- data.table::as.data.table(
    as.data.frame(deltas)
  )
  
  t <- data.table::as.data.table(
    as.data.frame(thresholds)
  )
  
  t <- t[
    ,
    .(
      object_id,
      assay,
      delta_g
    )
  ]
  
  data.table::setkeyv(
    d,
    c(
      "object_id",
      "assay"
    )
  )
  
  data.table::setkeyv(
    t,
    c(
      "object_id",
      "assay"
    )
  )
  
  d <- t[d]
  
  if (anyNA(d$delta_g)) {
    
    missing_ids <- unique(
      d$object_id[
        is.na(d$delta_g)
      ]
    )
    
    stop(
      "Missing thresholds for objects: ",
      paste(
        utils::head(
          missing_ids,
          10L
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  d[
    ,
    state := "flat"
  ]
  
  d[
    delta_hat > delta_g,
    state := "up"
  ]
  
  d[
    delta_hat < -delta_g,
    state := "down"
  ]
  
  data.table::setorderv(
    d,
    c(
      "object_id",
      "object_type",
      "assay",
      "group",
      "estimator",
      "transition_index"
    )
  )
  
  ots_dt <- d[
    ,
    .(
      topology_label = paste(
        state,
        collapse = "_"
      ),
      n_transitions = .N
    ),
    by = .(
      object_id,
      object_type,
      assay,
      group,
      estimator
    )
  ]
  
  S4Vectors::DataFrame(
    as.data.frame(ots_dt)
  )
}


#' Map molecular changes to Omics Trajectory Signatures
#'
#' Classifies feature-level molecular changes between consecutive visits and
#' combines the transition states into Omics Trajectory Signature labels.
#'
#' The method uses delta estimates and thresholds already stored in the
#' `LongOmicsTraj` object.
#'
#' @param object A `LongOmicsTraj` object.
#' @param visits Character vector giving the visits in chronological order.
#' @param overwrite Logical; overwrite existing OTS results when `TRUE`.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @return A `LongOmicsTraj` object with updated OTS assignments.
#' @examples
#' data("lot_glucold_ots")
#'
#' visits <- c(
#'   "Baseline",
#'   "6 Months",
#'   "30 Months"
#' )
#'
#' fitted <- lot_map_to_ots(
#'   object = lot_glucold_ots,
#'   visits = visits,
#'   overwrite = TRUE
#' )
#'
#' utils::head(
#'   as.data.frame(
#'     lot_ots(fitted)
#'   )
#' )
#' @rdname lot_map_to_ots
#' @export
methods::setMethod(
  "lot_map_to_ots",
  "LongOmicsTraj",
  function(
    object,
    visits,
    overwrite = FALSE,
    ...
  ) {
    
    if (
      nrow(object@ots) > 0L &&
      !isTRUE(overwrite)
    ) {
      message(
        "OTS labels already exist. ",
        "Use `overwrite = TRUE` to recompute."
      )
      
      return(object)
    }
    
    if (nrow(object@deltas) == 0L) {
      stop(
        "No deltas found. Run `lot_compute_deltas()` before ",
        "`lot_map_to_ots()`.",
        call. = FALSE
      )
    }
    
    if (nrow(object@thresholds) == 0L) {
      stop(
        "No thresholds found. Run `lot_compute_thresholds()` before ",
        "`lot_map_to_ots()`.",
        call. = FALSE
      )
    }
    
    object@ots <- .compute_ots(
      deltas = object@deltas,
      thresholds = object@thresholds
    )
    
    object@history <- c(
      object@history,
      list(
        list(
          step = "lot_map_to_ots",
          time = Sys.time(),
          visits = visits
        )
      )
    )
    
    methods::validObject(
      object
    )
    
    object
  }
)