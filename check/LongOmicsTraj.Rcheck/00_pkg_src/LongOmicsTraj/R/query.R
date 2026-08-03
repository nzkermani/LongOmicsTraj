#' Negate a topology query condition
#'
#' Creates an exclusion condition for use in `lot_query()`. Features whose
#' topology labels match any value in `x` are excluded.
#'
#' @param x Character vector of topology labels to exclude.
#'
#' @return An object of class `"lot_not"` representing an exclusion condition.
#'
#' @rdname lot_query
#' @export
NOT <- function(x) {
  
  if (!is.character(x)) {
    stop(
      "`NOT()` expects a character vector.",
      call. = FALSE
    )
  }
  
  structure(
    list(
      exclude = x
    ),
    class = "lot_not"
  )
}


#' Match a query condition internally
#'
#' Tests whether an observed topology label satisfies an inclusion,
#' exclusion or combined query condition.
#'
#' @param actual Observed topology label.
#' @param condition A character vector, a `"lot_not"` object or a list
#'   containing `include` and/or `exclude` elements.
#'
#' @return A single logical value.
#'
#' @keywords internal
.match_condition <- function(
    actual,
    condition
) {
  
  if (
    length(actual) == 0L ||
    is.na(actual[[1L]])
  ) {
    return(FALSE)
  }
  
  if (inherits(condition, "lot_not")) {
    return(
      !(actual %in% condition$exclude)
    )
  }
  
  if (is.character(condition)) {
    return(
      actual %in% condition
    )
  }
  
  if (is.list(condition)) {
    
    ok <- TRUE
    
    if (!is.null(condition$include)) {
      ok <- ok &&
        actual %in% condition$include
    }
    
    if (!is.null(condition$exclude)) {
      ok <- ok &&
        !(actual %in% condition$exclude)
    }
    
    return(ok)
  }
  
  stop(
    "Invalid condition. Use a character vector, `NOT()`, or ",
    "`list(include = ..., exclude = ...)`.",
    call. = FALSE
  )
}


#' Handle cross-group topology queries internally
#'
#' Evaluates topology conditions across multiple groups for each molecular
#' feature.
#'
#' @param df A data frame containing OTS assignments.
#' @param conditions A named list of query conditions. Names must correspond
#'   to group labels.
#'
#' @return A filtered data frame containing features satisfying all requested
#' group-specific conditions.
#'
#' @keywords internal
.handle_cross_arm_query <- function(
    df,
    conditions
) {
  
  required <- c(
    "object_id",
    "object_type",
    "assay",
    "group",
    "estimator",
    "topology_label"
  )
  
  missing_cols <- setdiff(
    required,
    colnames(df)
  )
  
  if (length(missing_cols) > 0L) {
    stop(
      "OTS table is missing required columns: ",
      paste(
        missing_cols,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  if (
    !is.list(conditions) ||
    is.null(names(conditions)) ||
    any(!nzchar(names(conditions)))
  ) {
    stop(
      "`conditions` must be a named list with group names as names.",
      call. = FALSE
    )
  }
  
  wide <- stats::reshape(
    df[
      ,
      c(
        "object_id",
        "group",
        "topology_label"
      ),
      drop = FALSE
    ],
    idvar = "object_id",
    timevar = "group",
    direction = "wide"
  )
  
  keep <- rep(
    TRUE,
    nrow(wide)
  )
  
  for (arm in names(conditions)) {
    
    column_name <- paste0(
      "topology_label.",
      arm
    )
    
    if (!column_name %in% colnames(wide)) {
      keep <- rep(
        FALSE,
        nrow(wide)
      )
      
      next
    }
    
    keep <- keep &
      vapply(
        wide[[column_name]],
        .match_condition,
        logical(1),
        condition = conditions[[arm]]
      )
  }
  
  ids <- wide$object_id[
    keep
  ]
  
  df[
    df$object_id %in% ids,
    ,
    drop = FALSE
  ]
}


#' Query Omics Trajectory Signature results
#'
#' Filters Omics Trajectory Signature assignments stored in a
#' `LongOmicsTraj` object.
#'
#' Queries may be restricted by molecular object type, assay, group,
#' estimator or topology label. The `conditions` argument supports
#' group-specific cross-arm queries.
#'
#' @param object A `LongOmicsTraj` object.
#' @param object_type Optional character vector of molecular object types.
#' @param assay Optional character vector of assays.
#' @param group Optional character vector of clinical or treatment groups.
#' @param topology Optional character vector of OTS topology labels.
#' @param conditions Optional named list defining group-specific topology
#'   conditions. Each element may be a character vector, a `NOT()` condition,
#'   or a list with `include` and/or `exclude` elements.
#' @param estimator Optional character vector of trajectory estimators.
#'
#' @return A `DataFrame` containing the filtered OTS assignments.
#'
#' @details
#' Simple filters are applied directly to the OTS table. Cross-group queries
#' supplied through `conditions` require the same molecular feature to satisfy
#' all named group-specific conditions.
#' @examples
#' data("lot_glucold_ots")
#'
#' flat <- lot_query(
#'   object = lot_glucold_ots,
#'   topology = "flat_flat"
#' )
#'
#' nrow(flat)
#'
#' utils::head(
#'   as.data.frame(flat)
#' )
#' @rdname lot_query
#' @export
methods::setMethod(
  "lot_query",
  "LongOmicsTraj",
  function(
    object,
    object_type = NULL,
    assay = NULL,
    group = NULL,
    topology = NULL,
    conditions = NULL,
    estimator = NULL
  ) {
    
    df <- as.data.frame(
      object@ots
    )
    
    if (nrow(df) == 0L) {
      stop(
        "No OTS labels found. Run `lot_map_to_ots()` first.",
        call. = FALSE
      )
    }
    
    if (!is.null(object_type)) {
      df <- df[
        df$object_type %in% object_type,
        ,
        drop = FALSE
      ]
    }
    
    if (!is.null(assay)) {
      df <- df[
        df$assay %in% assay,
        ,
        drop = FALSE
      ]
    }
    
    if (!is.null(estimator)) {
      df <- df[
        df$estimator %in% estimator,
        ,
        drop = FALSE
      ]
    }
    
    if (!is.null(group)) {
      df <- df[
        df$group %in% group,
        ,
        drop = FALSE
      ]
    }
    
    if (!is.null(topology)) {
      df <- df[
        df$topology_label %in% topology,
        ,
        drop = FALSE
      ]
    }
    
    if (!is.null(conditions)) {
      df <- .handle_cross_arm_query(
        df = df,
        conditions = conditions
      )
    }
    
    S4Vectors::DataFrame(
      df
    )
  }
)