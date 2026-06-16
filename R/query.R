#' Negate a topology query condition
#'
#' @param x Character vector of topology labels to exclude.
#'
#' @return A query condition used by lot_query().
#' @export
NOT <- function(x) {
  if (!is.character(x)) {
    stop("NOT() expects a character vector.", call. = FALSE)
  }

  structure(
    list(exclude = x),
    class = "lot_not"
  )
}
.match_condition <- function(actual, condition) {
  if (length(actual) == 0 || is.na(actual[1])) return(FALSE)

  if (inherits(condition, "lot_not")) {
    return(!(actual %in% condition$exclude))
  }

  if (is.character(condition)) {
    return(actual %in% condition)
  }

  if (is.list(condition)) {
    ok <- TRUE

    if (!is.null(condition$include)) {
      ok <- ok && actual %in% condition$include
    }

    if (!is.null(condition$exclude)) {
      ok <- ok && !(actual %in% condition$exclude)
    }

    return(ok)
  }

  stop(
    "Invalid condition. Use character vector, NOT(), or list(include=, exclude=).",
    call. = FALSE
  )
}

.handle_cross_arm_query <- function(df, conditions) {
  required <- c(
    "object_id", "object_type", "assay",
    "group", "estimator", "topology_label"
  )

  missing_cols <- setdiff(required, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      "OTS table is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  wide <- stats::reshape(
    df[, c("object_id", "group", "topology_label")],
    idvar = "object_id",
    timevar = "group",
    direction = "wide"
  )

  keep <- rep(TRUE, nrow(wide))

  for (arm in names(conditions)) {
    col <- paste0("topology_label.", arm)

    if (!col %in% colnames(wide)) {
      keep <- rep(FALSE, nrow(wide))
      next
    }

    keep <- keep & vapply(
      wide[[col]],
      .match_condition,
      logical(1),
      condition = conditions[[arm]]
    )
  }

  ids <- wide$object_id[keep]

  df[df$object_id %in% ids, , drop = FALSE]
}

setMethod(
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

    df <- as.data.frame(object@ots)

    if (nrow(df) == 0) {
      stop("No OTS labels found. Run lot_map_to_ots() first.", call. = FALSE)
    }

    if (!is.null(object_type)) {
      df <- df[df$object_type %in% object_type, , drop = FALSE]
    }

    if (!is.null(assay)) {
      df <- df[df$assay %in% assay, , drop = FALSE]
    }

    if (!is.null(estimator)) {
      df <- df[df$estimator %in% estimator, , drop = FALSE]
    }

    if (!is.null(group)) {
      df <- df[df$group %in% group, , drop = FALSE]
    }

    if (!is.null(topology)) {
      df <- df[df$topology_label %in% topology, , drop = FALSE]
    }

    if (!is.null(conditions)) {
      df <- .handle_cross_arm_query(df, conditions)
    }

    S4Vectors::DataFrame(df)
  }
)
