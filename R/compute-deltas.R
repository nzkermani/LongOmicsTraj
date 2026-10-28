#' Compute between-visit molecular changes
#'
#' Calculates feature-level changes between consecutive visits from the
#' visit-level molecular trajectories stored in a `LongOmicsTraj` object.
#'
#' @param object A `LongOmicsTraj` object.
#' @param assay Name of the assay to process.
#' @param visits Character vector giving the visits in chronological order.
#' @param overwrite Logical; overwrite existing delta results when `TRUE`.
#' @param center Character string controlling between-visit centring of the
#'   visit-level trajectories before changes are calculated. `"none"` (the
#'   default) uses the trajectories as estimated. `"median"` or `"mean"`
#'   subtracts, within each group and visit, the median or mean of the
#'   estimated values across all features of the assay. Centring removes
#'   visit-wide offsets that affect every feature in the same direction,
#'   such as differences in capture efficiency or tissue composition between
#'   ordered spatial windows, so that topology labels describe relative
#'   rather than global change.
#'
#' @return A `LongOmicsTraj` object with updated delta results. When
#'   centring is applied, the offset removed for each group and visit is
#'   stored in `results$visit_centring`.
#'
#' @details
#' Centring assumes that most features do not change between consecutive
#' visits. It is intended for settings where visit-wide shifts are technical
#' or compositional, for example depth windows in spatial omics. It should not
#' be used when a genuine global change across most features is expected.
#'
#' @examples
#' data("lot_glucold_ots")
#'
#' visits <- c(
#'   "Baseline",
#'   "6 Months",
#'   "30 Months"
#' )
#'
#' fitted <- lot_compute_deltas(
#'   object = lot_glucold_ots,
#'   assay = "transcriptomics",
#'   visits = visits,
#'   overwrite = TRUE
#' )
#'
#' deltas <- as.data.frame(
#'   lot_deltas(fitted)
#' )
#'
#' utils::head(
#'   deltas[
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
#' # Remove visit-wide offsets before computing changes
#' centred <- lot_compute_deltas(
#'   object = lot_glucold_ots,
#'   assay = "transcriptomics",
#'   visits = visits,
#'   overwrite = TRUE,
#'   center = "median"
#' )
#'
#' as.data.frame(centred@results$visit_centring)
#'
#' @rdname lot_compute_deltas
#' @export
setMethod(
  "lot_compute_deltas",
  "LongOmicsTraj",
  function(
    object,
    assay,
    visits,
    overwrite = FALSE,
    center = c("none", "median", "mean")
  ) {

    center <- match.arg(center)

    if (!overwrite && nrow(object@deltas) > 0) {
      stop("Deltas already exist. Use overwrite=TRUE.", call. = FALSE)
    }

    vmm <- as.data.frame(object@vmm)

    if (nrow(vmm) == 0) {
      stop("VMM is empty. Run a lot_from_* method first.", call. = FALSE)
    }

    # Optional between-visit centring (removes visit-wide offsets)
    if (center != "none") {
      fun <- if (center == "median") stats::median else base::mean
      key <- paste(vmm$group, vmm$visit, sep = "\r")
      offset <- tapply(vmm$estimated_value, key, fun, na.rm = TRUE)
      vmm$estimated_value <- vmm$estimated_value - offset[key]
      parts <- strsplit(names(offset), "\r", fixed = TRUE)
      object@results$visit_centring <- S4Vectors::DataFrame(
        group = vapply(parts, `[`, character(1), 1L),
        visit = vapply(parts, `[`, character(1), 2L),
        method = center,
        offset = as.numeric(offset)
      )
    }

    # Ensure visit ordering
    vmm$visit <- factor(vmm$visit, levels = visits, ordered = TRUE)

    # Order properly
    vmm <- vmm[order(vmm$object_id, vmm$group, vmm$visit), ]

    # Compute deltas
    delta_list <- split(vmm, list(vmm$object_id, vmm$group), drop = TRUE)

    res <- lapply(delta_list, function(df) {

      if (nrow(df) < 2) return(NULL)

      df <- df[order(df$visit), ]

      data.frame(
        object_id = df$object_id[-1],
        object_type = df$object_type[-1],
        assay = df$assay[-1],
        group = df$group[-1],
        estimator = df$estimator[-1],
        transition_index = seq_len(nrow(df) - 1),
        from_visit = df$visit[-nrow(df)],
        to_visit = df$visit[-1],
        transition = paste(df$visit[-nrow(df)], df$visit[-1], sep = "_"),
        delta_hat = diff(df$estimated_value),
        stringsAsFactors = FALSE
      )
    })

    deltas <- do.call(rbind, res)

    object@deltas <- S4Vectors::DataFrame(deltas)

    object
  }
)
