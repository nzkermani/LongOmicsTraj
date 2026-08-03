#' Compute molecular-change thresholds
#'
#' Computes feature-specific thresholds used to classify molecular changes
#' between consecutive visits as down, flat or up.
#'
#' Thresholds are calculated from the delta estimates already stored in the
#' `LongOmicsTraj` object. The function supports fixed, global MAD,
#' feature-specific MAD, hybrid and quantile-based thresholding.
#'
#' @param object A `LongOmicsTraj` object.
#' @param assay Name of the assay for which thresholds should be computed.
#' @param method Threshold-estimation method. One of `"fixed"`,
#'   `"global_mad"`, `"gene_mad"`, `"hybrid"` or `"quantile"`.
#' @param delta_min Optional minimum threshold value. When `NULL`, the median
#'   absolute delta across all objects and transitions is used.
#' @param k Numeric scaling factor applied to MAD-based thresholds.
#' @param prob Probability used for quantile-based threshold estimation.
#' @param visits Character vector giving the visits in chronological order.
#'   When omitted, `object@settings$visits` is used.
#' @param overwrite Logical; overwrite existing thresholds for the selected
#'   assay when `TRUE`.
#' @param verbose Logical; print progress messages.
#'
#' @return A `LongOmicsTraj` object with updated threshold results.
#'
#' @details
#' Delta values are grouped by molecular object and converted into a matrix
#' of adjacent-visit changes. The available thresholding methods are:
#'
#' \describe{
#'   \item{`fixed`}{
#'     Uses the same `delta_min` threshold for every object.
#'   }
#'   \item{`global_mad`}{
#'     Uses a single threshold based on the global median MAD.
#'   }
#'   \item{`gene_mad`}{
#'     Uses an object-specific MAD threshold.
#'   }
#'   \item{`hybrid`}{
#'     Uses the maximum of `delta_min` and the object-specific MAD threshold.
#'   }
#'   \item{`quantile`}{
#'     Uses the maximum of `delta_min` and an object-specific absolute-delta
#'     quantile.
#'   }
#' }
#'
#' The function supports both gene-level and cluster-level delta results.
#' Object type is inferred from identifiers beginning with `"cluster_"`.
#' @examples
#' data("lot_glucold_ots")
#'
#' visits <- c(
#'   "Baseline",
#'   "6 Months",
#'   "30 Months"
#' )
#'
#' fitted <- lot_compute_thresholds(
#'   object = lot_glucold_ots,
#'   assay = "transcriptomics",
#'   method = "hybrid",
#'   visits = visits,
#'   overwrite = TRUE
#' )
#'
#' thresholds <- as.data.frame(
#'   lot_thresholds(fitted)
#' )
#'
#' utils::head(
#'   thresholds[
#'     ,
#'     c(
#'       "object_id",
#'       "delta_g",
#'       "threshold_method"
#'     )
#'   ]
#' )
#' @rdname lot_compute_thresholds
#' @export
methods::setMethod(
  "lot_compute_thresholds",
  "LongOmicsTraj",
  function(
    object,
    assay,
    method = c(
      "fixed",
      "global_mad",
      "gene_mad",
      "hybrid",
      "quantile"
    ),
    delta_min = NULL,
    k = 0.25,
    prob = 0.75,
    visits,
    overwrite = FALSE,
    verbose = TRUE
  ) {
    
    method <- match.arg(method)
    
    if (missing(visits)) {
      visits <- object@settings$visits
    }
    
    if (
      is.null(visits) ||
      length(visits) < 2L
    ) {
      stop(
        "At least two ordered visits are required.",
        call. = FALSE
      )
    }
    
    old_thresholds <- as.data.frame(
      object@thresholds
    )
    
    if (
      nrow(old_thresholds) > 0L &&
      "assay" %in% colnames(old_thresholds) &&
      assay %in% old_thresholds$assay &&
      !isTRUE(overwrite)
    ) {
      stop(
        "Thresholds for assay '",
        assay,
        "' already exist. Use `overwrite = TRUE` to recompute.",
        call. = FALSE
      )
    }
    
    # -----------------------------------------------------------------------
    # Delta-driven threshold calculation
    # -----------------------------------------------------------------------
    
    deltas <- as.data.frame(
      object@deltas
    )
    
    if (nrow(deltas) == 0L) {
      stop(
        "No deltas found. Run `lot_compute_deltas()` first.",
        call. = FALSE
      )
    }
    
    split_deltas <- split(
      deltas$delta_hat,
      deltas$object_id
    )
    
    diff_mat <- do.call(
      rbind,
      split_deltas
    )
    
    diff_mat <- as.matrix(
      diff_mat
    )
    
    rownames(diff_mat) <- names(
      split_deltas
    )
    
    abs_diff_mat <- abs(
      diff_mat
    )
    
    # -----------------------------------------------------------------------
    # Threshold statistics
    # -----------------------------------------------------------------------
    
    tau_g <- matrixStats::rowMads(
      diff_mat,
      constant = 1.4826,
      na.rm = TRUE
    )
    
    q_g <- matrixStats::rowQuantiles(
      abs_diff_mat,
      probs = prob,
      na.rm = TRUE,
      drop = TRUE
    )
    
    delta_min_0 <- stats::median(
      abs_diff_mat,
      na.rm = TRUE
    )
    
    if (is.null(delta_min)) {
      delta_min <- delta_min_0
    }
    
    tau_baseline <- stats::median(
      tau_g[
        is.finite(tau_g) &
          tau_g > 0
      ],
      na.rm = TRUE
    )
    
    q_baseline <- stats::median(
      q_g[
        is.finite(q_g) &
          q_g > 0
      ],
      na.rm = TRUE
    )
    
    if (!is.finite(tau_baseline)) {
      tau_baseline <- 0
    }
    
    if (!is.finite(q_baseline)) {
      q_baseline <- 0
    }
    
    tau_adjusted <- pmax(
      tau_g,
      tau_baseline,
      na.rm = TRUE
    )
    
    tau_adjusted[
      !is.finite(tau_adjusted)
    ] <- tau_baseline
    
    q_adjusted <- pmax(
      q_g,
      q_baseline,
      na.rm = TRUE
    )
    
    q_adjusted[
      !is.finite(q_adjusted)
    ] <- q_baseline
    
    global_mad_delta <- max(
      delta_min,
      k * tau_baseline
    )
    
    delta_g <- switch(
      method,
      fixed = rep(
        delta_min,
        nrow(diff_mat)
      ),
      global_mad = rep(
        global_mad_delta,
        nrow(diff_mat)
      ),
      gene_mad = k * tau_adjusted,
      hybrid = pmax(
        delta_min,
        k * tau_adjusted
      ),
      quantile = pmax(
        delta_min,
        q_adjusted
      )
    )
    
    # -----------------------------------------------------------------------
    # Build threshold results
    # -----------------------------------------------------------------------
    
    object_ids <- rownames(
      diff_mat
    )
    
    object_type <- ifelse(
      grepl(
        "^cluster_",
        object_ids
      ),
      "cluster",
      "gene"
    )
    
    new_thresholds <- S4Vectors::DataFrame(
      object_id = object_ids,
      object_type = object_type,
      assay = assay,
      threshold_method = method,
      n_adjacent_pairs = ncol(diff_mat),
      tau_g = tau_g,
      tau_baseline = tau_baseline,
      tau_adjusted = tau_adjusted,
      q_abs_delta = q_g,
      q_baseline = q_baseline,
      q_adjusted = q_adjusted,
      delta_min = delta_min,
      delta_min_0 = delta_min_0,
      k = k,
      prob = prob,
      delta_g = delta_g
    )
    
    # -----------------------------------------------------------------------
    # Assign threshold results
    # -----------------------------------------------------------------------
    
    if (
      nrow(old_thresholds) > 0L &&
      isTRUE(overwrite)
    ) {
      old_thresholds <- old_thresholds[
        old_thresholds$assay != assay,
        ,
        drop = FALSE
      ]
    }
    
    if (nrow(old_thresholds) > 0L) {
      object@thresholds <- S4Vectors::DataFrame(
        rbind(
          old_thresholds,
          as.data.frame(new_thresholds)
        )
      )
    } else {
      object@thresholds <- new_thresholds
    }
    
    # -----------------------------------------------------------------------
    # Store settings and history
    # -----------------------------------------------------------------------
    
    object@settings[[
      paste0(
        "threshold_",
        assay
      )
    ]] <- list(
      method = method,
      delta_min = delta_min,
      delta_min_0 = delta_min_0,
      k = k,
      prob = prob,
      visits = visits,
      time = Sys.time()
    )
    
    object@history <- c(
      object@history,
      list(
        list(
          step = "lot_compute_thresholds",
          assay = assay,
          method = method,
          n_features = nrow(diff_mat),
          n_adjacent_pairs = ncol(diff_mat),
          time = Sys.time()
        )
      )
    )
    
    if (isTRUE(verbose)) {
      message(
        "Computed ",
        method,
        " thresholds for assay '",
        assay,
        "' using ",
        nrow(diff_mat),
        " objects and ",
        ncol(diff_mat),
        " transitions."
      )
    }
    
    methods::validObject(
      object
    )
    
    object
  }
)