setMethod(
  "lot_from_empirical",
  "LongOmicsTraj",
  function(
    object,
    assay,
    group_col = "group",
    visit_col = "visit",
    visits,
    estimator = "empirical",
    overwrite = FALSE,
    verbose = TRUE
  ) {

    if (!assay %in% names(MultiAssayExperiment::experiments(object@mae))) {
      stop("Assay not found in object@mae: ", assay, call. = FALSE)
    }

    if (nrow(object@vmm) > 0 && !overwrite) {
      stop("VMM already exists. Use overwrite = TRUE to replace.", call. = FALSE)
    }

    se <- MultiAssayExperiment::experiments(object@mae)[[assay]]
    mat <- SummarizedExperiment::assay(se)
    meta <- as.data.frame(SummarizedExperiment::colData(se))

    required_meta <- c(group_col, visit_col)
    missing_meta <- setdiff(required_meta, colnames(meta))

    if (length(missing_meta) > 0) {
      stop(
        "colData is missing required columns: ",
        paste(missing_meta, collapse = ", "),
        call. = FALSE
      )
    }

    if (missing(visits) || is.null(visits) || length(visits) < 2) {
      stop("At least two ordered visits are required.", call. = FALSE)
    }

    if (is.null(rownames(mat))) {
      rownames(mat) <- paste0("feature_", seq_len(nrow(mat)))
    }

    meta <- meta[meta[[visit_col]] %in% visits, , drop = FALSE]

    if (nrow(meta) == 0) {
      stop(
        "No samples matched the supplied visits. Check visit labels.",
        call. = FALSE
      )
    }

    if (is.null(rownames(meta))) {
      stop("colData must have rownames matching assay column names.", call. = FALSE)
    }

    missing_samples <- setdiff(rownames(meta), colnames(mat))

    if (length(missing_samples) > 0) {
      stop(
        "Some colData rownames are missing from assay matrix column names: ",
        paste(head(missing_samples, 10), collapse = ", "),
        call. = FALSE
      )
    }

    mat <- mat[, rownames(meta), drop = FALSE]

    groups <- unique(as.character(meta[[group_col]]))

    out_list <- vector("list", length(groups) * length(visits))
    idx <- 1L

    for (grp in groups) {
      for (v in visits) {
        keep <- as.character(meta[[group_col]]) == grp &
          as.character(meta[[visit_col]]) == v

        if (!any(keep)) next

        means <- matrixStats::rowMeans2(
          as.matrix(mat[, keep, drop = FALSE]),
          na.rm = TRUE
        )

        out_list[[idx]] <- data.frame(
          object_id = rownames(mat),
          object_type = "gene",
          assay = assay,
          group = grp,
          visit = v,
          estimated_value = means,
          estimator = estimator,
          stringsAsFactors = FALSE
        )

        idx <- idx + 1L
      }
    }

    out_list <- out_list[!vapply(out_list, is.null, logical(1))]

    if (length(out_list) == 0) {
      stop(
        "No empirical VMM rows were created. Check that visits match colData visit labels and that each group/visit has samples.",
        call. = FALSE
      )
    }

    vmm_new <- S4Vectors::DataFrame(
      do.call(rbind, out_list)
    )

    if (overwrite || nrow(object@vmm) == 0) {
      object@vmm <- vmm_new
    } else {
      object@vmm <- S4Vectors::DataFrame(
        rbind(as.data.frame(object@vmm), as.data.frame(vmm_new))
      )
    }

    object@settings[[paste0("empirical_", assay)]] <- list(
      assay = assay,
      estimator = estimator,
      visits = visits,
      group_col = group_col,
      visit_col = visit_col,
      time = Sys.time()
    )

    object@history <- c(
      object@history,
      list(
        list(
          step = "lot_from_empirical",
          assay = assay,
          estimator = estimator,
          n_features = nrow(mat),
          n_vmm_rows = nrow(vmm_new),
          time = Sys.time()
        )
      )
    )

    if (verbose) {
      message(
        "Created empirical VMM for assay '", assay,
        "' with ", nrow(vmm_new), " rows."
      )
    }

    validObject(object)
    object
  }
)
