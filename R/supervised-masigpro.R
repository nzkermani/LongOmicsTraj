#' @export
setMethod(
  "lot_from_masigpro",
  "LongOmicsTraj",
  function(
    object,
    assay,
    group_col = "group",
    visit_col = "visit",
    time_col = "time",
    subject_col = "subject_id",
    visits,
    degree = NULL,
    estimator = "masigpro",
    overwrite = FALSE,
    verbose = TRUE
  ) {
    
    if (!assay %in% names(MultiAssayExperiment::experiments(object@mae))) {
      stop("Assay not found: ", assay, call. = FALSE)
    }
    
    if (nrow(object@vmm) > 0 && !overwrite) {
      stop("VMM exists. Use overwrite=TRUE.", call. = FALSE)
    }
    
    se <- MultiAssayExperiment::experiments(object@mae)[[assay]]
    mat <- as.matrix(SummarizedExperiment::assay(se))
    meta <- as.data.frame(SummarizedExperiment::colData(se))
    
    # Align columns and format structural types
    mat   <- mat[, rownames(meta), drop = FALSE]
    Time  <- as.numeric(meta[[time_col]])
    Group <- as.factor(meta[[group_col]])
    
    # ADAPTIVE LOGIC: Calculate degree dynamically based on unique target visits
    if (is.null(degree)) {
      n_unique_visits <- length(unique(Time))
      degree <- max(1, n_unique_visits - 1) # Floor at 1 for safety
    }
    
    n_groups <- length(levels(Group))
    
    if (verbose) message("Fitting polynomial trajectories via multivariate regression...")
    
    # Dynamic formula adaptation based on group variation
    if (n_groups > 1) {
      fit <- lm(t(mat) ~ poly(Time, degree, raw = TRUE) * Group)
    } else {
      # Drop Group contrast calculations if everyone is in the same group ("sim")
      fit <- lm(t(mat) ~ poly(Time, degree, raw = TRUE))
    }
    
    # Extract complete coefficients matrix (Genes x Parameters)
    coefs <- t(coef(fit))
    coefs[is.na(coefs)] <- 0
    feature_ids <- rownames(coefs)
    
    # Setup prediction target mapping grid
    grid <- expand.grid(
      group = unique(meta[[group_col]]),
      visit = visits,
      stringsAsFactors = FALSE
    )
    
    visit_map <- meta %>%
      dplyr::select(!!visit_col, !!time_col) %>%
      distinct()
    
    grid <- dplyr::left_join(grid, visit_map, by = setNames("visit", visit_col))
    
    # Build the matching prediction design matrix
    pred_df <- data.frame(
      Time = as.numeric(grid[[time_col]]),
      Group = factor(grid$group, levels = levels(Group))
    )
    
    if (n_groups > 1) {
      pred_matrix <- model.matrix(~ poly(Time, degree, raw = TRUE) * Group, data = pred_df)
    } else {
      pred_matrix <- model.matrix(~ poly(Time, degree, raw = TRUE), data = pred_df)
    }
    
    # Fast project trajectories
    preds <- coefs %*% t(pred_matrix)
    
    # Standardize back to LongOmicsTraj VMM format
    vmm_list <- lapply(seq_len(nrow(preds)), function(i) {
      data.frame(
        object_id = feature_ids[i],
        object_type = "gene",
        assay = assay,
        group = grid$group,
        visit = grid$visit,
        estimated_value = as.numeric(preds[i, ]),
        estimator = estimator,
        stringsAsFactors = FALSE
      )
    })
    
    vmm_new <- S4Vectors::DataFrame(do.call(rbind, vmm_list))
    object@vmm = vmm_new
    
    object@history <- c(
      object@history,
      list(list(step = "lot_from_masigpro", assay = assay, n_features = nrow(mat), time = Sys.time()))
    )
    
    if (verbose) {
      message("Created maSigPro-equivalent VMM for assay '", assay, "' with ", nrow(vmm_new), " rows.")
    }
    
    validObject(object)
    object
  }
)
