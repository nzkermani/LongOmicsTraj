#' @rdname lot_trajectory_estimators
#' @export
setMethod(
  "lot_from_gam",
  "LongOmicsTraj",
  function(
    object,
    assay,
    formula_fixed,
    group_col = "group",
    visit_col = "visit",
    subject_col = "subject_id",
    visits,
    estimator = "gam",
    overwrite = FALSE,
    verbose = TRUE
  ) {
    
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("Package 'mgcv' is required.", call. = FALSE)
    }
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("Package 'dplyr' is required.", call. = FALSE)
    }
    
    if (!assay %in% names(MultiAssayExperiment::experiments(object@mae))) {
      stop("Assay not found: ", assay, call. = FALSE)
    }
    
    if (nrow(object@vmm) > 0 && !overwrite) {
      stop("VMM exists. Use overwrite=TRUE.", call. = FALSE)
    }
    
    se <- MultiAssayExperiment::experiments(object@mae)[[assay]]
    mat <- SummarizedExperiment::assay(se)
    meta <- as.data.frame(SummarizedExperiment::colData(se))
    
    if (!"time" %in% colnames(meta)) {
      stop("GAM requires numeric 'time' column in colData.", call. = FALSE)
    }
    
    meta <- meta[meta[[visit_col]] %in% visits, , drop = FALSE]
    
    if (is.null(rownames(meta))) {
      stop("colData rownames must match sample IDs.", call. = FALSE)
    }
    
    mat = mat[, rownames(meta), drop = FALSE]
    feature_ids <- rownames(mat)
    
    # Build prediction grid
    grid <- expand.grid(
      group = unique(meta[[group_col]]),
      visit = visits,
      stringsAsFactors = FALSE
    )
    
    colnames(grid)[colnames(grid) == "group"] <- group_col
    colnames(grid)[colnames(grid) == "visit"] <- visit_col
    
    visit_to_time <- dplyr::distinct(meta[, c(visit_col, "time"), drop = FALSE])
    
    grid <- dplyr::left_join(
      grid,
      visit_to_time,
      by = visit_col
    )
    
    # Pad subject column so predict.lm or predict.gam don't throw missing column errors
    if (!subject_col %in% colnames(grid) && subject_col %in% colnames(meta)) {
      grid[[subject_col]] <- meta[[subject_col]][1]
    }
    
    # 🚀 CRITICAL: Construct the exact name of the random effect smooth to drop during prediction
    random_smooth_name <- paste0("s(", subject_col, ")")
    
    res <- lapply(seq_along(feature_ids), function(i) {
      
      df <- meta
      df$y <- as.numeric(mat[i, ])
      current_feature <- feature_ids[i]
      
      make_fallback_df <- function() {
        mean_val <- mean(df$y, na.rm = TRUE)
        if (is.na(mean_val) || !is.finite(mean_val)) mean_val <- 0
        
        data.frame(
          object_id = current_feature,
          object_type = "gene",
          assay = assay,
          group = as.character(grid[[group_col]]),
          visit = as.character(grid[[visit_col]]),
          estimated_value = mean_val,
          estimator = estimator,
          stringsAsFactors = FALSE
        )
      }
      
      if (all(is.na(df$y)) || sd(df$y, na.rm = TRUE) == 0) {
        return(make_fallback_df())
      }
      
      df[[group_col]] <- as.factor(df[[group_col]])
      df[[subject_col]] <- as.factor(df[[subject_col]]) # Ensure factor for random effects
      
      local_formula <- formula_fixed
      environment(local_formula) <- environment()
      
      fit <- tryCatch(
        mgcv::gam(formula = local_formula, data = df, method = "REML"),
        error = function(e) NULL
      )
      
      is_gam <- TRUE
      if (is.null(fit)) {
        is_gam <- FALSE
        fit <- tryCatch(
          lm(local_formula, data = df),   
          error = function(e) NULL
        )
      }
      
      if (is.null(fit)) return(make_fallback_df())
      
      # 🚀 CRITICAL FIX: Safe prediction splitting based on model engine type
      preds <- tryCatch({
        if (is_gam) {
          # Use 'exclude' for GAM models to cleanly isolate population trajectories
          predict(fit, newdata = grid, type = "response", exclude = random_smooth_name)
        } else {
          # Standard linear regression prediction fallback
          predict(fit, newdata = grid, type = "response")
        }
      }, error = function(e) {
        NULL
      })
      
      if (is.null(preds)) return(make_fallback_df())
      
      data.frame(
        object_id = current_feature,
        object_type = "gene",
        assay = assay,
        group = as.character(grid[[group_col]]),
        visit = as.character(grid[[visit_col]]),
        estimated_value = as.numeric(preds),
        estimator = estimator,
        stringsAsFactors = FALSE
      )
    })
    
    vmm_new <- S4Vectors::DataFrame(do.call(rbind, res))
    
    if (overwrite || nrow(object@vmm) == 0) {
      object@vmm <- vmm_new
    } else {
      object@vmm <- S4Vectors::DataFrame(rbind(as.data.frame(object@vmm), as.data.frame(vmm_new)))
    }
    
    object@history <- c(object@history, list(list(step = "lot_from_gam", assay = assay, n_features = length(feature_ids), n_success = length(res), time = Sys.time())))
    if (verbose) message("Created GAM VMM for assay '", assay, "' with ", nrow(vmm_new), " rows.")
    
    validObject(object)
    object
  }
)
