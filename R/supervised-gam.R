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
    BPPARAM = BiocParallel::SerialParam(),
    overwrite = FALSE,
    verbose = TRUE
  ) {

    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("Package 'mgcv' is required.", call. = FALSE)
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

    # Ensure visit is numeric (GAM needs continuous time)
    if (!"time" %in% colnames(meta)) {
      stop("GAM requires numeric 'time' column in colData.", call. = FALSE)
    }

    meta <- meta[meta[[visit_col]] %in% visits, , drop = FALSE]

    if (is.null(rownames(meta))) {
      stop("colData rownames must match sample IDs.", call. = FALSE)
    }

    mat <- mat[, rownames(meta), drop = FALSE]

    feature_ids <- rownames(mat)

    # Build prediction grid
    grid <- expand.grid(
      group = unique(meta[[group_col]]),
      visit = visits,
      stringsAsFactors = FALSE
    )

    # Map visit → time
    visit_to_time <- meta %>%
      dplyr::select(!!visit_col, time) %>%
      distinct()

    grid <- dplyr::left_join(
      grid,
      visit_to_time,
      by = setNames("visit", visit_col)
    )

    fit_one <- function(i) {

      df <- meta
      df$y <- as.numeric(mat[i, ])

      #  remove bad genes early
      if (all(is.na(df$y)) || sd(df$y, na.rm = TRUE) == 0) {
        return(NULL)
      }

      # ensure factor
      df[[group_col]] <- as.factor(df[[group_col]])

      fit <- tryCatch(
        mgcv::gam(
          formula = formula_fixed,
          data = df,
          method = "REML"
        ),
        error = function(e) NULL
      )

      # fallback if GAM fails
      if (is.null(fit)) {
        fit <- tryCatch(
          lm(y ~ group + time, data = df),
          error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)
      }

      preds <- tryCatch(
        predict(fit, newdata = grid, type = "response"),
        error = function(e) NULL
      )

      if (is.null(preds)) return(NULL)

      data.frame(
        object_id = feature_ids[i],
        object_type = "gene",
        assay = assay,
        group = grid$group,
        visit = grid$visit,
        estimated_value = as.numeric(preds),
        estimator = estimator,
        stringsAsFactors = FALSE
      )
    }

    res <- BiocParallel::bplapply(
      seq_along(feature_ids),
      fit_one,
      BPPARAM = BPPARAM
    )

    res <- res[!vapply(res, is.null, logical(1))]

    if (length(res) == 0) {
      stop("No GAM models fitted.", call. = FALSE)
    }

    vmm_new <- S4Vectors::DataFrame(
      do.call(rbind, res)
    )

    if (overwrite || nrow(object@vmm) == 0) {
      object@vmm <- vmm_new
    } else {
      object@vmm <- S4Vectors::DataFrame(
        rbind(as.data.frame(object@vmm), as.data.frame(vmm_new))
      )
    }

    object@history <- c(
      object@history,
      list(
        list(
          step = "lot_from_gam",
          assay = assay,
          n_features = length(feature_ids),
          n_success = length(res),
          time = Sys.time()
        )
      )
    )

    if (verbose) {
      message(
        "Created GAM VMM for assay '", assay,
        "' with ", nrow(vmm_new), " rows."
      )
    }

    validObject(object)
    object
  }
)
