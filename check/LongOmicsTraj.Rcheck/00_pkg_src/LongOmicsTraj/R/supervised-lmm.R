#' @rdname lot_trajectory_estimators
#' @export
setMethod(
  "lot_from_lmm",
  "LongOmicsTraj",
  function(
    object,
    assay,
    formula_fixed,
    subject_col = "subject_id",
    group_col = "group",
    visit_col = "visit",
    visits,
    estimator = "lmm",
    BPPARAM = BiocParallel::SerialParam(),
    overwrite = FALSE,
    verbose = TRUE
  ) {

    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("Package 'lme4' is required.", call. = FALSE)
    }

    if (!requireNamespace("emmeans", quietly = TRUE)) {
      stop("Package 'emmeans' is required.", call. = FALSE)
    }

    if (!assay %in% names(MultiAssayExperiment::experiments(object@mae))) {
      stop("Assay not found in object@mae: ", assay, call. = FALSE)
    }

    if (nrow(object@vmm) > 0 && !overwrite) {
      stop("VMM already exists. Use overwrite = TRUE to replace.", call. = FALSE)
    }

    se <- MultiAssayExperiment::experiments(object@mae)[[assay]]
    mat <- SummarizedExperiment::assay(se)
    meta <- as.data.frame(SummarizedExperiment::colData(se))

    required_meta <- unique(c(subject_col, group_col, visit_col))
    missing_meta <- setdiff(required_meta, colnames(meta))

    if (length(missing_meta) > 0) {
      stop(
        "colData is missing required columns: ",
        paste(missing_meta, collapse = ", "),
        call. = FALSE
      )
    }

    if (is.null(rownames(mat))) {
      rownames(mat) <- paste0("feature_", seq_len(nrow(mat)))
    }

    meta <- meta[meta[[visit_col]] %in% visits, , drop = FALSE]
    meta[[visit_col]] <- factor(meta[[visit_col]], levels = visits, ordered = TRUE)

    missing_samples <- setdiff(rownames(meta), colnames(mat))

    if (length(missing_samples) > 0) {
      stop(
        "Some colData rownames are missing from assay matrix column names: ",
        paste(head(missing_samples, 10), collapse = ", "),
        call. = FALSE
      )
    }

    mat <- mat[, rownames(meta), drop = FALSE]

    feature_ids <- rownames(mat)

    fit_one <- function(i) {
      df <- meta
      df$y <- as.numeric(mat[i, ])

      fit <- tryCatch(
        lme4::lmer(
          formula = formula_fixed,
          data = df,
          REML = FALSE
        ),
        error = function(e) NULL
      )



      if (is.null(fit)) return(NULL)

      # fallback: if emmeans fails, use fixed-effect predictions
      emm <- tryCatch(
        as.data.frame(
          emmeans::emmeans(
            fit,
            specs = stats::as.formula(
              paste("~", group_col, "*", visit_col)
            )
          )
        ),
        error = function(e) {

          # fallback to model.matrix prediction
          design <- unique(df[, c(group_col, visit_col)])
          preds <- predict(fit, newdata = design, re.form = NA)

          design$emmean <- preds
          design
        }
      )

      if (is.null(emm)) return(NULL)

      data.frame(
        object_id = feature_ids[i],
        object_type = "gene",
        assay = assay,
        group = as.character(emm[[group_col]]),
        visit = as.character(emm[[visit_col]]),
        estimated_value = emm$emmean,
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
      stop("No LMM models were successfully fitted.", call. = FALSE)
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
          step = "lot_from_lmm",
          assay = assay,
          estimator = estimator,
          n_features = length(feature_ids),
          n_success = length(res),
          time = Sys.time()
        )
      )
    )

    if (verbose) {
      message(
        "Created LMM VMM for assay '", assay,
        "' with ", nrow(vmm_new), " rows."
      )
    }

    validObject(object)
    object
  }
)
