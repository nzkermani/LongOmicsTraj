#' @rdname lot_trajectory_estimators
#' @examples
#' \donttest{
#' if (
#'   requireNamespace("lme4", quietly = TRUE) &&
#'   requireNamespace("emmeans", quietly = TRUE) &&
#'   requireNamespace("BiocParallel", quietly = TRUE)
#' ) {
#'   data(lot_glucold)
#'
#'   lot_lmm <- lot_from_lmm(
#'     object = lot_glucold,
#'     assay = "transcriptomics",
#'     formula_fixed = y ~ group * visit + (1 | subject_id),
#'     visits = c(
#'       "Baseline",
#'       "6 Months",
#'       "30 Months"
#'     ),
#'     BPPARAM = BiocParallel::SerialParam(),
#'     overwrite = TRUE
#'   )
#'
#'   lot_lmm
#' }
#' }
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
    BPPARAM = NULL,
    overwrite = FALSE,
    verbose = TRUE
  ) {
    
    ## ------------------------------------------------------------------------
    ## Check required packages
    ## ------------------------------------------------------------------------
    
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop(
        "Package 'lme4' is required. Install it with ",
        "install.packages('lme4').",
        call. = FALSE
      )
    }
    
    if (!requireNamespace("emmeans", quietly = TRUE)) {
      stop(
        "Package 'emmeans' is required. Install it with ",
        "install.packages('emmeans').",
        call. = FALSE
      )
    }
    
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop(
        "Package 'BiocParallel' is required. Install it with ",
        "BiocManager::install('BiocParallel').",
        call. = FALSE
      )
    }
    
    ## Use serial processing unless the user supplies another backend
    if (is.null(BPPARAM)) {
      BPPARAM <- BiocParallel::SerialParam()
    }
    
    
    ## ------------------------------------------------------------------------
    ## Validate assay and existing results
    ## ------------------------------------------------------------------------
    
    if (!assay %in%
        names(MultiAssayExperiment::experiments(object@mae))) {
      stop(
        "Assay not found in object@mae: ",
        assay,
        call. = FALSE
      )
    }
    
    if (nrow(object@vmm) > 0L && !overwrite) {
      stop(
        "VMM already exists. Use overwrite = TRUE to replace.",
        call. = FALSE
      )
    }
    
    
    ## ------------------------------------------------------------------------
    ## Extract assay and metadata
    ## ------------------------------------------------------------------------
    
    se <- MultiAssayExperiment::experiments(object@mae)[[assay]]
    
    assay_names <- SummarizedExperiment::assayNames(se)
    
    if (length(assay_names) == 0L) {
      stop(
        "The selected SummarizedExperiment contains no assay data.",
        call. = FALSE
      )
    }
    
    ## Use the first assay matrix
    mat <- SummarizedExperiment::assay(
      se,
      assay_names[1]
    )
    
    meta <- as.data.frame(
      SummarizedExperiment::colData(se)
    )
    
    
    ## ------------------------------------------------------------------------
    ## Validate metadata
    ## ------------------------------------------------------------------------
    
    required_meta <- unique(
      c(
        subject_col,
        group_col,
        visit_col
      )
    )
    
    missing_meta <- setdiff(
      required_meta,
      colnames(meta)
    )
    
    if (length(missing_meta) > 0L) {
      stop(
        "colData is missing required columns: ",
        paste(missing_meta, collapse = ", "),
        call. = FALSE
      )
    }
    
    if (missing(visits) || length(visits) < 2L) {
      stop(
        "'visits' must contain at least two ordered visit labels.",
        call. = FALSE
      )
    }
    
    if (is.null(rownames(meta))) {
      stop(
        "colData must have sample identifiers as row names.",
        call. = FALSE
      )
    }
    
    if (is.null(colnames(mat))) {
      stop(
        "The assay matrix must have sample identifiers as column names.",
        call. = FALSE
      )
    }
    
    if (is.null(rownames(mat))) {
      rownames(mat) <- paste0(
        "feature_",
        seq_len(nrow(mat))
      )
    }
    
    
    ## ------------------------------------------------------------------------
    ## Restrict to requested visits and align samples
    ## ------------------------------------------------------------------------
    
    meta <- meta[
      meta[[visit_col]] %in% visits,
      ,
      drop = FALSE
    ]
    
    if (nrow(meta) == 0L) {
      stop(
        "No samples remained after filtering to the requested visits.",
        call. = FALSE
      )
    }
    
    meta[[visit_col]] <- factor(
      meta[[visit_col]],
      levels = visits,
      ordered = TRUE
    )
    
    missing_samples <- setdiff(
      rownames(meta),
      colnames(mat)
    )
    
    if (length(missing_samples) > 0L) {
      stop(
        "Some colData row names are missing from the assay matrix: ",
        paste(
          head(missing_samples, 10),
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    
    mat <- mat[
      ,
      rownames(meta),
      drop = FALSE
    ]
    
    feature_ids <- rownames(mat)
    
    
    ## ------------------------------------------------------------------------
    ## Fit one LMM per feature
    ## ------------------------------------------------------------------------
    
    fit_one <- function(i) {
      
      df <- meta
      df$y <- as.numeric(mat[i, ])
      
      ## Remove observations with missing response values
      keep <- is.finite(df$y)
      
      df <- df[
        keep,
        ,
        drop = FALSE
      ]
      
      if (nrow(df) == 0L) {
        return(NULL)
      }
      
      fit <- tryCatch(
        lme4::lmer(
          formula = formula_fixed,
          data = df,
          REML = FALSE
        ),
        error = function(e) {
          if (verbose) {
            message(
              "LMM failed for feature ",
              feature_ids[i],
              ": ",
              conditionMessage(e)
            )
          }
          
          NULL
        }
      )
      
      if (is.null(fit)) {
        return(NULL)
      }
      
      ## Try estimated marginal means first
      emm <- tryCatch(
        as.data.frame(
          emmeans::emmeans(
            fit,
            specs = stats::as.formula(
              paste(
                "~",
                group_col,
                "*",
                visit_col
              )
            )
          )
        ),
        error = function(e) {
          
          ## Fallback to fixed-effect predictions
          design <- unique(
            df[
              ,
              c(
                group_col,
                visit_col
              ),
              drop = FALSE
            ]
          )
          
          design[[visit_col]] <- factor(
            design[[visit_col]],
            levels = visits,
            ordered = TRUE
          )
          
          preds <- tryCatch(
            stats::predict(
              fit,
              newdata = design,
              re.form = NA,
              allow.new.levels = TRUE
            ),
            error = function(e2) NULL
          )
          
          if (is.null(preds)) {
            return(NULL)
          }
          
          design$emmean <- as.numeric(preds)
          
          design
        }
      )
      
      if (is.null(emm) || !"emmean" %in% colnames(emm)) {
        return(NULL)
      }
      
      data.frame(
        object_id = feature_ids[i],
        object_type = "gene",
        assay = assay,
        group = as.character(
          emm[[group_col]]
        ),
        visit = as.character(
          emm[[visit_col]]
        ),
        estimated_value = as.numeric(
          emm$emmean
        ),
        estimator = estimator,
        stringsAsFactors = FALSE
      )
    }
    
    
    ## ------------------------------------------------------------------------
    ## Run feature-level models
    ## ------------------------------------------------------------------------
    
    res <- BiocParallel::bplapply(
      seq_along(feature_ids),
      fit_one,
      BPPARAM = BPPARAM
    )
    
    res <- res[
      !vapply(
        res,
        is.null,
        logical(1)
      )
    ]
    
    if (length(res) == 0L) {
      stop(
        "No LMM models were successfully fitted.",
        call. = FALSE
      )
    }
    
    
    ## ------------------------------------------------------------------------
    ## Store VMM
    ## ------------------------------------------------------------------------
    
    vmm_new <- S4Vectors::DataFrame(
      do.call(
        rbind,
        res
      )
    )
    
    if (overwrite || nrow(object@vmm) == 0L) {
      
      object@vmm <- vmm_new
      
    } else {
      
      object@vmm <- S4Vectors::DataFrame(
        rbind(
          as.data.frame(object@vmm),
          as.data.frame(vmm_new)
        )
      )
    }
    
    
    ## ------------------------------------------------------------------------
    ## Record provenance
    ## ------------------------------------------------------------------------
    
    object@history <- c(
      object@history,
      list(
        list(
          step = "lot_from_lmm",
          assay = assay,
          estimator = estimator,
          formula = paste(
            deparse(formula_fixed),
            collapse = " "
          ),
          visits = visits,
          n_features = length(feature_ids),
          n_success = length(res),
          BPPARAM = class(BPPARAM)[1],
          time = Sys.time()
        )
      )
    )
    
    if (verbose) {
      message(
        "Created LMM VMM for assay '",
        assay,
        "' with ",
        nrow(vmm_new),
        " rows from ",
        length(res),
        " successfully fitted features."
      )
    }
    
    methods::validObject(object)
    
    object
  }
)