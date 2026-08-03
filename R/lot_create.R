#' Create a LongOmicsTraj object
#'
#' Creates a `LongOmicsTraj` object from a numeric molecular-data matrix
#' and sample-level clinical metadata.
#'
#' The input matrix should contain molecular features in rows and samples
#' in columns. The metadata should contain one row per sample. Sample IDs
#' are checked, aligned and stored in a `SummarizedExperiment`, which is
#' then placed inside a `MultiAssayExperiment`.
#'
#' @param expression A numeric matrix or matrix-like object. Rows represent
#'   molecular features and columns represent samples.
#' @param metadata A data frame containing one row per sample.
#' @param assay Character string giving the assay name. The default is
#'   `"transcriptomics"`.
#' @param sample_col Character string naming the sample identifier column
#'   in `metadata`. The default is `"sample_id"`.
#' @param subject_col Character string naming the participant identifier
#'   column in `metadata`. The default is `"subject_id"`.
#' @param visit_col Character string naming the visit column in `metadata`.
#'   The default is `"visit"`.
#' @param group_col Character string naming the clinical or treatment-group
#'   column in `metadata`. The default is `"group"`.
#' @param visit_levels Optional character vector giving the chronological
#'   order of visits. When supplied, the visit column is converted to an
#'   ordered factor.
#' @param time_values Optional named numeric vector mapping visit labels to
#'   numerical time. For example,
#'   `c("Baseline" = 0, "6 Months" = 6, "30 Months" = 30)`.
#' @param time_col Character string used for the numerical time column created
#'   from `time_values`. The default is `"time"`.
#' @param feature_type Character string describing the molecular features.
#'   The default is `"feature"`.
#'
#' @return A valid `LongOmicsTraj` object.
#'
#' @details
#' `lot_create()` performs the following steps:
#'
#' \enumerate{
#'   \item Checks that the expression data are numeric.
#'   \item Checks that feature and sample identifiers are available.
#'   \item Checks the required metadata columns.
#'   \item Aligns metadata rows to expression-matrix columns.
#'   \item Optionally defines visit order and numerical time.
#'   \item Creates a `SummarizedExperiment`.
#'   \item Creates a `MultiAssayExperiment`.
#'   \item Returns a valid `LongOmicsTraj` object.
#' }
#'
#' @examples
#' data("lot_glucold_expression")
#' data("lot_glucold_clinical")
#'
#' visits <- c(
#'   "Baseline",
#'   "6 Months",
#'   "30 Months"
#' )
#'
#' time_values <- c(
#'   "Baseline" = 0,
#'   "6 Months" = 6,
#'   "30 Months" = 30
#' )
#'
#' lot <- lot_create(
#'   expression = lot_glucold_expression,
#'   metadata = lot_glucold_clinical,
#'   assay = "transcriptomics",
#'   visit_levels = visits,
#'   time_values = time_values
#' )
#'
#' methods::is(
#'   lot,
#'   "LongOmicsTraj"
#' )
#'
#' names(
#'   MultiAssayExperiment::experiments(
#'     lot_mae(lot)
#'   )
#' )
#'
#' @export
lot_create <- function(
    expression,
    metadata,
    assay = "transcriptomics",
    sample_col = "sample_id",
    subject_col = "subject_id",
    visit_col = "visit",
    group_col = "group",
    visit_levels = NULL,
    time_values = NULL,
    time_col = "time",
    feature_type = "feature"
) {
  
  # ---------------------------------------------------------------------------
  # Validate basic arguments
  # ---------------------------------------------------------------------------
  
  if (missing(expression)) {
    stop(
      "`expression` must be supplied.",
      call. = FALSE
    )
  }
  
  if (missing(metadata)) {
    stop(
      "`metadata` must be supplied.",
      call. = FALSE
    )
  }
  
  if (
    !is.character(assay) ||
    length(assay) != 1L ||
    is.na(assay) ||
    !nzchar(assay)
  ) {
    stop(
      "`assay` must be one non-empty character string.",
      call. = FALSE
    )
  }
  
  if (
    !is.character(feature_type) ||
    length(feature_type) != 1L ||
    is.na(feature_type) ||
    !nzchar(feature_type)
  ) {
    stop(
      "`feature_type` must be one non-empty character string.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Prepare the expression matrix
  # ---------------------------------------------------------------------------
  
  expression <- as.matrix(expression)
  
  if (!is.numeric(expression)) {
    stop(
      "`expression` must be numeric.",
      call. = FALSE
    )
  }
  
  storage.mode(expression) <- "numeric"
  
  if (nrow(expression) == 0L) {
    stop(
      "`expression` contains no molecular features.",
      call. = FALSE
    )
  }
  
  if (ncol(expression) == 0L) {
    stop(
      "`expression` contains no samples.",
      call. = FALSE
    )
  }
  
  if (
    is.null(rownames(expression)) ||
    anyNA(rownames(expression)) ||
    any(!nzchar(rownames(expression)))
  ) {
    stop(
      "`expression` must have non-missing feature names as row names.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(rownames(expression))) {
    duplicated_features <- unique(
      rownames(expression)[duplicated(rownames(expression))]
    )
    
    stop(
      "Duplicated feature names were found in `expression`: ",
      paste(utils::head(duplicated_features, 10L), collapse = ", "),
      if (length(duplicated_features) > 10L) " ..." else "",
      call. = FALSE
    )
  }
  
  if (
    is.null(colnames(expression)) ||
    anyNA(colnames(expression)) ||
    any(!nzchar(colnames(expression)))
  ) {
    stop(
      "`expression` must have non-missing sample IDs as column names.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(colnames(expression))) {
    duplicated_samples <- unique(
      colnames(expression)[duplicated(colnames(expression))]
    )
    
    stop(
      "Duplicated sample IDs were found in `expression`: ",
      paste(utils::head(duplicated_samples, 10L), collapse = ", "),
      if (length(duplicated_samples) > 10L) " ..." else "",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(expression), na.rm = TRUE)) {
    stop(
      "`expression` contains non-finite values such as Inf or -Inf.",
      call. = FALSE
    )
  }
  
  # Missing values can be valid in some omics datasets, so report but do not stop.
  if (anyNA(expression)) {
    warning(
      "`expression` contains missing values. ",
      "Confirm that the selected estimator can handle them.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Prepare and validate metadata
  # ---------------------------------------------------------------------------
  
  metadata <- as.data.frame(
    metadata,
    stringsAsFactors = FALSE
  )
  
  if (nrow(metadata) == 0L) {
    stop(
      "`metadata` contains no samples.",
      call. = FALSE
    )
  }
  
  required_columns <- c(
    sample_col,
    subject_col,
    visit_col,
    group_col
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(metadata)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "`metadata` is missing the following required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  metadata[[sample_col]] <- as.character(
    metadata[[sample_col]]
  )
  
  metadata[[subject_col]] <- as.character(
    metadata[[subject_col]]
  )
  
  metadata[[visit_col]] <- as.character(
    metadata[[visit_col]]
  )
  
  metadata[[group_col]] <- as.character(
    metadata[[group_col]]
  )
  
  if (
    anyNA(metadata[[sample_col]]) ||
    any(!nzchar(metadata[[sample_col]]))
  ) {
    stop(
      "`metadata$", sample_col,
      "` contains missing or empty sample IDs.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(metadata[[sample_col]])) {
    duplicated_metadata_samples <- unique(
      metadata[[sample_col]][duplicated(metadata[[sample_col]])]
    )
    
    stop(
      "Duplicated sample IDs were found in `metadata$", sample_col, "`: ",
      paste(
        utils::head(duplicated_metadata_samples, 10L),
        collapse = ", "
      ),
      if (length(duplicated_metadata_samples) > 10L) " ..." else "",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Compare expression and metadata sample IDs
  # ---------------------------------------------------------------------------
  
  expression_samples <- colnames(expression)
  metadata_samples <- metadata[[sample_col]]
  
  missing_from_metadata <- setdiff(
    expression_samples,
    metadata_samples
  )
  
  missing_from_expression <- setdiff(
    metadata_samples,
    expression_samples
  )
  
  if (
    length(missing_from_metadata) > 0L ||
    length(missing_from_expression) > 0L
  ) {
    
    error_lines <- c(
      "Expression and metadata sample IDs do not match.",
      paste0(
        "Expression matrix: ",
        length(expression_samples),
        " samples."
      ),
      paste0(
        "Metadata: ",
        length(metadata_samples),
        " rows."
      )
    )
    
    if (length(missing_from_metadata) > 0L) {
      error_lines <- c(
        error_lines,
        paste0(
          "Missing from metadata: ",
          paste(
            utils::head(missing_from_metadata, 10L),
            collapse = ", "
          ),
          if (length(missing_from_metadata) > 10L) " ..." else ""
        )
      )
    }
    
    if (length(missing_from_expression) > 0L) {
      error_lines <- c(
        error_lines,
        paste0(
          "Missing from expression matrix: ",
          paste(
            utils::head(missing_from_expression, 10L),
            collapse = ", "
          ),
          if (length(missing_from_expression) > 10L) " ..." else ""
        )
      )
    }
    
    stop(
      paste(error_lines, collapse = "\n"),
      call. = FALSE
    )
  }
  
  # Reorder metadata to match the expression matrix exactly.
  metadata <- metadata[
    match(
      expression_samples,
      metadata[[sample_col]]
    ),
    ,
    drop = FALSE
  ]
  
  rownames(metadata) <- metadata[[sample_col]]
  
  stopifnot(
    identical(
      colnames(expression),
      metadata[[sample_col]]
    )
  )
  
  # ---------------------------------------------------------------------------
  # Define visit order
  # ---------------------------------------------------------------------------
  
  if (!is.null(visit_levels)) {
    
    if (
      !is.character(visit_levels) ||
      length(visit_levels) == 0L ||
      anyNA(visit_levels) ||
      any(!nzchar(visit_levels))
    ) {
      stop(
        "`visit_levels` must be a non-empty character vector.",
        call. = FALSE
      )
    }
    
    observed_visits <- unique(
      metadata[[visit_col]]
    )
    
    unknown_visits <- setdiff(
      observed_visits,
      visit_levels
    )
    
    if (length(unknown_visits) > 0L) {
      stop(
        "The following observed visits are absent from `visit_levels`: ",
        paste(unknown_visits, collapse = ", "),
        call. = FALSE
      )
    }
    
    metadata[[visit_col]] <- factor(
      metadata[[visit_col]],
      levels = visit_levels,
      ordered = TRUE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Create numerical time
  # ---------------------------------------------------------------------------
  
  if (!is.null(time_values)) {
    
    if (
      !is.numeric(time_values) ||
      is.null(names(time_values)) ||
      anyNA(names(time_values)) ||
      any(!nzchar(names(time_values)))
    ) {
      stop(
        "`time_values` must be a named numeric vector.",
        call. = FALSE
      )
    }
    
    visit_labels <- as.character(
      metadata[[visit_col]]
    )
    
    missing_time_values <- setdiff(
      unique(visit_labels),
      names(time_values)
    )
    
    if (length(missing_time_values) > 0L) {
      stop(
        "`time_values` does not define numerical time for: ",
        paste(missing_time_values, collapse = ", "),
        call. = FALSE
      )
    }
    
    metadata[[time_col]] <- unname(
      time_values[visit_labels]
    )
    
    if (anyNA(metadata[[time_col]])) {
      stop(
        "Numerical time could not be assigned to every sample.",
        call. = FALSE
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Create feature metadata
  # ---------------------------------------------------------------------------
  
  feature_data <- S4Vectors::DataFrame(
    object_id = rownames(expression),
    object_type = rep(
      feature_type,
      nrow(expression)
    ),
    row.names = rownames(expression)
  )
  
  # ---------------------------------------------------------------------------
  # Create the SummarizedExperiment
  # ---------------------------------------------------------------------------
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = stats::setNames(
      list(expression),
      assay
    ),
    rowData = feature_data,
    colData = S4Vectors::DataFrame(
      metadata,
      row.names = metadata[[sample_col]]
    )
  )
  
  # ---------------------------------------------------------------------------
  # Create the MultiAssayExperiment
  # ---------------------------------------------------------------------------
  
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments =
      MultiAssayExperiment::ExperimentList(
        stats::setNames(
          list(se),
          assay
        )
      )
  )
  
  # ---------------------------------------------------------------------------
  # Create and validate the LongOmicsTraj object
  # ---------------------------------------------------------------------------
  
  object <- methods::new(
    "LongOmicsTraj",
    mae = mae
  )
  
  methods::validObject(object)
  
  object
}