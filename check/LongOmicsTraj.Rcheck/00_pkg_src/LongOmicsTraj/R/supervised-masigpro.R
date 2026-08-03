#' Estimate trajectories using a maSigPro-style polynomial model
#'
#' Fits polynomial regression models to longitudinal molecular features and
#' stores predicted group-by-visit trajectories in the visit mean matrix.
#'
#' @param object A `LongOmicsTraj` object.
#' @param assay Name of the assay to analyse.
#' @param group_col Metadata column defining clinical or treatment groups.
#' @param visit_col Metadata column defining study visits.
#' @param time_col Metadata column containing numeric visit times.
#' @param subject_col Metadata column containing subject identifiers.
#'   Currently retained for interface consistency.
#' @param visits Character vector giving visits in chronological order.
#' @param degree Optional polynomial degree. When `NULL`, the degree is set to
#'   the number of unique time points minus one, with a minimum of one.
#' @param estimator Name stored for the trajectory estimator.
#' @param overwrite Logical; overwrite an existing visit mean matrix.
#' @param verbose Logical; print progress messages.
#'
#' @return A `LongOmicsTraj` object with updated visit mean trajectories.
#'
#' @rdname lot_trajectory_estimators
#' @export
methods::setMethod(
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
    
    experiments <- MultiAssayExperiment::experiments(
      object@mae
    )
    
    if (!assay %in% names(experiments)) {
      stop(
        "Assay not found: ",
        assay,
        call. = FALSE
      )
    }
    
    if (
      nrow(object@vmm) > 0L &&
      !isTRUE(overwrite)
    ) {
      stop(
        "VMM exists. Use `overwrite = TRUE`.",
        call. = FALSE
      )
    }
    
    se <- experiments[[assay]]
    
    mat <- as.matrix(
      SummarizedExperiment::assay(
        se
      )
    )
    
    meta <- as.data.frame(
      SummarizedExperiment::colData(
        se
      )
    )
    
    required_metadata <- c(
      group_col,
      visit_col,
      time_col,
      subject_col
    )
    
    missing_metadata <- setdiff(
      required_metadata,
      colnames(meta)
    )
    
    if (length(missing_metadata) > 0L) {
      stop(
        "Metadata are missing: ",
        paste(
          missing_metadata,
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    
    if (is.null(rownames(meta))) {
      stop(
        "Metadata must have sample row names.",
        call. = FALSE
      )
    }
    
    if (!all(rownames(meta) %in% colnames(mat))) {
      stop(
        "Metadata samples could not be matched to assay columns.",
        call. = FALSE
      )
    }
    
    mat <- mat[
      ,
      rownames(meta),
      drop = FALSE
    ]
    
    time_values <- as.numeric(
      meta[[time_col]]
    )
    
    group_values <- factor(
      meta[[group_col]]
    )
    
    if (anyNA(time_values)) {
      stop(
        "The time column must contain numeric, non-missing values.",
        call. = FALSE
      )
    }
    
    if (
      missing(visits) ||
      length(visits) < 2L ||
      anyNA(visits)
    ) {
      stop(
        "`visits` must contain at least two ordered visit labels.",
        call. = FALSE
      )
    }
    
    if (is.null(degree)) {
      n_unique_times <- length(
        unique(
          time_values
        )
      )
      
      degree <- max(
        1L,
        n_unique_times - 1L
      )
    }
    
    if (
      length(degree) != 1L ||
      !is.numeric(degree) ||
      is.na(degree) ||
      !is.finite(degree) ||
      degree < 1L ||
      degree != as.integer(degree)
    ) {
      stop(
        "`degree` must be one positive integer.",
        call. = FALSE
      )
    }
    
    degree <- as.integer(
      degree
    )
    
    n_groups <- nlevels(
      group_values
    )
    
    if (isTRUE(verbose)) {
      message(
        "Fitting polynomial trajectories via multivariate regression..."
      )
    }
    
    model_data <- data.frame(
      Time = time_values,
      Group = group_values
    )
    
    if (n_groups > 1L) {
      fit <- stats::lm(
        t(mat) ~ stats::poly(
          Time,
          degree,
          raw = TRUE
        ) * Group,
        data = model_data
      )
    } else {
      fit <- stats::lm(
        t(mat) ~ stats::poly(
          Time,
          degree,
          raw = TRUE
        ),
        data = model_data
      )
    }
    
    coefficients <- t(
      stats::coef(
        fit
      )
    )
    
    coefficients[
      is.na(coefficients)
    ] <- 0
    
    feature_ids <- rownames(
      coefficients
    )
    
    prediction_grid <- expand.grid(
      group = unique(
        as.character(
          meta[[group_col]]
        )
      ),
      visit = visits,
      stringsAsFactors = FALSE
    )
    
    visit_map <- meta[
      ,
      c(
        visit_col,
        time_col
      ),
      drop = FALSE
    ]
    
    visit_map <- unique(
      visit_map
    )
    
    names(visit_map) <- c(
      "visit",
      "time_value"
    )
    
    prediction_grid <- merge(
      prediction_grid,
      visit_map,
      by = "visit",
      all.x = TRUE,
      sort = FALSE
    )
    
    visit_order <- match(
      prediction_grid$visit,
      visits
    )
    
    prediction_grid <- prediction_grid[
      order(
        prediction_grid$group,
        visit_order
      ),
      ,
      drop = FALSE
    ]
    
    if (anyNA(prediction_grid$time_value)) {
      missing_visits <- unique(
        prediction_grid$visit[
          is.na(
            prediction_grid$time_value
          )
        ]
      )
      
      stop(
        "No numeric time mapping was found for visits: ",
        paste(
          missing_visits,
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    
    prediction_data <- data.frame(
      Time = as.numeric(
        prediction_grid$time_value
      ),
      Group = factor(
        prediction_grid$group,
        levels = levels(
          group_values
        )
      )
    )
    
    if (n_groups > 1L) {
      prediction_matrix <- stats::model.matrix(
        ~ stats::poly(
          Time,
          degree,
          raw = TRUE
        ) * Group,
        data = prediction_data
      )
    } else {
      prediction_matrix <- stats::model.matrix(
        ~ stats::poly(
          Time,
          degree,
          raw = TRUE
        ),
        data = prediction_data
      )
    }
    
    predicted_values <- coefficients %*%
      t(
        prediction_matrix
      )
    
    vmm_list <- lapply(
      seq_len(
        nrow(predicted_values)
      ),
      function(i) {
        data.frame(
          object_id = feature_ids[[i]],
          object_type = "gene",
          assay = assay,
          group = prediction_grid$group,
          visit = prediction_grid$visit,
          estimated_value = as.numeric(
            predicted_values[
              i,
              ,
              drop = TRUE
            ]
          ),
          estimator = estimator,
          stringsAsFactors = FALSE
        )
      }
    )
    
    vmm_new <- S4Vectors::DataFrame(
      do.call(
        rbind,
        vmm_list
      )
    )
    
    object@vmm <- vmm_new
    
    object@history <- c(
      object@history,
      list(
        list(
          step = "lot_from_masigpro",
          assay = assay,
          estimator = estimator,
          degree = degree,
          n_features = nrow(mat),
          time = Sys.time()
        )
      )
    )
    
    if (isTRUE(verbose)) {
      message(
        "Created maSigPro-equivalent VMM for assay '",
        assay,
        "' with ",
        nrow(vmm_new),
        " rows."
      )
    }
    
    methods::validObject(
      object
    )
    
    object
  }
)