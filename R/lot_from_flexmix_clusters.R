#' Cluster molecular features using FlexMix trajectory models
#'
#' Fits finite mixtures of longitudinal regression trajectories to molecular
#' features within one clinical or treatment group. Candidate cluster
#' resolutions are evaluated using likelihood-based model statistics and
#' LongOmicsTraj topology purity.
#'
#' Trajectories can be obtained from model-estimated values already stored in
#' `object@vmm`, or calculated as mean molecular values at each visit from the
#' original assay matrix.
#'
#' @param object A `LongOmicsTraj` object containing molecular data and
#'   feature-level OTS results.
#' @param assay Name of the experiment in the `MultiAssayExperiment`.
#' @param group Clinical or treatment group to analyse.
#' @param visits Character vector giving the chronological visit order.
#' @param k Integer vector containing candidate numbers of clusters.
#' @param assay_data Optional assay-matrix name. The first available matrix is
#'   used by default.
#' @param group_col Metadata column containing clinical groups.
#' @param visit_col Metadata column containing visits.
#' @param sample_col Metadata column containing sample identifiers.
#' @param trajectory_source Source of the gene trajectories: `"estimated"`
#'   uses `object@vmm`, while `"mean"` calculates visit-level assay means.
#' @param trajectory_estimator Optional estimator used to select trajectories
#'   from `object@vmm`, such as `"masigpro"`, `"gam"`, `"lmm"` or
#'   `"empirical"`.
#' @param model Trajectory regression model: `"quadratic"` or `"linear"`.
#' @param scale_features Logical; standardise each molecular-feature trajectory
#'   before clustering.
#' @param purity_threshold Minimum dominant-OTS proportion required for a
#'   cluster to pass the topology-purity criterion.
#' @param required_passing_fraction Minimum fraction of clusters required to
#'   pass the topology-purity criterion.
#' @param min_cluster_size Minimum allowed number of molecular features in a
#'   cluster.
#' @param nrep Number of repeated FlexMix initialisations for each candidate
#'   value of `k`.
#' @param seed Random seed.
#' @param store_fits Logical; store the fitted FlexMix models in
#'   `object@results`.
#' @param overwrite Logical; overwrite existing FlexMix clustering results.
#' @param verbose Logical; print progress messages.
#'
#' @return A `LongOmicsTraj` object containing model-comparison statistics,
#'   the selected cluster resolution, feature assignments, cluster
#'   trajectories, OTS distributions and topology-purity summaries.
#'
#' @details
#' The selected solution is the smallest candidate `k` for which:
#'
#' \itemize{
#'   \item the required fraction of clusters reaches `purity_threshold`; and
#'   \item all clusters contain at least `min_cluster_size` features.
#' }
#'
#' If no candidate satisfies both criteria, the solution with the lowest BIC
#' is selected.
#'
#' When `trajectory_source = "estimated"`, LongOmicsTraj trajectories stored
#' in `object@vmm` are used. This is generally preferred because trajectories
#' may already account for repeated measurements, smoothing or other model
#' assumptions.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("lot_glucold_ots")
#'
#' lot_flexmix <- lot_from_flexmix_clusters(
#'   object = lot_glucold_ots,
#'   assay = "transcriptomics",
#'   group = "ICS 6 months then withdrawal",
#'   visits = c(
#'     "Baseline",
#'     "6 Months",
#'     "30 Months"
#'   ),
#'   trajectory_source = "estimated",
#'   trajectory_estimator = "masigpro",
#'   k = 2:5,
#'   purity_threshold = 0.90,
#'   required_passing_fraction = 0.80,
#'   min_cluster_size = 5
#' )
#' }
lot_from_flexmix_clusters <- function(
    object,
    assay,
    group,
    visits,
    k = 2:5,
    assay_data = NULL,
    group_col = "group",
    visit_col = "visit",
    sample_col = "sample_id",
    trajectory_source = c("estimated", "mean"),
    trajectory_estimator = NULL,
    model = c("quadratic", "linear"),
    scale_features = TRUE,
    purity_threshold = 0.90,
    required_passing_fraction = 0.80,
    min_cluster_size = 5L,
    nrep = 10L,
    seed = 42L,
    store_fits = FALSE,
    overwrite = FALSE,
    verbose = TRUE
) {
  
  trajectory_source <- match.arg(
    trajectory_source
  )
  
  model <- match.arg(
    model
  )
  
  # ---------------------------------------------------------------------------
  # Validate dependencies and object
  # ---------------------------------------------------------------------------
  
  if (!requireNamespace("flexmix", quietly = TRUE)) {
    stop(
      "Package `flexmix` is required. Install it with ",
      "`install.packages(\"flexmix\")`.",
      call. = FALSE
    )
  }
  
  if (!methods::is(object, "LongOmicsTraj")) {
    stop(
      "`object` must be a LongOmicsTraj object.",
      call. = FALSE
    )
  }
  
  if (
    missing(assay) ||
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
    missing(group) ||
    !is.character(group) ||
    length(group) != 1L ||
    is.na(group) ||
    !nzchar(group)
  ) {
    stop(
      "`group` must be one non-empty character string.",
      call. = FALSE
    )
  }
  
  if (
    missing(visits) ||
    !is.character(visits) ||
    length(visits) < 2L ||
    anyNA(visits) ||
    any(!nzchar(visits)) ||
    anyDuplicated(visits)
  ) {
    stop(
      "`visits` must contain at least two unique, non-missing labels.",
      call. = FALSE
    )
  }
  
  if (
    !is.numeric(k) ||
    length(k) == 0L ||
    anyNA(k) ||
    any(!is.finite(k)) ||
    any(k < 2L) ||
    any(k != as.integer(k))
  ) {
    stop(
      "`k` must contain integers greater than or equal to 2.",
      call. = FALSE
    )
  }
  
  k <- sort(
    unique(
      as.integer(k)
    )
  )
  
  if (
    length(purity_threshold) != 1L ||
    !is.numeric(purity_threshold) ||
    is.na(purity_threshold) ||
    !is.finite(purity_threshold) ||
    purity_threshold < 0 ||
    purity_threshold > 1
  ) {
    stop(
      "`purity_threshold` must be between 0 and 1.",
      call. = FALSE
    )
  }
  
  if (
    length(required_passing_fraction) != 1L ||
    !is.numeric(required_passing_fraction) ||
    is.na(required_passing_fraction) ||
    !is.finite(required_passing_fraction) ||
    required_passing_fraction < 0 ||
    required_passing_fraction > 1
  ) {
    stop(
      "`required_passing_fraction` must be between 0 and 1.",
      call. = FALSE
    )
  }
  
  if (
    length(min_cluster_size) != 1L ||
    !is.numeric(min_cluster_size) ||
    is.na(min_cluster_size) ||
    min_cluster_size < 1L ||
    min_cluster_size != as.integer(min_cluster_size)
  ) {
    stop(
      "`min_cluster_size` must be one positive integer.",
      call. = FALSE
    )
  }
  
  if (
    length(nrep) != 1L ||
    !is.numeric(nrep) ||
    is.na(nrep) ||
    nrep < 1L ||
    nrep != as.integer(nrep)
  ) {
    stop(
      "`nrep` must be one positive integer.",
      call. = FALSE
    )
  }
  
  if (
    length(seed) != 1L ||
    !is.numeric(seed) ||
    is.na(seed) ||
    !is.finite(seed)
  ) {
    stop(
      "`seed` must be one finite numeric value.",
      call. = FALSE
    )
  }
  
  logical_arguments <- list(
    scale_features = scale_features,
    store_fits = store_fits,
    overwrite = overwrite,
    verbose = verbose
  )
  
  for (argument_name in names(logical_arguments)) {
    argument_value <- logical_arguments[[argument_name]]
    
    if (
      !is.logical(argument_value) ||
      length(argument_value) != 1L ||
      is.na(argument_value)
    ) {
      stop(
        "`",
        argument_name,
        "` must be `TRUE` or `FALSE`.",
        call. = FALSE
      )
    }
  }
  
  min_cluster_size <- as.integer(
    min_cluster_size
  )
  
  nrep <- as.integer(
    nrep
  )
  
  seed <- as.integer(
    seed
  )
  
  # ---------------------------------------------------------------------------
  # Protect existing FlexMix results
  # ---------------------------------------------------------------------------
  
  flexmix_result_names <- c(
    "flexmix_model_comparison",
    "flexmix_selected_k",
    "flexmix_selection_reason",
    "flexmix_gene_assignments",
    "flexmix_cluster_trajectories",
    "flexmix_ots_distribution",
    "flexmix_cluster_purity",
    "flexmix_overall_purity",
    "flexmix_trajectory_data",
    "flexmix_fits"
  )
  
  existing_flexmix_results <- intersect(
    flexmix_result_names,
    names(object@results)
  )
  
  if (
    length(existing_flexmix_results) > 0L &&
    !isTRUE(overwrite)
  ) {
    stop(
      "FlexMix clustering results already exist. ",
      "Use `overwrite = TRUE` to replace them.",
      call. = FALSE
    )
  }
  
  if (isTRUE(overwrite)) {
    for (result_name in existing_flexmix_results) {
      object@results[[result_name]] <- NULL
    }
  }
  
  # ---------------------------------------------------------------------------
  # Extract experiment, assay matrix and metadata
  # ---------------------------------------------------------------------------
  
  experiments <- MultiAssayExperiment::experiments(
    object@mae
  )
  
  experiment_names <- names(
    experiments
  )
  
  if (!assay %in% experiment_names) {
    stop(
      "Assay `",
      assay,
      "` was not found. Available assays: ",
      paste(
        experiment_names,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  se <- experiments[[assay]]
  
  available_matrices <- SummarizedExperiment::assayNames(
    se
  )
  
  if (length(available_matrices) == 0L) {
    stop(
      "The selected experiment contains no assay matrices.",
      call. = FALSE
    )
  }
  
  if (is.null(assay_data)) {
    assay_data <- available_matrices[[1L]]
  }
  
  if (
    !is.character(assay_data) ||
    length(assay_data) != 1L ||
    is.na(assay_data) ||
    !nzchar(assay_data) ||
    !assay_data %in% available_matrices
  ) {
    stop(
      "Assay matrix `",
      assay_data,
      "` was not found. Available matrices: ",
      paste(
        available_matrices,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  expression_matrix <- SummarizedExperiment::assay(
    se,
    assay_data
  )
  
  if (!is.numeric(expression_matrix)) {
    stop(
      "The selected assay matrix must be numeric.",
      call. = FALSE
    )
  }
  
  if (
    is.null(rownames(expression_matrix)) ||
    anyNA(rownames(expression_matrix)) ||
    any(!nzchar(rownames(expression_matrix))) ||
    anyDuplicated(rownames(expression_matrix))
  ) {
    stop(
      "The assay matrix must have unique, non-missing feature identifiers.",
      call. = FALSE
    )
  }
  
  metadata <- as.data.frame(
    SummarizedExperiment::colData(se)
  )
  
  required_metadata <- unique(
    c(
      sample_col,
      group_col,
      visit_col
    )
  )
  
  missing_metadata <- setdiff(
    required_metadata,
    colnames(metadata)
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
  
  metadata[[sample_col]] <- as.character(
    metadata[[sample_col]]
  )
  
  if (
    anyNA(metadata[[sample_col]]) ||
    any(!nzchar(metadata[[sample_col]])) ||
    anyDuplicated(metadata[[sample_col]])
  ) {
    stop(
      "The metadata sample identifiers must be unique and non-missing.",
      call. = FALSE
    )
  }
  
  metadata <- metadata[
    match(
      colnames(expression_matrix),
      metadata[[sample_col]]
    ),
    ,
    drop = FALSE
  ]
  
  if (
    anyNA(metadata[[sample_col]]) ||
    !identical(
      colnames(expression_matrix),
      metadata[[sample_col]]
    )
  ) {
    stop(
      "Expression samples could not be aligned with the metadata.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Select the requested clinical group and visits
  # ---------------------------------------------------------------------------
  
  group_values <- as.character(
    metadata[[group_col]]
  )
  
  visit_values <- as.character(
    metadata[[visit_col]]
  )
  
  keep_samples <- group_values == group &
    visit_values %in% visits
  
  metadata_group <- metadata[
    keep_samples,
    ,
    drop = FALSE
  ]
  
  if (nrow(metadata_group) == 0L) {
    stop(
      "No samples were found for group `",
      group,
      "` at the requested visits.",
      call. = FALSE
    )
  }
  
  metadata_group[[visit_col]] <- factor(
    as.character(
      metadata_group[[visit_col]]
    ),
    levels = visits,
    ordered = TRUE
  )
  
  visit_counts <- table(
    metadata_group[[visit_col]]
  )
  
  missing_visits <- names(
    visit_counts[
      visit_counts == 0L
    ]
  )
  
  if (length(missing_visits) > 0L) {
    stop(
      "No samples were available for visits: ",
      paste(
        missing_visits,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  visit_time <- stats::setNames(
    seq_along(visits) - 1,
    visits
  )
  
  # ---------------------------------------------------------------------------
  # Construct one longitudinal trajectory per feature
  # ---------------------------------------------------------------------------
  
  if (trajectory_source == "estimated") {
    
    # -------------------------------------------------------------------------
    # Use model-estimated trajectories stored in object@vmm
    # -------------------------------------------------------------------------
    
    vmm_data <- as.data.frame(
      object@vmm
    )
    
    if (nrow(vmm_data) == 0L) {
      stop(
        "No estimated trajectories were found in `object@vmm`. ",
        "Run a trajectory estimator first or use ",
        "`trajectory_source = \"mean\"`.",
        call. = FALSE
      )
    }
    
    required_vmm_columns <- c(
      "object_id",
      "assay",
      "group",
      "visit",
      "estimated_value"
    )
    
    missing_vmm_columns <- setdiff(
      required_vmm_columns,
      colnames(vmm_data)
    )
    
    if (length(missing_vmm_columns) > 0L) {
      stop(
        "The estimated trajectories are missing: ",
        paste(
          missing_vmm_columns,
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    
    keep_vmm <- as.character(vmm_data$assay) == assay &
      as.character(vmm_data$group) == group &
      as.character(vmm_data$visit) %in% visits
    
    vmm_data <- vmm_data[
      keep_vmm,
      ,
      drop = FALSE
    ]
    
    if (nrow(vmm_data) == 0L) {
      stop(
        "No estimated trajectories were found for assay `",
        assay,
        "`, group `",
        group,
        "` and the requested visits.",
        call. = FALSE
      )
    }
    
    if ("estimator" %in% colnames(vmm_data)) {
      
      available_estimators <- sort(
        unique(
          as.character(
            vmm_data$estimator[
              !is.na(vmm_data$estimator)
            ]
          )
        )
      )
      
      available_estimators <- available_estimators[
        nzchar(available_estimators)
      ]
      
      if (!is.null(trajectory_estimator)) {
        
        if (
          !is.character(trajectory_estimator) ||
          length(trajectory_estimator) != 1L ||
          is.na(trajectory_estimator) ||
          !nzchar(trajectory_estimator)
        ) {
          stop(
            "`trajectory_estimator` must be one non-empty character string.",
            call. = FALSE
          )
        }
        
        if (!trajectory_estimator %in% available_estimators) {
          stop(
            "Trajectory estimator `",
            trajectory_estimator,
            "` was not found. Available estimators: ",
            paste(
              available_estimators,
              collapse = ", "
            ),
            call. = FALSE
          )
        }
        
        vmm_data <- vmm_data[
          as.character(vmm_data$estimator) ==
            trajectory_estimator,
          ,
          drop = FALSE
        ]
        
      } else if (length(available_estimators) > 1L) {
        
        stop(
          "Multiple trajectory estimators are available: ",
          paste(
            available_estimators,
            collapse = ", "
          ),
          ". Select one using `trajectory_estimator`.",
          call. = FALSE
        )
        
      } else if (length(available_estimators) == 1L) {
        
        trajectory_estimator <- available_estimators[[1L]]
      }
      
    } else if (!is.null(trajectory_estimator)) {
      
      stop(
        "`object@vmm` does not contain an `estimator` column.",
        call. = FALSE
      )
    }
    
    trajectory_data <- data.frame(
      gene_id = as.character(
        vmm_data$object_id
      ),
      visit = as.character(
        vmm_data$visit
      ),
      mean_expression = as.numeric(
        vmm_data$estimated_value
      ),
      stringsAsFactors = FALSE
    )
    
    trajectory_data <- trajectory_data |>
      dplyr::group_by(
        .data$gene_id,
        .data$visit
      ) |>
      dplyr::summarise(
        mean_expression = mean(
          .data$mean_expression,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
    
    trajectory_data$time <- unname(
      visit_time[
        trajectory_data$visit
      ]
    )
    
    trajectory_data$n_samples <- NA_integer_
    
  } else {
    
    # -------------------------------------------------------------------------
    # Calculate mean assay value for each feature and visit
    # -------------------------------------------------------------------------
    
    trajectory_list <- lapply(
      visits,
      function(current_visit) {
        
        sample_ids <- metadata_group[[sample_col]][
          as.character(
            metadata_group[[visit_col]]
          ) == current_visit
        ]
        
        data.frame(
          gene_id = rownames(
            expression_matrix
          ),
          visit = current_visit,
          time = unname(
            visit_time[[current_visit]]
          ),
          mean_expression = rowMeans(
            expression_matrix[
              ,
              sample_ids,
              drop = FALSE
            ],
            na.rm = TRUE
          ),
          n_samples = length(
            sample_ids
          ),
          stringsAsFactors = FALSE
        )
      }
    )
    
    trajectory_data <- dplyr::bind_rows(
      trajectory_list
    )
  }
  
  trajectory_data$visit <- factor(
    as.character(
      trajectory_data$visit
    ),
    levels = visits,
    ordered = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # Remove incomplete and invariant feature trajectories
  # ---------------------------------------------------------------------------
  
  valid_features <- trajectory_data |>
    dplyr::group_by(
      .data$gene_id
    ) |>
    dplyr::summarise(
      n_observed_visits = sum(
        is.finite(.data$mean_expression)
      ),
      trajectory_sd = stats::sd(
        .data$mean_expression,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::filter(
      .data$n_observed_visits == length(.env$visits),
      is.finite(.data$trajectory_sd),
      .data$trajectory_sd > 0
    ) |>
    dplyr::pull(
      .data$gene_id
    )
  
  trajectory_data <- trajectory_data |>
    dplyr::filter(
      .data$gene_id %in% .env$valid_features
    )
  
  n_features <- dplyr::n_distinct(
    trajectory_data$gene_id
  )
  
  if (n_features == 0L) {
    stop(
      "No complete, variable feature trajectories were available.",
      call. = FALSE
    )
  }
  
  if (n_features < max(k)) {
    stop(
      "The largest requested `k` exceeds the number of retained features (",
      n_features,
      ").",
      call. = FALSE
    )
  }
  
  if (
    max(k) * min_cluster_size >
    n_features
  ) {
    warning(
      "For some candidate values of `k`, it is mathematically impossible ",
      "for every cluster to contain at least `min_cluster_size` features.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Standardise feature trajectories
  # ---------------------------------------------------------------------------
  
  if (isTRUE(scale_features)) {
    
    trajectory_data <- trajectory_data |>
      dplyr::group_by(
        .data$gene_id
      ) |>
      dplyr::mutate(
        clustering_value = as.numeric(
          base::scale(
            .data$mean_expression
          )
        )
      ) |>
      dplyr::ungroup()
    
  } else {
    
    trajectory_data$clustering_value <-
      trajectory_data$mean_expression
  }
  
  if (any(!is.finite(trajectory_data$clustering_value))) {
    stop(
      "Non-finite clustering values were produced.",
      call. = FALSE
    )
  }
  
  unique_times <- sort(
    unique(
      trajectory_data$time
    )
  )
  
  time_mean <- mean(
    unique_times
  )
  
  time_sd <- stats::sd(
    unique_times
  )
  
  if (
    !is.finite(time_sd) ||
    time_sd <= 0
  ) {
    stop(
      "Visit times could not be standardised.",
      call. = FALSE
    )
  }
  
  trajectory_data$time_scaled <-
    (
      trajectory_data$time -
        time_mean
    ) / time_sd
  
  # ---------------------------------------------------------------------------
  # Join existing feature-level OTS labels
  # ---------------------------------------------------------------------------
  
  ots_data <- as.data.frame(
    object@ots
  )
  
  if (nrow(ots_data) == 0L) {
    stop(
      "Feature-level OTS results are required to calculate topology purity.",
      call. = FALSE
    )
  }
  
  required_ots_columns <- c(
    "object_id",
    "group",
    "topology_label"
  )
  
  missing_ots_columns <- setdiff(
    required_ots_columns,
    colnames(ots_data)
  )
  
  if (length(missing_ots_columns) > 0L) {
    stop(
      "Feature-level OTS results are missing: ",
      paste(
        missing_ots_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  ots_keep <- as.character(
    ots_data$group
  ) == group
  
  if (
    !is.null(trajectory_estimator) &&
    "estimator" %in% colnames(ots_data)
  ) {
    ots_keep <- ots_keep &
      as.character(ots_data$estimator) ==
      trajectory_estimator
  }
  
  feature_ots <- ots_data[
    ots_keep,
    ,
    drop = FALSE
  ]
  
  feature_ots <- feature_ots |>
    dplyr::transmute(
      gene_id = as.character(
        .data$object_id
      ),
      topology_label = as.character(
        .data$topology_label
      )
    ) |>
    dplyr::distinct()
  
  duplicated_ots <- feature_ots |>
    dplyr::count(
      .data$gene_id,
      name = "n_labels"
    ) |>
    dplyr::filter(
      .data$n_labels > 1L
    )
  
  if (nrow(duplicated_ots) > 0L) {
    stop(
      "Some features have more than one OTS label for the selected group",
      if (!is.null(trajectory_estimator)) {
        paste0(
          " and estimator `",
          trajectory_estimator,
          "`"
        )
      } else {
        ""
      },
      ".",
      call. = FALSE
    )
  }
  
  trajectory_data$gene_id <- as.character(
    trajectory_data$gene_id
  )
  
  trajectory_data <- trajectory_data |>
    dplyr::left_join(
      feature_ots,
      by = "gene_id"
    )
  
  missing_ots_count <- dplyr::n_distinct(
    trajectory_data$gene_id[
      is.na(trajectory_data$topology_label)
    ]
  )
  
  if (
    missing_ots_count > 0L &&
    isTRUE(verbose)
  ) {
    message(
      missing_ots_count,
      " retained feature",
      if (missing_ots_count == 1L) "" else "s",
      " had no matching OTS label and will be excluded from purity ",
      "calculations."
    )
  }
  
  trajectory_data$gene_id <- factor(
    trajectory_data$gene_id
  )
  
  # ---------------------------------------------------------------------------
  # Internal function: fit one candidate resolution
  # ---------------------------------------------------------------------------
  
  fit_one_k <- function(k_value) {
    
    if (isTRUE(verbose)) {
      message(
        "Fitting FlexMix trajectory model for k = ",
        k_value,
        "."
      )
    }
    
    if (
      model == "quadratic" &&
      length(visits) < 3L
    ) {
      stop(
        "A quadratic model requires at least three visits.",
        call. = FALSE
      )
    }
    
    model_driver <- if (model == "quadratic") {
      
      flexmix::FLXMRglm(
        clustering_value ~
          time_scaled +
          I(time_scaled^2),
        family = "gaussian"
      )
      
    } else {
      
      flexmix::FLXMRglm(
        clustering_value ~
          time_scaled,
        family = "gaussian"
      )
    }
    

    fit <- withr::with_seed(
      seed + k_value,
      flexmix::stepFlexmix(
        . ~ . | gene_id,
        data = trajectory_data,
        k = k_value,
        nrep = nrep,
        model = model_driver,
        control = list(
          iter.max = 500,
          tolerance = 1e-8,
          verbose = 0
        )
      )
    )
    
    if (!methods::is(fit, "flexmix")) {
      fit <- flexmix::getModel(
        fit,
        "BIC"
      )
    }
    
    if (is.null(fit)) {
      stop(
        "FlexMix did not return a valid model for k = ",
        k_value,
        ".",
        call. = FALSE
      )
    }
    
    cluster_vector <- flexmix::clusters(
      fit
    )
    
    if (length(cluster_vector) != nrow(trajectory_data)) {
      stop(
        "FlexMix returned ",
        length(cluster_vector),
        " assignments for ",
        nrow(trajectory_data),
        " trajectory observations.",
        call. = FALSE
      )
    }
    
    plot_data <- trajectory_data
    
    plot_data$cluster <- factor(
      cluster_vector
    )
    
    cluster_check <- plot_data |>
      dplyr::group_by(
        .data$gene_id
      ) |>
      dplyr::summarise(
        n_clusters = dplyr::n_distinct(
          .data$cluster
        ),
        .groups = "drop"
      )
    
    if (any(cluster_check$n_clusters != 1L)) {
      stop(
        "One or more features were assigned to multiple clusters.",
        call. = FALSE
      )
    }
    
    assignments <- plot_data |>
      dplyr::distinct(
        .data$gene_id,
        .data$cluster,
        .data$topology_label
      ) |>
      dplyr::mutate(
        group = .env$group,
        k = .env$k_value
      )
    
    purity_assignments <- assignments |>
      dplyr::filter(
        !is.na(.data$topology_label)
      )
    
    if (nrow(purity_assignments) == 0L) {
      stop(
        "No matching OTS labels were available for purity calculation.",
        call. = FALSE
      )
    }
    
    purity_result <- lot_cluster_purity(
      assignments = purity_assignments,
      purity_threshold = purity_threshold
    )
    
    cluster_trajectories <- plot_data |>
      dplyr::group_by(
        .data$cluster,
        .data$visit,
        .data$time
      ) |>
      dplyr::summarise(
        mean_clustering_value = mean(
          .data$clustering_value,
          na.rm = TRUE
        ),
        mean_trajectory_value = mean(
          .data$mean_expression,
          na.rm = TRUE
        ),
        n_genes = dplyr::n_distinct(
          .data$gene_id
        ),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        group = .env$group,
        k = .env$k_value
      )
    
    ots_distribution <- purity_assignments |>
      dplyr::count(
        .data$cluster,
        .data$group,
        .data$topology_label,
        name = "n_genes"
      ) |>
      dplyr::group_by(
        .data$cluster,
        .data$group
      ) |>
      dplyr::mutate(
        proportion = .data$n_genes /
          sum(.data$n_genes)
      ) |>
      dplyr::ungroup()
    
    cluster_sizes <- assignments |>
      dplyr::distinct(
        .data$gene_id,
        .data$cluster
      ) |>
      dplyr::count(
        .data$cluster,
        name = "n_genes"
      )
    
    minimum_observed_cluster_size <- min(
      cluster_sizes$n_genes
    )
    
    log_likelihood_object <- stats4::logLik(
      fit
    )
    
    log_likelihood <- as.numeric(
      log_likelihood_object
    )
    
    n_parameters <- attr(
      log_likelihood_object,
      "df"
    )
    
    if (
      is.null(n_parameters) ||
      length(n_parameters) != 1L ||
      is.na(n_parameters) ||
      !is.finite(n_parameters)
    ) {
      stop(
        "Could not determine the number of fitted parameters for k = ",
        k_value,
        ".",
        call. = FALSE
      )
    }
    
    n_observations <- nrow(
      trajectory_data
    )
    
    aic_value <- -2 * log_likelihood +
      2 * n_parameters
    
    bic_value <- -2 * log_likelihood +
      log(n_observations) *
      n_parameters
    
    list(
      k = k_value,
      fit = fit,
      assignments = assignments,
      cluster_sizes = cluster_sizes,
      cluster_trajectories = cluster_trajectories,
      ots_distribution = ots_distribution,
      purity_by_cluster =
        purity_result$cluster_group_purity,
      purity_overall =
        purity_result$overall,
      minimum_cluster_size =
        minimum_observed_cluster_size,
      bic = bic_value,
      aic = aic_value,
      logLik = log_likelihood,
      n_parameters = as.numeric(
        n_parameters
      ),
      n_observations = n_observations
    )
  }
  
  # ---------------------------------------------------------------------------
  # Fit all candidate values of k
  # ---------------------------------------------------------------------------
  
  fitted_results <- lapply(
    k,
    fit_one_k
  )
  
  names(fitted_results) <- as.character(
    k
  )
  
  # ---------------------------------------------------------------------------
  # Build model-comparison table
  # ---------------------------------------------------------------------------
  
  trajectory_estimator_value <- if (
    is.null(trajectory_estimator)
  ) {
    NA_character_
  } else {
    trajectory_estimator
  }
  
  model_comparison <- data.frame(
    k = k,
    trajectory_source = rep(
      trajectory_source,
      length(k)
    ),
    trajectory_estimator = rep(
      trajectory_estimator_value,
      length(k)
    ),
    BIC = vapply(
      fitted_results,
      function(result) {
        result$bic
      },
      numeric(1)
    ),
    AIC = vapply(
      fitted_results,
      function(result) {
        result$aic
      },
      numeric(1)
    ),
    logLik = vapply(
      fitted_results,
      function(result) {
        result$logLik
      },
      numeric(1)
    ),
    n_parameters = vapply(
      fitted_results,
      function(result) {
        result$n_parameters
      },
      numeric(1)
    ),
    minimum_purity = vapply(
      fitted_results,
      function(result) {
        result$purity_overall$minimum_purity
      },
      numeric(1)
    ),
    mean_purity = vapply(
      fitted_results,
      function(result) {
        result$purity_overall$mean_purity
      },
      numeric(1)
    ),
    weighted_purity = vapply(
      fitted_results,
      function(result) {
        result$purity_overall$weighted_purity
      },
      numeric(1)
    ),
    passing_fraction = vapply(
      fitted_results,
      function(result) {
        result$purity_overall$passing_fraction
      },
      numeric(1)
    ),
    minimum_cluster_size = vapply(
      fitted_results,
      function(result) {
        result$minimum_cluster_size
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
  
  model_comparison <- model_comparison |>
    dplyr::mutate(
      passes_topology_rule =
        .data$passing_fraction >=
        .env$required_passing_fraction,
      passes_size_rule =
        .data$minimum_cluster_size >=
        .env$min_cluster_size,
      acceptable =
        .data$passes_topology_rule &
        .data$passes_size_rule
    )
  
  # ---------------------------------------------------------------------------
  # Select the cluster resolution
  # ---------------------------------------------------------------------------
  
  acceptable_models <- model_comparison |>
    dplyr::filter(
      .data$acceptable
    ) |>
    dplyr::arrange(
      .data$k
    )
  
  if (nrow(acceptable_models) > 0L) {
    
    selected_k <- acceptable_models$k[[1L]]
    
    selection_reason <- paste0(
      "Smallest k satisfying the topology-purity ",
      "and minimum-cluster-size criteria."
    )
    
  } else {
    
    selected_k <- model_comparison |>
      dplyr::arrange(
        .data$BIC
      ) |>
      dplyr::slice_head(
        n = 1L
      ) |>
      dplyr::pull(
        .data$k
      )
    
    selection_reason <- paste0(
      "No candidate satisfied both topology and size criteria; ",
      "the lowest-BIC solution was selected."
    )
  }
  
  selected_result <- fitted_results[[
    as.character(selected_k)
  ]]
  
  # ---------------------------------------------------------------------------
  # Store selected and comparison results
  # ---------------------------------------------------------------------------
  
  object@results$flexmix_model_comparison <-
    S4Vectors::DataFrame(
      model_comparison
    )
  
  object@results$flexmix_selected_k <-
    selected_k
  
  object@results$flexmix_selection_reason <-
    selection_reason
  
  object@results$flexmix_gene_assignments <-
    S4Vectors::DataFrame(
      selected_result$assignments
    )
  
  object@results$flexmix_cluster_trajectories <-
    S4Vectors::DataFrame(
      selected_result$cluster_trajectories
    )
  
  object@results$flexmix_ots_distribution <-
    S4Vectors::DataFrame(
      selected_result$ots_distribution
    )
  
  object@results$flexmix_cluster_purity <-
    S4Vectors::DataFrame(
      selected_result$purity_by_cluster
    )
  
  object@results$flexmix_overall_purity <-
    S4Vectors::DataFrame(
      selected_result$purity_overall
    )
  
  object@results$flexmix_trajectory_data <-
    S4Vectors::DataFrame(
      trajectory_data
    )
  
  if (isTRUE(store_fits)) {
    object@results$flexmix_fits <-
      fitted_results
  }
  
  # ---------------------------------------------------------------------------
  # Record history
  # ---------------------------------------------------------------------------
  
  object@history <- c(
    object@history,
    list(
      list(
        step = "lot_from_flexmix_clusters",
        assay = assay,
        assay_data = assay_data,
        group = group,
        visits = visits,
        candidate_k = k,
        selected_k = selected_k,
        trajectory_source = trajectory_source,
        trajectory_estimator =
          trajectory_estimator,
        model = model,
        scale_features = scale_features,
        purity_threshold = purity_threshold,
        required_passing_fraction =
          required_passing_fraction,
        min_cluster_size = min_cluster_size,
        nrep = nrep,
        seed = seed,
        n_features = n_features,
        selection_reason = selection_reason,
        time = Sys.time()
      )
    )
  )
  
  if (isTRUE(verbose)) {
    message(
      "Selected k = ",
      selected_k,
      ". ",
      selection_reason
    )
  }
  
  methods::validObject(
    object
  )
  
  object
}