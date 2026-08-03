#' Plot longitudinal trajectories for selected molecular features
#'
#' Plots observed longitudinal measurements for one or more molecular
#' features across clinical groups. Individual-subject profiles, group
#' means, uncertainty intervals and assigned Omics Trajectory Signatures
#' can be displayed.
#'
#' @param object A `LongOmicsTraj` object.
#' @param features Character vector of feature identifiers, such as gene symbols.
#' @param assay Name of the experiment in the `MultiAssayExperiment`.
#' @param assay_data Optional name of the matrix within the
#'   `SummarizedExperiment`. By default, the first available matrix is used.
#' @param groups Optional character vector of groups to display. By default,
#'   all available groups are included.
#' @param estimator Optional estimator used to retrieve OTS labels.
#' @param show_subjects Logical; if `TRUE`, individual-subject trajectories
#'   are displayed as faint lines.
#' @param summary One of `"mean_se"`, `"mean_ci"` or `"mean"`.
#' @param facet_scales Either `"free_y"` or `"fixed"`.
#' @param ncol Number of facet columns.
#'
#' @return A `ggplot` object.
#'
#' @export
#'
#' @examples
#' data("lot_glucold_ots")
#'
#' lot_plot_trajectory(
#'   lot_glucold_ots,
#'   features = c("ATP2A3", "CAT"),
#'   groups = c(
#'     "ICS 30 months",
#'     "ICS 6 months then withdrawal",
#'     "Placebo"
#'   ),
#'   estimator = "masigpro",
#'   show_subjects = TRUE
#' )
lot_plot_trajectory <- function(
    object,
    features,
    assay = "transcriptomics",
    assay_data = NULL,
    groups = NULL,
    estimator = NULL,
    show_subjects = FALSE,
    summary = c("mean_se", "mean_ci", "mean"),
    facet_scales = c("free_y", "fixed"),
    ncol = NULL
) {
  
  summary <- match.arg(summary)
  facet_scales <- match.arg(facet_scales)
  
  # ---------------------------------------------------------------------------
  # Validate the input object
  # ---------------------------------------------------------------------------
  
  if (!methods::is(object, "LongOmicsTraj")) {
    stop(
      "`object` must be a LongOmicsTraj object.",
      call. = FALSE
    )
  }
  
  if (
    missing(features) ||
    length(features) == 0L ||
    anyNA(features)
  ) {
    stop(
      "Provide at least one valid feature in `features`.",
      call. = FALSE
    )
  }
  
  features <- unique(
    as.character(features)
  )
  
  # ---------------------------------------------------------------------------
  # Extract the requested experiment
  # ---------------------------------------------------------------------------
  
  mae <- object@mae
  
  experiment_names <- names(
    MultiAssayExperiment::experiments(mae)
  )
  
  if (!assay %in% experiment_names) {
    stop(
      "Assay `", assay, "` was not found. Available assays: ",
      paste(experiment_names, collapse = ", "),
      call. = FALSE
    )
  }
  
  se <- MultiAssayExperiment::experiments(
    mae
  )[[assay]]
  
  available_matrices <- SummarizedExperiment::assayNames(
    se
  )
  
  if (length(available_matrices) == 0L) {
    stop(
      "The selected assay contains no data matrices.",
      call. = FALSE
    )
  }
  
  if (is.null(assay_data)) {
    assay_data <- available_matrices[[1L]]
  }
  
  if (!assay_data %in% available_matrices) {
    stop(
      "Data matrix `", assay_data, "` was not found. Available matrices: ",
      paste(available_matrices, collapse = ", "),
      call. = FALSE
    )
  }
  
  expression_matrix <- SummarizedExperiment::assay(
    se,
    assay_data
  )
  
  metadata <- as.data.frame(
    SummarizedExperiment::colData(se)
  )
  
  # ---------------------------------------------------------------------------
  # Validate metadata
  # ---------------------------------------------------------------------------
  
  required_metadata <- c(
    "sample_id",
    "subject_id",
    "visit",
    "group"
  )
  
  missing_metadata <- setdiff(
    required_metadata,
    colnames(metadata)
  )
  
  if (length(missing_metadata) > 0L) {
    stop(
      "Metadata are missing: ",
      paste(missing_metadata, collapse = ", "),
      call. = FALSE
    )
  }
  
  metadata$sample_id <- as.character(
    metadata$sample_id
  )
  
  if (!identical(
    colnames(expression_matrix),
    metadata$sample_id
  )) {
    
    metadata <- metadata[
      match(
        colnames(expression_matrix),
        metadata$sample_id
      ),
      ,
      drop = FALSE
    ]
  }
  
  if (
    anyNA(metadata$sample_id) ||
    !identical(
      colnames(expression_matrix),
      metadata$sample_id
    )
  ) {
    stop(
      "Expression samples could not be matched to the metadata.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Validate requested features
  # ---------------------------------------------------------------------------
  
  missing_features <- setdiff(
    features,
    rownames(expression_matrix)
  )
  
  if (length(missing_features) > 0L) {
    warning(
      "These features were not found and will be omitted: ",
      paste(missing_features, collapse = ", "),
      call. = FALSE
    )
  }
  
  features <- intersect(
    features,
    rownames(expression_matrix)
  )
  
  if (length(features) == 0L) {
    stop(
      "None of the requested features were found.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Convert expression data to long format
  # ---------------------------------------------------------------------------
  
  expression_long <- as.data.frame(
    expression_matrix[
      features,
      ,
      drop = FALSE
    ]
  ) |>
    tibble::rownames_to_column(
      "feature"
    ) |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of("feature"),
      names_to = "sample_id",
      values_to = "value"
    ) |>
    dplyr::left_join(
      metadata |>
        dplyr::select(
          dplyr::all_of(
            c(
              "sample_id",
              "subject_id",
              "visit",
              "group"
            )
          )
        ),
      by = "sample_id"
    )
  
  # Preserve the visit order stored in the metadata.
  if (is.factor(metadata$visit)) {
    expression_long$visit <- factor(
      expression_long$visit,
      levels = levels(metadata$visit),
      ordered = is.ordered(metadata$visit)
    )
  }
  
  # ---------------------------------------------------------------------------
  # Filter clinical groups
  # ---------------------------------------------------------------------------
  
  available_groups <- unique(
    as.character(expression_long$group)
  )
  
  if (!is.null(groups)) {
    
    groups <- unique(
      as.character(groups)
    )
    
    missing_groups <- setdiff(
      groups,
      available_groups
    )
    
    if (length(missing_groups) > 0L) {
      warning(
        "These groups were not found and will be omitted: ",
        paste(missing_groups, collapse = ", "),
        call. = FALSE
      )
    }
    
    groups <- intersect(
      groups,
      available_groups
    )
    
    if (length(groups) == 0L) {
      stop(
        "None of the requested groups were found.",
        call. = FALSE
      )
    }
    
    expression_long <- expression_long |>
      dplyr::filter(
        .data$group %in% .env$groups
      )
    
    expression_long$group <- factor(
      expression_long$group,
      levels = groups
    )
  }
  
  if (nrow(expression_long) == 0L) {
    stop(
      "No observations remained after filtering.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Calculate group summaries
  # ---------------------------------------------------------------------------
  
  summary_data <- expression_long |>
    dplyr::group_by(
      .data$feature,
      .data$group,
      .data$visit
    ) |>
    dplyr::summarise(
      n = sum(!is.na(.data$value)),
      mean_value = mean(
        .data$value,
        na.rm = TRUE
      ),
      sd_value = stats::sd(
        .data$value,
        na.rm = TRUE
      ),
      se_value = dplyr::if_else(
        .data$n > 1L,
        .data$sd_value / sqrt(.data$n),
        0
      ),
      .groups = "drop"
    )
  
  if (summary == "mean_se") {
    
    summary_data <- summary_data |>
      dplyr::mutate(
        lower = .data$mean_value - .data$se_value,
        upper = .data$mean_value + .data$se_value
      )
    
  } else if (summary == "mean_ci") {
    
    summary_data <- summary_data |>
      dplyr::mutate(
        lower = .data$mean_value - 1.96 * .data$se_value,
        upper = .data$mean_value + 1.96 * .data$se_value
      )
    
  } else {
    
    summary_data <- summary_data |>
      dplyr::mutate(
        lower = .data$mean_value,
        upper = .data$mean_value
      )
  }
  
  # ---------------------------------------------------------------------------
  # Add OTS labels to facets when one group is selected
  # ---------------------------------------------------------------------------
  
  feature_labels <- stats::setNames(
    features,
    features
  )
  
  if (
    !is.null(groups) &&
    length(groups) == 1L &&
    nrow(as.data.frame(object@ots)) > 0L
  ) {
    
    ots_data <- as.data.frame(
      object@ots
    )
    
    if (!is.null(estimator)) {
      
      if (!"estimator" %in% colnames(ots_data)) {
        stop(
          "The OTS results do not contain an `estimator` column.",
          call. = FALSE
        )
      }
      
      available_estimators <- sort(
        unique(
          as.character(
            ots_data$estimator[
              !is.na(ots_data$estimator)
            ]
          )
        )
      )
      
      if (!estimator %in% available_estimators) {
        stop(
          "Estimator `", estimator, "` was not found. Available estimators: ",
          paste(available_estimators, collapse = ", "),
          call. = FALSE
        )
      }
      
      ots_data <- ots_data[
        ots_data$estimator == estimator,
        ,
        drop = FALSE
      ]
    }
    
    ots_labels <- ots_data |>
      dplyr::filter(
        .data$object_id %in% .env$features,
        .data$group == .env$groups[[1L]]
      ) |>
      dplyr::distinct(
        .data$object_id,
        .data$topology_label
      )
    
    if (nrow(ots_labels) > 0L) {
      
      feature_labels[
        ots_labels$object_id
      ] <- paste0(
        ots_labels$object_id,
        "\nOTS: ",
        ots_labels$topology_label
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Build the plot
  # ---------------------------------------------------------------------------
  
  p <- ggplot2::ggplot()
  
  if (isTRUE(show_subjects)) {
    
    p <- p +
      ggplot2::geom_line(
        data = expression_long,
        ggplot2::aes(
          x = .data$visit,
          y = .data$value,
          group = interaction(
            .data$group,
            .data$subject_id
          ),
          colour = .data$group
        ),
        alpha = 0.12,
        linewidth = 0.3,
        show.legend = FALSE
      )
  }
  
  if (summary != "mean") {
    
    p <- p +
      ggplot2::geom_errorbar(
        data = summary_data,
        ggplot2::aes(
          x = .data$visit,
          ymin = .data$lower,
          ymax = .data$upper,
          colour = .data$group
        ),
        width = 0.08,
        linewidth = 0.45
      )
  }
  
  number_of_groups <- dplyr::n_distinct(
    expression_long$group
  )
  
  subtitle_text <- paste0(
    length(features),
    if (length(features) == 1L) {
      " molecular feature across "
    } else {
      " molecular features across "
    },
    number_of_groups,
    if (number_of_groups == 1L) {
      " group"
    } else {
      " groups"
    }
  )
  
  p +
    ggplot2::geom_line(
      data = summary_data,
      ggplot2::aes(
        x = .data$visit,
        y = .data$mean_value,
        group = .data$group,
        colour = .data$group
      ),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = summary_data,
      ggplot2::aes(
        x = .data$visit,
        y = .data$mean_value,
        colour = .data$group
      ),
      size = 2
    ) +
    ggplot2::facet_wrap(
      facets = ggplot2::vars(
        .data$feature
      ),
      scales = facet_scales,
      ncol = ncol,
      labeller = ggplot2::as_labeller(
        feature_labels
      )
    ) +
    ggplot2::labs(
      title = "Longitudinal molecular trajectories",
      subtitle = subtitle_text,
      x = "Visit",
      y = "Expression",
      colour = "Group"
    ) +
    ggplot2::theme_bw(
      base_size = 11
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 30,
        hjust = 1
      ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}