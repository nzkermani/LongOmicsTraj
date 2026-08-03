#' Plot one clustered molecular-feature trajectory profile
#'
#' Creates a spaghetti plot for a single molecular-feature cluster.
#' Individual feature trajectories are coloured by their Omics Trajectory
#' Signature, while the representative cluster trajectory is displayed in
#' black.
#'
#' @param object A `LongOmicsTraj` object containing cluster relationships,
#'   gene-level OTS results, cluster trajectories and OTS metrics.
#' @param cluster_id Identifier of the cluster to display.
#' @param assay_name Name of the assay to display.
#' @param visits_vector Character vector giving visits in chronological order.
#'
#' @return A `ggplot` object.
#'
#' @export
plot_cluster_spaghetti_single <- function(
    object,
    cluster_id,
    assay_name,
    visits_vector
) {
  
  # ---------------------------------------------------------------------------
  # Validate inputs
  # ---------------------------------------------------------------------------
  
  if (!methods::is(object, "LongOmicsTraj")) {
    stop(
      "`object` must be a LongOmicsTraj object.",
      call. = FALSE
    )
  }
  
  if (
    length(cluster_id) != 1L ||
    is.na(cluster_id)
  ) {
    stop(
      "`cluster_id` must contain one non-missing cluster identifier.",
      call. = FALSE
    )
  }
  
  cluster_id <- as.character(
    cluster_id
  )
  
  if (!nzchar(cluster_id)) {
    stop(
      "`cluster_id` must not be empty.",
      call. = FALSE
    )
  }
  
  if (
    length(assay_name) != 1L ||
    !is.character(assay_name) ||
    is.na(assay_name) ||
    !nzchar(assay_name)
  ) {
    stop(
      "`assay_name` must be one non-empty character string.",
      call. = FALSE
    )
  }
  
  if (
    !is.character(visits_vector) ||
    length(visits_vector) < 2L ||
    anyNA(visits_vector) ||
    any(!nzchar(visits_vector)) ||
    anyDuplicated(visits_vector)
  ) {
    stop(
      "`visits_vector` must contain at least two unique visit labels.",
      call. = FALSE
    )
  }
  
  experiments <- MultiAssayExperiment::experiments(
    object@mae
  )
  
  if (!assay_name %in% names(experiments)) {
    stop(
      "Assay `",
      assay_name,
      "` was not found.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Relationships
  # ---------------------------------------------------------------------------
  
  relationship_data <- as.data.frame(
    object@relationships
  )
  
  required_relationship_columns <- c(
    "assay",
    "child_id",
    "parent_id"
  )
  
  missing_relationship_columns <- setdiff(
    required_relationship_columns,
    colnames(relationship_data)
  )
  
  if (length(missing_relationship_columns) > 0L) {
    stop(
      "Relationship data are missing: ",
      paste(
        missing_relationship_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  relationship_data <- relationship_data[
    as.character(relationship_data$assay) == assay_name,
    ,
    drop = FALSE
  ]
  
  if (nrow(relationship_data) == 0L) {
    stop(
      "No cluster relationships were found for assay `",
      assay_name,
      "`.",
      call. = FALSE
    )
  }
  
  relationship_data$child_id <- toupper(
    as.character(relationship_data$child_id)
  )
  
  relationship_data$parent_id <- as.character(
    relationship_data$parent_id
  )
  
  genes_in_cluster <- relationship_data |>
    dplyr::filter(
      .data$parent_id == .env$cluster_id
    ) |>
    dplyr::pull(
      .data$child_id
    ) |>
    unique()
  
  if (length(genes_in_cluster) == 0L) {
    stop(
      "No features were assigned to cluster `",
      cluster_id,
      "`.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Expression matrix and metadata
  # ---------------------------------------------------------------------------
  
  se <- experiments[[assay_name]]
  
  raw_matrix <- SummarizedExperiment::assay(
    se
  )
  
  metadata <- as.data.frame(
    SummarizedExperiment::colData(se)
  )
  
  if (!is.numeric(raw_matrix)) {
    stop(
      "The assay matrix must be numeric.",
      call. = FALSE
    )
  }
  
  required_metadata_columns <- c(
    "visit",
    "group"
  )
  
  missing_metadata_columns <- setdiff(
    required_metadata_columns,
    colnames(metadata)
  )
  
  if (length(missing_metadata_columns) > 0L) {
    stop(
      "Metadata are missing: ",
      paste(
        missing_metadata_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  if (is.null(rownames(raw_matrix))) {
    stop(
      "The assay matrix must have feature row names.",
      call. = FALSE
    )
  }
  
  rownames(raw_matrix) <- toupper(
    rownames(raw_matrix)
  )
  
  if ("sample_id" %in% colnames(metadata)) {
    metadata$sample_id <- as.character(
      metadata$sample_id
    )
  } else {
    metadata$sample_id <- rownames(
      metadata
    )
  }
  
  metadata <- metadata[
    as.character(metadata$visit) %in% visits_vector,
    ,
    drop = FALSE
  ]
  
  if (nrow(metadata) == 0L) {
    stop(
      "No samples remained after filtering the requested visits.",
      call. = FALSE
    )
  }
  
  groups <- unique(
    as.character(metadata$group)
  )
  
  # ---------------------------------------------------------------------------
  # Build the trajectory matrix
  # ---------------------------------------------------------------------------
  
  trajectory_list <- lapply(
    groups,
    function(current_group) {
      
      group_metadata <- metadata[
        as.character(metadata$group) == current_group,
        ,
        drop = FALSE
      ]
      
      sapply(
        visits_vector,
        function(current_visit) {
          
          sample_ids <- group_metadata$sample_id[
            as.character(group_metadata$visit) ==
              current_visit
          ]
          
          if (length(sample_ids) == 0L) {
            return(
              rep(
                NA_real_,
                nrow(raw_matrix)
              )
            )
          }
          
          rowMeans(
            raw_matrix[
              ,
              sample_ids,
              drop = FALSE
            ],
            na.rm = TRUE
          )
        }
      )
    }
  )
  
  trajectory_matrix <- do.call(
    cbind,
    trajectory_list
  )
  
  colnames(trajectory_matrix) <- as.vector(
    outer(
      groups,
      visits_vector,
      paste,
      sep = "_"
    )
  )
  
  trajectory_matrix <- trajectory_matrix[
    rownames(trajectory_matrix) %in% genes_in_cluster,
    ,
    drop = FALSE
  ]
  
  if (nrow(trajectory_matrix) == 0L) {
    stop(
      "No features from cluster `",
      cluster_id,
      "` were found in the assay matrix.",
      call. = FALSE
    )
  }
  
  trajectory_scaled <- t(
    base::scale(
      t(trajectory_matrix)
    )
  )
  
  complete_features <- rowSums(
    is.finite(trajectory_scaled)
  ) == ncol(trajectory_scaled)
  
  trajectory_scaled <- trajectory_scaled[
    complete_features,
    ,
    drop = FALSE
  ]
  
  if (nrow(trajectory_scaled) == 0L) {
    stop(
      "No complete, variable feature trajectories remained.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Convert trajectories to long format
  # ---------------------------------------------------------------------------
  
  gene_long <- as.data.frame(
    trajectory_scaled
  ) |>
    tibble::rownames_to_column(
      var = "gene_id"
    ) |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of("gene_id"),
      names_to = "combined",
      values_to = "expression"
    ) |>
    tidyr::separate(
      col = .data$combined,
      into = c(
        "group",
        "visit"
      ),
      sep = "_",
      extra = "merge",
      fill = "right"
    ) |>
    dplyr::mutate(
      cluster = .env$cluster_id
    )
  
  gene_long$visit <- factor(
    gene_long$visit,
    levels = visits_vector,
    ordered = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # Add topology labels
  # ---------------------------------------------------------------------------
  
  if (
    is.null(object@results$gene_ots) ||
    nrow(object@results$gene_ots) == 0L
  ) {
    stop(
      "No gene-level OTS results were found in ",
      "`object@results$gene_ots`.",
      call. = FALSE
    )
  }
  
  gene_ots <- as.data.frame(
    object@results$gene_ots
  )
  
  required_ots_columns <- c(
    "object_id",
    "topology_label"
  )
  
  missing_ots_columns <- setdiff(
    required_ots_columns,
    colnames(gene_ots)
  )
  
  if (length(missing_ots_columns) > 0L) {
    stop(
      "Gene-level OTS results are missing: ",
      paste(
        missing_ots_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  gene_ots$object_id <- toupper(
    as.character(gene_ots$object_id)
  )
  
  gene_long <- gene_long |>
    dplyr::left_join(
      gene_ots[
        ,
        c(
          "object_id",
          "topology_label"
        ),
        drop = FALSE
      ],
      by = c(
        "gene_id" = "object_id"
      )
    )
  
  # ---------------------------------------------------------------------------
  # Cluster-centre trajectory
  # ---------------------------------------------------------------------------
  
  vmm_data <- as.data.frame(
    object@vmm
  )
  
  required_vmm_columns <- c(
    "assay",
    "object_id",
    "visit",
    "estimated_value"
  )
  
  missing_vmm_columns <- setdiff(
    required_vmm_columns,
    colnames(vmm_data)
  )
  
  if (length(missing_vmm_columns) > 0L) {
    stop(
      "VMM results are missing: ",
      paste(
        missing_vmm_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  vmm_data <- vmm_data |>
    dplyr::filter(
      .data$assay == .env$assay_name,
      as.character(.data$object_id) ==
        .env$cluster_id
    ) |>
    dplyr::transmute(
      cluster = as.character(
        .data$object_id
      ),
      visit = .data$visit,
      expression = .data$estimated_value
    )
  
  if (nrow(vmm_data) == 0L) {
    stop(
      "No representative trajectory was found for cluster `",
      cluster_id,
      "`.",
      call. = FALSE
    )
  }
  
  vmm_data$visit <- factor(
    vmm_data$visit,
    levels = visits_vector,
    ordered = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # Build topology labels
  # ---------------------------------------------------------------------------
  
  topology_distribution <- gene_long |>
    dplyr::filter(
      !is.na(.data$topology_label)
    ) |>
    dplyr::count(
      .data$topology_label,
      name = "n"
    ) |>
    dplyr::mutate(
      pct = round(
        100 * .data$n /
          sum(.data$n),
        1
      ),
      label = paste0(
        .data$topology_label,
        " (",
        .data$pct,
        "%)"
      )
    )
  
  gene_long <- gene_long |>
    dplyr::left_join(
      topology_distribution[
        ,
        c(
          "topology_label",
          "label"
        ),
        drop = FALSE
      ],
      by = "topology_label"
    )
  
  # ---------------------------------------------------------------------------
  # Cluster metrics
  # ---------------------------------------------------------------------------
  
  if (
    is.null(object@results$ots_metrics) ||
    nrow(object@results$ots_metrics) == 0L
  ) {
    stop(
      "No OTS metrics were found in `object@results$ots_metrics`.",
      call. = FALSE
    )
  }
  
  metrics <- as.data.frame(
    object@results$ots_metrics
  )
  
  required_metric_columns <- c(
    "cluster_id",
    "purity",
    "entropy"
  )
  
  missing_metric_columns <- setdiff(
    required_metric_columns,
    colnames(metrics)
  )
  
  if (length(missing_metric_columns) > 0L) {
    stop(
      "OTS metrics are missing: ",
      paste(
        missing_metric_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  cluster_metrics <- metrics |>
    dplyr::filter(
      as.character(.data$cluster_id) ==
        .env$cluster_id
    )
  
  purity_value <- if (nrow(cluster_metrics) > 0L) {
    round(
      cluster_metrics$purity[[1L]],
      2
    )
  } else {
    NA_real_
  }
  
  entropy_value <- if (nrow(cluster_metrics) > 0L) {
    round(
      cluster_metrics$entropy[[1L]],
      2
    )
  } else {
    NA_real_
  }
  
  # ---------------------------------------------------------------------------
  # Build the plot
  # ---------------------------------------------------------------------------
  
  plot_object <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = gene_long,
      mapping = ggplot2::aes(
        x = .data$visit,
        y = .data$expression,
        group = .data$gene_id,
        colour = .data$label
      ),
      alpha = 0.4,
      linewidth = 0.5
    ) +
    ggplot2::geom_line(
      data = vmm_data,
      mapping = ggplot2::aes(
        x = .data$visit,
        y = .data$expression,
        group = .data$cluster
      ),
      colour = "black",
      linewidth = 1.5
    ) +
    ggplot2::geom_point(
      data = vmm_data,
      mapping = ggplot2::aes(
        x = .data$visit,
        y = .data$expression
      ),
      colour = "black",
      size = 3
    ) +
    ggplot2::scale_colour_brewer(
      palette = "Set2",
      na.value = "grey70"
    ) +
    ggplot2::labs(
      title = paste(
        "Cluster:",
        cluster_id
      ),
      subtitle = paste0(
        "n = ",
        dplyr::n_distinct(
          gene_long$gene_id
        ),
        " | Purity = ",
        purity_value,
        " | Entropy = ",
        entropy_value
      ),
      x = "Visit",
      y = "Scaled expression",
      colour = "Topology (%)"
    ) +
    ggplot2::theme_minimal(
      base_size = 13
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "right"
    )
  
  plot_object
}