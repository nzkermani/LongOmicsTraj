#' Plot clustered molecular-feature trajectories
#'
#' Creates one spaghetti plot per molecular-feature cluster. Individual
#' feature trajectories are coloured by their Omics Trajectory Signature,
#' while the cluster-level representative trajectory is shown in black.
#'
#' @param object A `LongOmicsTraj` object containing cluster relationships,
#'   gene-level OTS results, cluster trajectories and OTS metrics.
#' @param assay_name Name of the assay to display.
#' @param visits_vector Character vector giving visits in chronological order.
#' @param max_clusters Maximum number of clusters to plot.
#'
#' @return Invisibly returns a named list of `ggplot` objects.
#'
#' @export
plot_cluster_spaghetti <- function(
    object,
    assay_name,
    visits_vector,
    max_clusters = 3L
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
  
  if (
    length(max_clusters) != 1L ||
    !is.numeric(max_clusters) ||
    is.na(max_clusters) ||
    max_clusters < 1L ||
    max_clusters != as.integer(max_clusters)
  ) {
    stop(
      "`max_clusters` must be one positive integer.",
      call. = FALSE
    )
  }
  
  max_clusters <- as.integer(
    max_clusters
  )
  
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
  
  rel_df <- as.data.frame(
    object@relationships
  )
  
  required_relationship_columns <- c(
    "assay",
    "child_id",
    "parent_id"
  )
  
  missing_relationship_columns <- setdiff(
    required_relationship_columns,
    colnames(rel_df)
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
  
  rel_df <- rel_df[
    as.character(rel_df$assay) == assay_name,
    ,
    drop = FALSE
  ]
  
  if (nrow(rel_df) == 0L) {
    stop(
      "No cluster relationships were found for assay `",
      assay_name,
      "`.",
      call. = FALSE
    )
  }
  
  rel_df$child_id <- toupper(
    as.character(rel_df$child_id)
  )
  
  rel_df$parent_id <- as.character(
    rel_df$parent_id
  )
  
  # ---------------------------------------------------------------------------
  # Expression matrix and metadata
  # ---------------------------------------------------------------------------
  
  se <- experiments[[assay_name]]
  
  raw_mat <- SummarizedExperiment::assay(
    se
  )
  
  meta_df <- as.data.frame(
    SummarizedExperiment::colData(se)
  )
  
  if (!is.numeric(raw_mat)) {
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
    colnames(meta_df)
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
  
  if (is.null(rownames(raw_mat))) {
    stop(
      "The assay matrix must have feature row names.",
      call. = FALSE
    )
  }
  
  rownames(raw_mat) <- toupper(
    rownames(raw_mat)
  )
  
  if ("sample_id" %in% colnames(meta_df)) {
    meta_df$sample_id <- as.character(
      meta_df$sample_id
    )
  } else {
    meta_df$sample_id <- rownames(
      meta_df
    )
  }
  
  meta_df <- meta_df[
    as.character(meta_df$visit) %in% visits_vector,
    ,
    drop = FALSE
  ]
  
  groups <- unique(
    as.character(meta_df$group)
  )
  
  if (length(groups) == 0L) {
    stop(
      "No groups remained after filtering the requested visits.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Build feature-level trajectory matrix
  # ---------------------------------------------------------------------------
  
  traj_list <- lapply(
    groups,
    function(current_group) {
      
      meta_group <- meta_df[
        as.character(meta_df$group) == current_group,
        ,
        drop = FALSE
      ]
      
      sapply(
        visits_vector,
        function(current_visit) {
          
          samples <- meta_group$sample_id[
            as.character(meta_group$visit) ==
              current_visit
          ]
          
          if (length(samples) == 0L) {
            return(
              rep(
                NA_real_,
                nrow(raw_mat)
              )
            )
          }
          
          rowMeans(
            raw_mat[
              ,
              samples,
              drop = FALSE
            ],
            na.rm = TRUE
          )
        }
      )
    }
  )
  
  traj_mat <- do.call(
    cbind,
    traj_list
  )
  
  colnames(traj_mat) <- as.vector(
    outer(
      groups,
      visits_vector,
      paste,
      sep = "_"
    )
  )
  
  keep_genes <- unique(
    rel_df$child_id
  )
  
  traj_mat <- traj_mat[
    rownames(traj_mat) %in% keep_genes,
    ,
    drop = FALSE
  ]
  
  if (nrow(traj_mat) == 0L) {
    stop(
      "No clustered features were found in the assay matrix.",
      call. = FALSE
    )
  }
  
  traj_scaled <- t(
    base::scale(
      t(traj_mat)
    )
  )
  
  keep_variable <- rowSums(
    is.finite(traj_scaled)
  ) == ncol(traj_scaled)
  
  traj_scaled <- traj_scaled[
    keep_variable,
    ,
    drop = FALSE
  ]
  
  if (nrow(traj_scaled) == 0L) {
    stop(
      "No complete, variable feature trajectories remained.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Convert trajectories to long format
  # ---------------------------------------------------------------------------
  
  gene_long <- as.data.frame(
    traj_scaled
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
    dplyr::left_join(
      rel_df[
        ,
        c(
          "child_id",
          "parent_id"
        ),
        drop = FALSE
      ],
      by = c(
        "gene_id" = "child_id"
      )
    ) |>
    dplyr::rename(
      cluster = .data$parent_id
    )
  
  gene_long$visit <- factor(
    gene_long$visit,
    levels = visits_vector,
    ordered = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # Add feature-level topology labels
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
  # Cluster-centre trajectories
  # ---------------------------------------------------------------------------
  
  vmm_df <- as.data.frame(
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
    colnames(vmm_df)
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
  
  vmm_df <- vmm_df |>
    dplyr::filter(
      .data$assay == .env$assay_name
    ) |>
    dplyr::transmute(
      cluster = as.character(
        .data$object_id
      ),
      visit = .data$visit,
      expression = .data$estimated_value
    )
  
  vmm_df$visit <- factor(
    vmm_df$visit,
    levels = visits_vector,
    ordered = TRUE
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
  
  # ---------------------------------------------------------------------------
  # Build one plot per cluster
  # ---------------------------------------------------------------------------
  
  unique_clusters <- sort(
    unique(
      as.character(gene_long$cluster)
    )
  )
  
  unique_clusters <- utils::head(
    unique_clusters,
    max_clusters
  )
  
  plot_list <- vector(
    mode = "list",
    length = length(unique_clusters)
  )
  
  names(plot_list) <- unique_clusters
  
  for (cluster_value in unique_clusters) {
    
    cluster_genes <- gene_long |>
      dplyr::filter(
        .data$cluster == .env$cluster_value
      )
    
    cluster_vmm <- vmm_df |>
      dplyr::filter(
        .data$cluster == .env$cluster_value
      )
    
    topology_distribution <- cluster_genes |>
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
    
    cluster_genes <- cluster_genes |>
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
    
    cluster_metrics <- metrics |>
      dplyr::filter(
        as.character(.data$cluster_id) ==
          .env$cluster_value
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
    
    plot_object <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = cluster_genes,
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
        data = cluster_vmm,
        mapping = ggplot2::aes(
          x = .data$visit,
          y = .data$expression,
          group = .data$cluster
        ),
        colour = "black",
        linewidth = 1.5
      ) +
      ggplot2::geom_point(
        data = cluster_vmm,
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
          cluster_value
        ),
        subtitle = paste0(
          "n = ",
          dplyr::n_distinct(
            cluster_genes$gene_id
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
    
    print(
      plot_object
    )
    
    plot_list[[cluster_value]] <- plot_object
  }
  
  invisible(
    plot_list
  )
}