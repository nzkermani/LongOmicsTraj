#' Plot trajectory topology profiles for fitted clusters
#'
#' Displays the representative longitudinal trajectory and the
#' within-cluster Omics Trajectory Signature composition for each selected
#' FlexMix cluster.
#'
#' The left-hand panel shows the representative cluster trajectory and,
#' optionally, the individual molecular-feature trajectories. The right-hand
#' panel shows the proportion of features assigned to each OTS class.
#'
#' @param object A `LongOmicsTraj` object containing FlexMix clustering
#'   results.
#' @param clusters Optional vector of cluster identifiers. By default, all
#'   clusters in the selected FlexMix solution are shown.
#' @param group Optional clinical group. When omitted, the group stored in the
#'   selected FlexMix results is used.
#' @param value_type Trajectory values to display: `"scaled"` uses the
#'   clustering values, while `"raw"` uses the original estimated trajectory
#'   values.
#' @param show_features Logical; display individual feature trajectories.
#' @param max_feature_lines Maximum number of individual feature trajectories
#'   drawn per cluster.
#' @param purity_threshold Reference purity threshold shown in the OTS panel.
#' @param sort_by Ordering of clusters: `"cluster"`, `"purity"`, `"size"` or
#'   `"dominant_topology"`.
#' @param seed Random seed used when subsampling feature trajectories.
#'
#' @return A `patchwork` plot object.
#'
#' @export
#'
#' @examples
#' data("lot_glucold_ots")
#'
#' example_object <- lot_glucold_ots
#'
#' example_results <- lot_results(
#'   example_object
#' )
#'
#' example_results$flexmix_gene_assignments <-
#'   S4Vectors::DataFrame(
#'     gene_id = c(
#'       "gene_1",
#'       "gene_2",
#'       "gene_3",
#'       "gene_4"
#'     ),
#'     cluster = c(
#'       "1",
#'       "1",
#'       "2",
#'       "2"
#'     ),
#'     group = "treated"
#'   )
#'
#' example_results$flexmix_cluster_trajectories <-
#'   S4Vectors::DataFrame(
#'     cluster = rep(
#'       c("1", "2"),
#'       each = 3
#'     ),
#'     group = "treated",
#'     visit = rep(
#'       c(
#'         "Baseline",
#'         "Month 6",
#'         "Month 12"
#'       ),
#'       times = 2
#'     ),
#'     mean_clustering_value = c(
#'       0, 1, 0,
#'       0, -1, 0.5
#'     ),
#'     mean_trajectory_value = c(
#'       5, 6, 5,
#'       5, 4, 5.5
#'     )
#'   )
#'
#' example_results$flexmix_ots_distribution <-
#'   S4Vectors::DataFrame(
#'     cluster = c(
#'       "1",
#'       "1",
#'       "2",
#'       "2"
#'     ),
#'     group = "treated",
#'     topology_label = c(
#'       "up_down",
#'       "flat_flat",
#'       "down_up",
#'       "flat_flat"
#'     ),
#'     n_genes = c(
#'       3L,
#'       1L,
#'       3L,
#'       1L
#'     ),
#'     proportion = c(
#'       0.75,
#'       0.25,
#'       0.75,
#'       0.25
#'     )
#'   )
#'
#' example_results$flexmix_cluster_purity <-
#'   S4Vectors::DataFrame(
#'     cluster = c(
#'       "1",
#'       "2"
#'     ),
#'     group = "treated",
#'     n_genes = c(
#'       2L,
#'       2L
#'     ),
#'     dominant_topology = c(
#'       "up_down",
#'       "down_up"
#'     ),
#'     purity = c(
#'       0.75,
#'       0.75
#'     )
#'   )
#'
#' example_results$flexmix_trajectory_data <-
#'   S4Vectors::DataFrame(
#'     gene_id = rep(
#'       c(
#'         "gene_1",
#'         "gene_2",
#'         "gene_3",
#'         "gene_4"
#'       ),
#'       each = 3
#'     ),
#'     visit = rep(
#'       c(
#'         "Baseline",
#'         "Month 6",
#'         "Month 12"
#'       ),
#'       times = 4
#'     ),
#'     clustering_value = c(
#'       0, 1.1, 0.1,
#'       0, 0.9, -0.1,
#'       0, -1.1, 0.6,
#'       0, -0.9, 0.4
#'     ),
#'     mean_expression = c(
#'       5, 6.1, 5.1,
#'       5, 5.9, 4.9,
#'       5, 3.9, 5.6,
#'       5, 4.1, 5.4
#'     )
#'   )
#'
#' methods::slot(
#'   example_object,
#'   "results"
#' ) <- example_results
#'
#' plot <- lot_plot_cluster_topology(
#'   object = example_object,
#'   group = "treated",
#'   show_features = FALSE
#' )
#'
#' plot
lot_plot_cluster_topology <- function(
    object,
    clusters = NULL,
    group = NULL,
    value_type = c("scaled", "raw"),
    show_features = TRUE,
    max_feature_lines = 50L,
    purity_threshold = 0.90,
    sort_by = c(
      "cluster",
      "purity",
      "size",
      "dominant_topology"
    ),
    seed = 1L
) {
  
  value_type <- match.arg(
    value_type
  )
  
  sort_by <- match.arg(
    sort_by
  )
  
  if (!methods::is(object, "LongOmicsTraj")) {
    stop(
      "`object` must be a LongOmicsTraj object.",
      call. = FALSE
    )
  }
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Package `patchwork` is required for this plot.",
      call. = FALSE
    )
  }
  
  if (
    length(max_feature_lines) != 1L ||
    !is.numeric(max_feature_lines) ||
    is.na(max_feature_lines) ||
    max_feature_lines < 0L ||
    max_feature_lines != as.integer(max_feature_lines)
  ) {
    stop(
      "`max_feature_lines` must be one non-negative integer.",
      call. = FALSE
    )
  }
  
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
  
  required_results <- c(
    "flexmix_gene_assignments",
    "flexmix_cluster_trajectories",
    "flexmix_ots_distribution",
    "flexmix_cluster_purity",
    "flexmix_trajectory_data"
  )
  
  missing_results <- setdiff(
    required_results,
    names(object@results)
  )
  
  if (length(missing_results) > 0L) {
    stop(
      "FlexMix results are missing: ",
      paste(
        missing_results,
        collapse = ", "
      ),
      ". Run `lot_from_flexmix_clusters()` first.",
      call. = FALSE
    )
  }
  
  assignments <- as.data.frame(
    object@results$flexmix_gene_assignments
  )
  
  cluster_trajectories <- as.data.frame(
    object@results$flexmix_cluster_trajectories
  )
  
  ots_distribution <- as.data.frame(
    object@results$flexmix_ots_distribution
  )
  
  purity_data <- as.data.frame(
    object@results$flexmix_cluster_purity
  )
  
  feature_trajectories <- as.data.frame(
    object@results$flexmix_trajectory_data
  )
  
  available_groups <- unique(
    as.character(purity_data$group)
  )
  
  if (is.null(group)) {
    
    if (length(available_groups) != 1L) {
      stop(
        "Multiple groups are available. Select one using `group`.",
        call. = FALSE
      )
    }
    
    group <- available_groups[[1L]]
  }
  
  purity_data <- purity_data |>
    dplyr::filter(
      as.character(.data$group) == .env$group
    )
  
  assignments <- assignments |>
    dplyr::filter(
      as.character(.data$group) == .env$group
    )
  
  cluster_trajectories <- cluster_trajectories |>
    dplyr::filter(
      as.character(.data$group) == .env$group
    )
  
  ots_distribution <- ots_distribution |>
    dplyr::filter(
      as.character(.data$group) == .env$group
    )
  
  if (nrow(purity_data) == 0L) {
    stop(
      "No FlexMix cluster results were found for group `",
      group,
      "`.",
      call. = FALSE
    )
  }
  
  available_clusters <- unique(
    as.character(purity_data$cluster)
  )
  
  if (!is.null(clusters)) {
    
    requested_clusters <- as.character(
      clusters
    )
    
    missing_clusters <- setdiff(
      requested_clusters,
      available_clusters
    )
    
    if (length(missing_clusters) > 0L) {
      warning(
        "These clusters were not found and will be omitted: ",
        paste(
          missing_clusters,
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    
    available_clusters <- intersect(
      requested_clusters,
      available_clusters
    )
  }
  
  if (length(available_clusters) == 0L) {
    stop(
      "No clusters remained for plotting.",
      call. = FALSE
    )
  }
  
  purity_data <- purity_data |>
    dplyr::filter(
      as.character(.data$cluster) %in%
        .env$available_clusters
    )
  
  assignments <- assignments |>
    dplyr::filter(
      as.character(.data$cluster) %in%
        .env$available_clusters
    )
  
  cluster_trajectories <- cluster_trajectories |>
    dplyr::filter(
      as.character(.data$cluster) %in%
        .env$available_clusters
    )
  
  ots_distribution <- ots_distribution |>
    dplyr::filter(
      as.character(.data$cluster) %in%
        .env$available_clusters
    )
  
  if (sort_by == "purity") {
    
    purity_data <- purity_data |>
      dplyr::arrange(
        dplyr::desc(.data$purity)
      )
    
  } else if (sort_by == "size") {
    
    purity_data <- purity_data |>
      dplyr::arrange(
        dplyr::desc(.data$n_genes)
      )
    
  } else if (sort_by == "dominant_topology") {
    
    purity_data <- purity_data |>
      dplyr::arrange(
        .data$dominant_topology,
        dplyr::desc(.data$purity)
      )
    
  } else {
    
    suppressWarnings(
      cluster_number <- as.numeric(
        as.character(purity_data$cluster)
      )
    )
    
    if (all(is.finite(cluster_number))) {
      
      purity_data$cluster_number <-
        cluster_number
      
      purity_data <- purity_data |>
        dplyr::arrange(
          .data$cluster_number
        )
      
    } else {
      
      purity_data <- purity_data |>
        dplyr::arrange(
          as.character(.data$cluster)
        )
    }
  }
  
  cluster_order <- as.character(
    purity_data$cluster
  )
  
  purity_lookup <- purity_data |>
    dplyr::transmute(
      cluster = as.character(.data$cluster),
      cluster_label = paste0(
        "Cluster ",
        .data$cluster,
        "  \u00b7  n = ",
        .data$n_genes,
        "  \u00b7  ",
        .data$dominant_topology,
        "  \u00b7  purity ",
        sprintf(
          "%.1f%%",
          100 * .data$purity
        )
      )
    )
  
  cluster_label_levels <- purity_lookup$cluster_label
  
  cluster_trajectories <- cluster_trajectories |>
    dplyr::mutate(
      cluster = as.character(.data$cluster)
    ) |>
    dplyr::left_join(
      purity_lookup,
      by = "cluster"
    )
  
  ots_distribution <- ots_distribution |>
    dplyr::mutate(
      cluster = as.character(.data$cluster)
    ) |>
    dplyr::left_join(
      purity_lookup,
      by = "cluster"
    )
  
  assignments <- assignments |>
    dplyr::mutate(
      gene_id = as.character(.data$gene_id),
      cluster = as.character(.data$cluster)
    )
  
  feature_trajectories <- feature_trajectories |>
    dplyr::mutate(
      gene_id = as.character(.data$gene_id)
    ) |>
    dplyr::inner_join(
      assignments |>
        dplyr::select(
          dplyr::all_of(
            c(
              "gene_id",
              "cluster"
            )
          )
        ) |>
        dplyr::distinct(),
      by = "gene_id"
    ) |>
    dplyr::left_join(
      purity_lookup,
      by = "cluster"
    )
  
  cluster_trajectories$cluster_label <- factor(
    cluster_trajectories$cluster_label,
    levels = cluster_label_levels
  )
  
  ots_distribution$cluster_label <- factor(
    ots_distribution$cluster_label,
    levels = cluster_label_levels
  )
  
  feature_trajectories$cluster_label <- factor(
    feature_trajectories$cluster_label,
    levels = cluster_label_levels
  )
  
  trajectory_y <- if (value_type == "scaled") {
    "mean_clustering_value"
  } else {
    "mean_trajectory_value"
  }
  
  feature_y <- if (value_type == "scaled") {
    "clustering_value"
  } else {
    "mean_expression"
  }
  
  if (
    !trajectory_y %in%
    colnames(cluster_trajectories)
  ) {
    stop(
      "Cluster trajectory results do not contain `",
      trajectory_y,
      "`.",
      call. = FALSE
    )
  }
  
  trajectory_plot <- ggplot2::ggplot()
  
  if (
    isTRUE(show_features) &&
    max_feature_lines > 0L
  ) {
    
    selected_features <- withr::with_seed(
      as.integer(seed),
      feature_trajectories |>
        dplyr::distinct(
          .data$cluster,
          .data$gene_id
        ) |>
        dplyr::group_by(
          .data$cluster
        ) |>
        dplyr::group_modify(
          function(.x, .y) {
            
            n_to_sample <- min(
              max_feature_lines,
              nrow(.x)
            )
            
            dplyr::slice_sample(
              .x,
              n = n_to_sample
            )
          }
        ) |>
        dplyr::ungroup()
    )
    
    
    feature_plot_data <- feature_trajectories |>
      dplyr::semi_join(
        selected_features,
        by = c(
          "cluster",
          "gene_id"
        )
      )
    
    trajectory_plot <- trajectory_plot +
      ggplot2::geom_line(
        data = feature_plot_data,
        ggplot2::aes(
          x = .data$visit,
          y = .data[[feature_y]],
          group = .data$gene_id
        ),
        linewidth = 0.3,
        alpha = 0.12
      )
  }
  
  trajectory_plot <- trajectory_plot +
    ggplot2::geom_line(
      data = cluster_trajectories,
      ggplot2::aes(
        x = .data$visit,
        y = .data[[trajectory_y]],
        group = .data$cluster
      ),
      linewidth = 1.15
    ) +
    ggplot2::geom_point(
      data = cluster_trajectories,
      ggplot2::aes(
        x = .data$visit,
        y = .data[[trajectory_y]]
      ),
      size = 2
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(
        .data[["cluster_label"]]
      ),
      scales = "free_y"
    ) +
    ggplot2::labs(
      title = "Representative cluster trajectories",
      x = "Visit",
      y = if (value_type == "scaled") {
        "Scaled trajectory"
      } else {
        "Estimated trajectory"
      }
    ) +
    ggplot2::theme_bw(
      base_size = 10
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      strip.text.y = ggplot2::element_text(
        angle = 0,
        hjust = 0
      ),
      axis.text.x = ggplot2::element_text(
        angle = 30,
        hjust = 1
      ),
      panel.grid.minor =
        ggplot2::element_blank()
    )
  
  composition_plot <- ggplot2::ggplot(
    ots_distribution,
    ggplot2::aes(
      x = .data$proportion,
      y = stats::reorder(
        .data$topology_label,
        .data$proportion
      )
    )
  ) +
    ggplot2::geom_col(
      width = 0.7
    ) +
    ggplot2::geom_vline(
      xintercept = purity_threshold,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(
          .data$n_genes,
          " (",
          sprintf(
            "%.1f%%",
            100 * .data$proportion
          ),
          ")"
        )
      ),
      hjust = -0.08,
      size = 3
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(
        .data[["cluster_label"]]
      ),
      scales = "free_y"
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(),
      limits = c(
        0,
        1.12
      ),
      breaks = c(
        0,
        0.25,
        0.50,
        0.75,
        1
      )
    ) +
    ggplot2::labs(
      title = "Member trajectory topologies",
      x = "Within-cluster proportion",
      y = "OTS"
    ) +
    ggplot2::theme_bw(
      base_size = 10
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      strip.text.y =
        ggplot2::element_blank(),
      panel.grid.minor =
        ggplot2::element_blank()
    )
  
  combined_plot <- trajectory_plot +
    composition_plot +
    patchwork::plot_layout(
      widths = c(
        1,
        1.25
      )
    ) +
    patchwork::plot_annotation(
      title = paste0(
        "LongOmicsTraj cluster topology profile: ",
        group
      ),
      subtitle = paste0(
        "Representative trajectories are shown beside their member ",
        "OTS composition; the dashed line marks ",
        sprintf(
          "%.0f%%",
          100 * purity_threshold
        ),
        " purity."
      )
    )
  
  combined_plot
}
