#' Plot trajectory arms within molecular-feature clusters
#'
#' Visualises the distribution of feature-level Omics Trajectory
#' Signatures within each molecular-feature cluster. Each panel represents
#' one cluster and clinical group. Arrows emerge from a shared cluster node
#' and represent the distinct trajectory classes found inside that cluster.
#'
#' Arrow width is proportional to the fraction of clustered features
#' assigned to each trajectory class.
#'
#' @param object A `LongOmicsTraj` object containing
#'   `cluster_ots_distribution` results.
#' @param clusters Optional character vector of cluster identifiers to show.
#'   By default, all clusters are included.
#' @param groups Optional character vector of clinical or treatment groups
#'   to show. By default, all groups are included.
#' @param top_n Optional maximum number of trajectory arms displayed within
#'   each cluster-group panel. By default, all trajectory classes are shown.
#' @param min_proportion Minimum within-cluster trajectory proportion required
#'   for display. The default is `0`.
#' @param show_counts Logical; include feature counts in trajectory labels.
#' @param show_percent Logical; include percentages in trajectory labels.
#' @param scale_width Logical; scale arrow width by trajectory proportion.
#' @param ncol Number of facet columns.
#'
#' @return A `ggplot` object.
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
#' example_results$cluster_ots_distribution <-
#'   S4Vectors::DataFrame(
#'     cluster_id = c(
#'       "cluster_1",
#'       "cluster_1",
#'       "cluster_2",
#'       "cluster_2"
#'     ),
#'     group = "treated",
#'     topology_label = c(
#'       "up_down",
#'       "flat_flat",
#'       "down_up",
#'       "flat_down"
#'     ),
#'     count = c(
#'       8L,
#'       2L,
#'       7L,
#'       3L
#'     ),
#'     proportion = c(
#'       0.8,
#'       0.2,
#'       0.7,
#'       0.3
#'     )
#'   )
#'
#' methods::slot(
#'   example_object,
#'   "results"
#' ) <- example_results
#'
#' plot <- lot_plot_cluster_trajectory_arms(
#'   object = example_object,
#'   groups = "treated"
#' )
#'
#' plot
lot_plot_cluster_trajectory_arms <- function(
    object,
    clusters = NULL,
    groups = NULL,
    top_n = NULL,
    min_proportion = 0,
    show_counts = TRUE,
    show_percent = TRUE,
    scale_width = TRUE,
    ncol = NULL
) {
  
  # ---------------------------------------------------------------------------
  # Validate input
  # ---------------------------------------------------------------------------
  
  if (!methods::is(object, "LongOmicsTraj")) {
    stop(
      "`object` must be a LongOmicsTraj object.",
      call. = FALSE
    )
  }
  
  if (
    length(min_proportion) != 1L ||
    is.na(min_proportion) ||
    !is.numeric(min_proportion) ||
    min_proportion < 0 ||
    min_proportion > 1
  ) {
    stop(
      "`min_proportion` must be a number between 0 and 1.",
      call. = FALSE
    )
  }
  
  if (!is.null(top_n)) {
    
    if (
      length(top_n) != 1L ||
      is.na(top_n) ||
      !is.numeric(top_n) ||
      top_n < 1L ||
      top_n != as.integer(top_n)
    ) {
      stop(
        "`top_n` must be one positive integer or `NULL`.",
        call. = FALSE
      )
    }
    
    top_n <- as.integer(top_n)
  }
  
  if (
    is.null(object@results$cluster_ots_distribution)
  ) {
    stop(
      "No within-cluster OTS distribution was found. ",
      "Run `lot_from_clusters()` on an object containing feature-level ",
      "OTS results first.",
      call. = FALSE
    )
  }
  
  plot_data <- as.data.frame(
    object@results$cluster_ots_distribution
  )
  
  if (nrow(plot_data) == 0L) {
    stop(
      "The within-cluster OTS distribution is empty.",
      call. = FALSE
    )
  }
  
  required_columns <- c(
    "cluster_id",
    "group",
    "topology_label",
    "count",
    "proportion"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(plot_data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Cluster OTS results are missing: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Filter clusters
  # ---------------------------------------------------------------------------
  
  available_clusters <- unique(
    as.character(plot_data$cluster_id)
  )
  
  if (!is.null(clusters)) {
    
    clusters <- unique(
      as.character(clusters)
    )
    
    missing_clusters <- setdiff(
      clusters,
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
    
    clusters <- intersect(
      clusters,
      available_clusters
    )
    
    if (length(clusters) == 0L) {
      stop(
        "None of the requested clusters were found.",
        call. = FALSE
      )
    }
    
    plot_data <- plot_data |>
      dplyr::filter(
        .data$cluster_id %in% .env$clusters
      )
    
    plot_data$cluster_id <- factor(
      plot_data$cluster_id,
      levels = clusters
    )
  }
  
  # ---------------------------------------------------------------------------
  # Filter groups
  # ---------------------------------------------------------------------------
  
  available_groups <- unique(
    as.character(plot_data$group)
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
        paste(
          missing_groups,
          collapse = ", "
        ),
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
    
    plot_data <- plot_data |>
      dplyr::filter(
        .data$group %in% .env$groups
      )
    
    plot_data$group <- factor(
      plot_data$group,
      levels = groups
    )
  }
  
  # ---------------------------------------------------------------------------
  # Apply minimum proportion and top-N filters
  # ---------------------------------------------------------------------------
  
  plot_data <- plot_data |>
    dplyr::filter(
      !is.na(.data$cluster_id),
      !is.na(.data$group),
      !is.na(.data$topology_label),
      !is.na(.data$proportion),
      .data$proportion >= .env$min_proportion
    ) |>
    dplyr::group_by(
      .data$cluster_id,
      .data$group
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$proportion),
      .data$topology_label,
      .by_group = TRUE
    )
  
  if (!is.null(top_n)) {
    plot_data <- plot_data |>
      dplyr::slice_head(
        n = top_n
      )
  }
  
  plot_data <- plot_data |>
    dplyr::ungroup()
  
  if (nrow(plot_data) == 0L) {
    stop(
      "No trajectory arms remained after filtering.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Create panel identifiers and arm positions
  # ---------------------------------------------------------------------------
  
  plot_data <- plot_data |>
    dplyr::mutate(
      panel_id = paste(
        .data$cluster_id,
        .data$group,
        sep = " | "
      )
    ) |>
    dplyr::group_by(
      .data$panel_id
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$proportion),
      .data$topology_label,
      .by_group = TRUE
    ) |>
    dplyr::mutate(
      n_arms = dplyr::n(),
      arm_index = dplyr::row_number(),
      end_y = if (dplyr::n() == 1L) {
        0
      } else {
        seq(
          from = 1,
          to = -1,
          length.out = dplyr::n()
        )
      }
    ) |>
    dplyr::ungroup()
  
  # ---------------------------------------------------------------------------
  # Build labels
  # ---------------------------------------------------------------------------
  
  plot_data$arm_label <- as.character(
    plot_data$topology_label
  )
  
  if (isTRUE(show_counts)) {
    plot_data$arm_label <- paste0(
      plot_data$arm_label,
      "\n",
      "n = ",
      plot_data$count
    )
  }
  
  if (isTRUE(show_percent)) {
    percentage_label <- paste0(
      round(
        100 * plot_data$proportion,
        1
      ),
      "%"
    )
    
    if (isTRUE(show_counts)) {
      plot_data$arm_label <- paste0(
        plot_data$arm_label,
        " (",
        percentage_label,
        ")"
      )
    } else {
      plot_data$arm_label <- paste0(
        plot_data$arm_label,
        "\n",
        percentage_label
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Calculate arrow width
  # ---------------------------------------------------------------------------
  
  if (isTRUE(scale_width)) {
    
    plot_data$arm_width <- 0.7 +
      5 * plot_data$proportion
    
  } else {
    
    plot_data$arm_width <- 1.4
  }
  
  node_data <- plot_data |>
    dplyr::distinct(
      .data$panel_id,
      .data$cluster_id,
      .data$group
    ) |>
    dplyr::mutate(
      node_x = 0,
      node_y = 0
    )
  
  # ---------------------------------------------------------------------------
  # Build plot
  # ---------------------------------------------------------------------------
  
  ggplot2::ggplot() +
    
    ggplot2::geom_segment(
      data = plot_data,
      ggplot2::aes(
        x = 0.08,
        y = 0,
        xend = 0.78,
        yend = .data$end_y,
        colour = .data$topology_label,
        linewidth = .data$arm_width
      ),
      lineend = "round",
      arrow = grid::arrow(
        length = grid::unit(
          0.16,
          "inches"
        ),
        type = "closed"
      )
    ) +
    
    ggplot2::geom_point(
      data = node_data,
      ggplot2::aes(
        x = .data$node_x,
        y = .data$node_y
      ),
      size = 8,
      shape = 21,
      fill = "grey20",
      colour = "grey20"
    ) +
    
    ggplot2::geom_text(
      data = node_data,
      ggplot2::aes(
        x = .data$node_x,
        y = .data$node_y - 0.28,
        label = .data$cluster_id
      ),
      fontface = "bold",
      size = 3.3
    ) +
    
    ggplot2::geom_label(
      data = plot_data,
      ggplot2::aes(
        x = 0.92,
        y = .data$end_y,
        label = .data$arm_label,
        fill = .data$topology_label
      ),
      hjust = 0,
      colour = "white",
      fontface = "bold",
      size = 3.1,
      linewidth  = 0,
      label.padding = grid::unit(
        0.16,
        "lines"
      ),
      show.legend = FALSE
    ) +
    
    ggplot2::facet_wrap(
      facets = ggplot2::vars(
        .data$panel_id
      ),
      ncol = ncol
    ) +
    
    ggplot2::scale_x_continuous(
      limits = c(
        -0.18,
        1.65
      ),
      expand = c(
        0,
        0
      )
    ) +
    
    ggplot2::scale_y_continuous(
      limits = c(
        -1.3,
        1.3
      ),
      expand = c(
        0,
        0
      )
    ) +
    
    ggplot2::scale_linewidth_identity() +
    
    ggplot2::labs(
      title = "Trajectory arms within molecular-feature clusters",
      subtitle = paste0(
        "One cluster may contain multiple longitudinal OTS patterns; ",
        "arrow width represents the within-cluster proportion."
      ),
      colour = "Trajectory class"
    ) +
    
    ggplot2::coord_cartesian(
      clip = "off"
    ) +
    
    ggplot2::theme_void(
      base_size = 11
    ) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 13
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10,
        margin = ggplot2::margin(
          b = 12
        )
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 10
      ),
      panel.spacing = grid::unit(
        1.3,
        "lines"
      ),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(
        10,
        75,
        10,
        10
      )
    )
}