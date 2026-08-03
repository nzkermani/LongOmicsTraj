#' Plot topology deviations within clusters
#'
#' Creates a PAMR-inspired nail plot showing how strongly each Omics
#' Trajectory Signature is enriched or depleted within each cluster compared
#' with its overall frequency.
#'
#' For topology `t` in cluster `c`, deviation is calculated as:
#'
#' \deqn{
#' D_{ct} = p_{ct} - p_t
#' }
#'
#' where `p_ct` is the within-cluster topology proportion and `p_t` is the
#' overall topology proportion.
#'
#' @param object A `LongOmicsTraj` object containing selected FlexMix
#'   clustering results.
#' @param clusters Optional vector of clusters to show.
#' @param group Optional clinical group.
#' @param min_abs_deviation Minimum absolute deviation required for display.
#' @param top_n Optional maximum number of topologies displayed per cluster.
#' @param sort_topologies Logical; order topologies by deviation within each
#'   cluster.
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
#'       "up_down",
#'       "flat_flat"
#'     ),
#'     n_genes = c(
#'       8L,
#'       2L,
#'       3L,
#'       7L
#'     ),
#'     proportion = c(
#'       0.8,
#'       0.2,
#'       0.3,
#'       0.7
#'     )
#'   )
#'
#' methods::slot(
#'   example_object,
#'   "results"
#' ) <- example_results
#'
#' plot <- lot_plot_topology_deviation(
#'   object = example_object,
#'   group = "treated"
#' )
#'
#' plot
lot_plot_topology_deviation <- function(
    object,
    clusters = NULL,
    group = NULL,
    min_abs_deviation = 0,
    top_n = NULL,
    sort_topologies = TRUE
) {
  
  if (!methods::is(object, "LongOmicsTraj")) {
    stop(
      "`object` must be a LongOmicsTraj object.",
      call. = FALSE
    )
  }
  
  if (
    !"flexmix_ots_distribution" %in%
    names(object@results)
  ) {
    stop(
      "No FlexMix OTS distribution was found. ",
      "Run `lot_from_flexmix_clusters()` first.",
      call. = FALSE
    )
  }
  
  if (
    length(min_abs_deviation) != 1L ||
    !is.numeric(min_abs_deviation) ||
    is.na(min_abs_deviation) ||
    !is.finite(min_abs_deviation) ||
    min_abs_deviation < 0 ||
    min_abs_deviation > 1
  ) {
    stop(
      "`min_abs_deviation` must be between 0 and 1.",
      call. = FALSE
    )
  }
  
  if (!is.null(top_n)) {
    
    if (
      length(top_n) != 1L ||
      !is.numeric(top_n) ||
      is.na(top_n) ||
      top_n < 1L ||
      top_n != as.integer(top_n)
    ) {
      stop(
        "`top_n` must be one positive integer or `NULL`.",
        call. = FALSE
      )
    }
    
    top_n <- as.integer(
      top_n
    )
  }
  
  ots_distribution <- as.data.frame(
    object@results$flexmix_ots_distribution
  )
  
  required_columns <- c(
    "cluster",
    "group",
    "topology_label",
    "n_genes",
    "proportion"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(ots_distribution)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "FlexMix OTS results are missing: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  available_groups <- unique(
    as.character(
      ots_distribution$group
    )
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
  
  plot_data <- ots_distribution |>
    dplyr::filter(
      as.character(group) == .env$group
    ) |>
    dplyr::mutate(
      cluster = as.character(cluster),
      topology_label = as.character(
        topology_label
      )
    )
  
  if (!is.null(clusters)) {
    
    requested_clusters <- as.character(
      clusters
    )
    
    plot_data <- plot_data |>
      dplyr::filter(
        cluster %in%
          requested_clusters
      )
  }
  
  if (nrow(plot_data) == 0L) {
    stop(
      "No OTS observations remained after filtering.",
      call. = FALSE
    )
  }
  
  overall_distribution <- plot_data |>
    dplyr::group_by(
      topology_label
    ) |>
    dplyr::summarise(
      overall_count = sum(
        n_genes
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      overall_proportion =
        overall_count /
        sum(overall_count)
    )
  
  all_clusters <- unique(
    plot_data$cluster
  )
  
  all_topologies <- unique(
    plot_data$topology_label
  )
  
  complete_grid <- expand.grid(
    cluster = all_clusters,
    topology_label = all_topologies,
    stringsAsFactors = FALSE
  )
  
  plot_data <- complete_grid |>
    dplyr::left_join(
      plot_data |>
        dplyr::select(
          cluster,
          topology_label,
          n_genes,
          proportion
        ),
      by = c(
        "cluster",
        "topology_label"
      )
    ) |>
    dplyr::mutate(
      n_genes = dplyr::coalesce(
        n_genes,
        0L
      ),
      proportion = dplyr::coalesce(
        proportion,
        0
      )
    ) |>
    dplyr::left_join(
      overall_distribution,
      by = "topology_label"
    ) |>
    dplyr::mutate(
      deviation =
        proportion -
        overall_proportion,
      direction = dplyr::if_else(
        deviation >= 0,
        "Enriched",
        "Depleted"
      )
    ) |>
    dplyr::filter(
      abs(deviation) >=
        min_abs_deviation
    )
  
  if (!is.null(top_n)) {
    plot_data <- plot_data |>
      dplyr::group_by(
        cluster
      ) |>
      dplyr::slice_max(
        order_by = abs(deviation),
        n = top_n,
        with_ties = FALSE
      ) |>
      dplyr::ungroup()
  }
  
  if (nrow(plot_data) == 0L) {
    stop(
      "No topology deviations passed the requested filter.",
      call. = FALSE
    )
  }
  
  cluster_sizes <- ots_distribution |>
    dplyr::filter(
      as.character(group) == .env$group
    ) |>
    dplyr::group_by(
      cluster
    ) |>
    dplyr::summarise(
      n_features = sum(
        n_genes
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      cluster = as.character(
        cluster
      ),
      cluster_label = paste0(
        "Cluster ",
        cluster,
        " u00b7 n = ",
        n_features
      )
    )
  
  plot_data <- plot_data |>
    dplyr::left_join(
      cluster_sizes,
      by = "cluster"
    )
  
  if (isTRUE(sort_topologies)) {
    
    plot_data <- plot_data |>
      dplyr::group_by(
        cluster
      ) |>
      dplyr::arrange(
        deviation,
        .by_group = TRUE
      ) |>
      dplyr::mutate(
        topology_plot_label = factor(
          topology_label,
          levels = unique(
            topology_label
          )
        )
      ) |>
      dplyr::ungroup()
    
  } else {
    
    plot_data$topology_plot_label <- factor(
      plot_data$topology_label,
      levels = sort(
        unique(
          plot_data$topology_label
        )
      )
    )
  }
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      y = topology_plot_label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.5
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = deviation,
        yend = topology_plot_label,
        linetype = direction
      ),
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = deviation,
        shape = direction
      ),
      size = 2.4
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = deviation,
        label = scales::percent(
          deviation,
          accuracy = 0.1
        )
      ),
      hjust = dplyr::if_else(
        plot_data$deviation >= 0,
        -0.15,
        1.15
      ),
      size = 3
    ) +
    ggplot2::facet_wrap(
      ~cluster_label,
      scales = "free_y"
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(),
      expand = ggplot2::expansion(
        mult = c(
          0.15,
          0.15
        )
      )
    ) +
    ggplot2::labs(
      title = "Topology-deviation nail plot",
      subtitle = paste0(
        "Positive deviations indicate OTS enrichment within a cluster; ",
        "negative deviations indicate depletion relative to the full ",
        "clustered feature set."
      ),
      x = "Cluster proportion u2212 overall proportion",
      y = "Omics Trajectory Signature",
      linetype = NULL,
      shape = NULL
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
      panel.grid.minor =
        ggplot2::element_blank(),
      panel.grid.major.y =
        ggplot2::element_blank(),
      legend.position = "bottom"
    )
}
