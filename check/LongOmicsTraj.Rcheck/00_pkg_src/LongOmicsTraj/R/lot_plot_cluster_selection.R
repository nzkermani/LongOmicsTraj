#' Plot topology-aware cluster selection
#'
#' Displays likelihood-based and topology-aware statistics across candidate
#' numbers of FlexMix clusters.
#'
#' @param object A `LongOmicsTraj` object containing FlexMix model-comparison
#'   results.
#' @param metrics Metrics to display. Supported values are `"BIC"`,
#'   `"weighted_purity"`, `"passing_fraction"`,
#'   `"minimum_cluster_size"`, `"minimum_purity"` and `"mean_purity"`.
#' @param purity_threshold Reference purity threshold.
#' @param required_passing_fraction Required fraction of clusters passing the
#'   purity rule.
#' @param min_cluster_size Minimum acceptable cluster size.
#' @param show_selected Logical; mark the selected cluster resolution.
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
#' example_results$flexmix_model_comparison <- S4Vectors::DataFrame(
#'   k = 2:4,
#'   BIC = c(540, 505, 512),
#'   weighted_purity = c(0.76, 0.91, 0.88),
#'   passing_fraction = c(0.50, 1.00, 0.75),
#'   minimum_cluster_size = c(18, 12, 6),
#'   minimum_purity = c(0.68, 0.87, 0.72),
#'   mean_purity = c(0.78, 0.92, 0.86)
#' )
#'
#' example_results$flexmix_selected_k <- 3L
#'
#' methods::slot(
#'   example_object,
#'   "results"
#' ) <- example_results
#'
#' plot <- lot_plot_cluster_selection(
#'   object = example_object
#' )
#'
#' plot
lot_plot_cluster_selection <- function(
    object,
    metrics = c(
      "BIC",
      "weighted_purity",
      "passing_fraction",
      "minimum_cluster_size"
    ),
    purity_threshold = 0.90,
    required_passing_fraction = 0.80,
    min_cluster_size = 5L,
    show_selected = TRUE
) {
  
  if (!methods::is(object, "LongOmicsTraj")) {
    stop(
      "`object` must be a LongOmicsTraj object.",
      call. = FALSE
    )
  }
  
  if (
    !"flexmix_model_comparison" %in%
    names(object@results)
  ) {
    stop(
      "No FlexMix model-comparison results were found. ",
      "Run `lot_from_flexmix_clusters()` first.",
      call. = FALSE
    )
  }
  
  comparison <- as.data.frame(
    object@results$flexmix_model_comparison
  )
  
  supported_metrics <- c(
    "BIC",
    "weighted_purity",
    "passing_fraction",
    "minimum_cluster_size",
    "minimum_purity",
    "mean_purity"
  )
  
  if (
    !is.character(metrics) ||
    length(metrics) == 0L ||
    anyNA(metrics) ||
    any(!nzchar(metrics))
  ) {
    stop(
      "`metrics` must contain one or more non-empty metric names.",
      call. = FALSE
    )
  }
  
  metrics <- unique(
    metrics
  )
  
  unsupported_metrics <- setdiff(
    metrics,
    supported_metrics
  )
  
  if (length(unsupported_metrics) > 0L) {
    stop(
      "Unsupported metrics: ",
      paste(
        unsupported_metrics,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  missing_metrics <- setdiff(
    metrics,
    colnames(comparison)
  )
  
  if (length(missing_metrics) > 0L) {
    stop(
      "Model-comparison results are missing: ",
      paste(
        missing_metrics,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  if (
    !"k" %in% colnames(comparison) ||
    !is.numeric(comparison$k)
  ) {
    stop(
      "Model-comparison results must contain a numeric `k` column.",
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
    !is.finite(min_cluster_size) ||
    min_cluster_size < 1L ||
    min_cluster_size != as.integer(min_cluster_size)
  ) {
    stop(
      "`min_cluster_size` must be one positive integer.",
      call. = FALSE
    )
  }
  
  if (
    !is.logical(show_selected) ||
    length(show_selected) != 1L ||
    is.na(show_selected)
  ) {
    stop(
      "`show_selected` must be `TRUE` or `FALSE`.",
      call. = FALSE
    )
  }
  
  min_cluster_size <- as.integer(
    min_cluster_size
  )
  
  metric_labels <- c(
    BIC = "BIC \u00b7 lower is better",
    weighted_purity =
      "Weighted topology purity \u00b7 higher is better",
    passing_fraction =
      "Clusters passing purity rule \u00b7 higher is better",
    minimum_cluster_size =
      "Smallest cluster \u00b7 must remain adequate",
    minimum_purity =
      "Minimum cluster purity \u00b7 higher is better",
    mean_purity =
      "Mean cluster purity \u00b7 higher is better"
  )
  
  plot_data <- dplyr::bind_rows(
    lapply(
      metrics,
      function(metric_name) {
        data.frame(
          k = comparison$k,
          metric = metric_labels[[metric_name]],
          metric_name = metric_name,
          value = comparison[[metric_name]],
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  plot_data$metric <- factor(
    plot_data$metric,
    levels = unname(
      metric_labels[
        metrics
      ]
    )
  )
  
  threshold_data <- data.frame(
    metric = factor(
      character(0),
      levels = levels(
        plot_data$metric
      )
    ),
    threshold = numeric(0),
    stringsAsFactors = FALSE
  )
  
  if ("passing_fraction" %in% metrics) {
    threshold_data <- rbind(
      threshold_data,
      data.frame(
        metric = factor(
          metric_labels[[
            "passing_fraction"
          ]],
          levels = levels(
            plot_data$metric
          )
        ),
        threshold = required_passing_fraction,
        stringsAsFactors = FALSE
      )
    )
  }
  
  if ("minimum_cluster_size" %in% metrics) {
    threshold_data <- rbind(
      threshold_data,
      data.frame(
        metric = factor(
          metric_labels[[
            "minimum_cluster_size"
          ]],
          levels = levels(
            plot_data$metric
          )
        ),
        threshold = min_cluster_size,
        stringsAsFactors = FALSE
      )
    )
  }
  
  if ("minimum_purity" %in% metrics) {
    threshold_data <- rbind(
      threshold_data,
      data.frame(
        metric = factor(
          metric_labels[[
            "minimum_purity"
          ]],
          levels = levels(
            plot_data$metric
          )
        ),
        threshold = purity_threshold,
        stringsAsFactors = FALSE
      )
    )
  }
  
  selected_k <- NULL
  
  if (
    "flexmix_selected_k" %in%
    names(object@results)
  ) {
    selected_k <- as.integer(
      object@results$flexmix_selected_k
    )
    
    if (
      length(selected_k) != 1L ||
      is.na(selected_k)
    ) {
      selected_k <- NULL
    }
  }
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$k,
      y = .data$value
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      size = 2
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(
        .data$metric
      ),
      scales = "free_y"
    ) +
    ggplot2::scale_x_continuous(
      breaks = sort(
        unique(
          comparison$k
        )
      )
    ) +
    ggplot2::labs(
      title = "Topology-aware cluster-resolution path",
      subtitle = paste0(
        "Likelihood fit, trajectory homogeneity and cluster size ",
        "are evaluated together."
      ),
      x = "Candidate number of clusters",
      y = NULL
    ) +
    ggplot2::theme_bw(
      base_size = 11
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      strip.text.y = ggplot2::element_text(
        angle = 0,
        hjust = 0
      ),
      panel.grid.minor =
        ggplot2::element_blank()
    )
  
  if (nrow(threshold_data) > 0L) {
    p <- p +
      ggplot2::geom_hline(
        data = threshold_data,
        ggplot2::aes(
          yintercept = .data$threshold
        ),
        inherit.aes = FALSE,
        linetype = "dashed",
        linewidth = 0.5
      )
  }
  
  if (
    isTRUE(show_selected) &&
    !is.null(selected_k)
  ) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = selected_k,
        linetype = "dotted",
        linewidth = 0.8
      ) +
      ggplot2::annotate(
        geom = "text",
        x = selected_k,
        y = Inf,
        label = paste0(
          "Selected k = ",
          selected_k
        ),
        vjust = 1.4,
        hjust = -0.08,
        size = 3.2
      )
  }
  
  p
}