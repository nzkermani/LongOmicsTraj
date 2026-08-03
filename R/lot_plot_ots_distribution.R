#' Plot the distribution of Omics Trajectory Signatures
#'
#' Summarises the proportion or number of molecular features assigned
#' to each Omics Trajectory Signature within each group.
#'
#' @param object A `LongOmicsTraj` object containing computed OTS results.
#' @param estimator Optional estimator name used to filter the OTS results.
#' @param type Either `"proportion"` or `"count"`.
#' @param palette Name of an RColorBrewer qualitative palette.
#'
#' @return A `ggplot` object.
#'
#' @export
#'
#' @examples
#' data("lot_glucold_ots")
#'
#' plot <- lot_plot_ots_distribution(
#'   object = lot_glucold_ots,
#'   estimator = "masigpro",
#'   type = "proportion"
#' )
#'
#' plot
lot_plot_ots_distribution <- function(
    object,
    estimator = NULL,
    type = c(
      "proportion",
      "count"
    ),
    palette = "Set3"
) {
  
  type <- match.arg(
    type
  )
  
  if (!methods::is(object, "LongOmicsTraj")) {
    stop(
      "`object` must be a LongOmicsTraj object.",
      call. = FALSE
    )
  }
  
  if (
    !is.character(palette) ||
    length(palette) != 1L ||
    is.na(palette) ||
    !nzchar(palette)
  ) {
    stop(
      "`palette` must be one non-empty character string.",
      call. = FALSE
    )
  }
  
  ots_data <- as.data.frame(
    object@ots
  )
  
  if (nrow(ots_data) == 0L) {
    stop(
      "No OTS results are stored in `object`. ",
      "Run `lot_map_to_ots()` first.",
      call. = FALSE
    )
  }
  
  required_columns <- c(
    "object_id",
    "group",
    "topology_label"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(ots_data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "The OTS results are missing: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
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
    
    available_estimators <- available_estimators[
      nzchar(available_estimators)
    ]
    
    if (!estimator %in% available_estimators) {
      stop(
        "Estimator `",
        estimator,
        "` was not found. Available estimators: ",
        paste(
          available_estimators,
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    
    ots_data <- ots_data[
      as.character(ots_data$estimator) == estimator,
      ,
      drop = FALSE
    ]
  }
  
  ots_data <- ots_data[
    !is.na(ots_data$group) &
      !is.na(ots_data$topology_label),
    ,
    drop = FALSE
  ]
  
  if (nrow(ots_data) == 0L) {
    stop(
      "No complete OTS results were available for plotting.",
      call. = FALSE
    )
  }
  
  plot_data <- ots_data |>
    dplyr::count(
      .data$group,
      .data$topology_label,
      name = "count"
    )
  
  if (type == "proportion") {
    
    plot_data <- plot_data |>
      dplyr::group_by(
        .data$group
      ) |>
      dplyr::mutate(
        value = .data$count /
          sum(.data$count)
      ) |>
      dplyr::ungroup()
    
    y_scale <- ggplot2::scale_y_continuous(
      labels = scales::percent_format(),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.02
        )
      )
    )
    
    y_label <- "Proportion of molecular features"
    
  } else {
    
    plot_data$value <- plot_data$count
    
    y_scale <- ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.05
        )
      )
    )
    
    y_label <- "Number of molecular features"
  }
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$group,
      y = .data$value,
      fill = .data$topology_label
    )
  ) +
    ggplot2::geom_col(
      width = 0.72,
      colour = "white",
      linewidth = 0.25
    ) +
    y_scale +
    ggplot2::scale_fill_brewer(
      palette = palette,
      name = "OTS"
    ) +
    ggplot2::labs(
      title = "Distribution of Omics Trajectory Signatures",
      subtitle = if (is.null(estimator)) {
        NULL
      } else {
        paste(
          "Estimator:",
          estimator
        )
      },
      x = "Group",
      y = y_label
    ) +
    ggplot2::theme_bw(
      base_size = 11
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 35,
        hjust = 1
      ),
      panel.grid.minor =
        ggplot2::element_blank(),
      legend.position = "bottom"
    )
}