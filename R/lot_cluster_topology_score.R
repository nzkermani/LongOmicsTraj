#' Calculate a topology-aware cluster score
#'
#' Summarises topology purity, dominant-topology diversity and duplication
#' for a fitted clustering solution.
#'
#' This function complements likelihood-based criteria such as BIC by asking
#' whether the fitted clusters represent distinct longitudinal trajectory
#' behaviours. A solution with many clusters but few unique dominant
#' topologies has low topology coverage and may reflect over-partitioning.
#'
#' @param purity A `lot_cluster_purity` object returned by
#'   [lot_cluster_purity()], or a data frame containing cluster-level purity
#'   results.
#' @param k Optional requested number of clusters.
#' @param bic Optional BIC value for the fitted clustering solution.
#' @param aic Optional AIC value for the fitted clustering solution.
#' @param minimum_cluster_size Optional size of the smallest fitted cluster.
#'
#' @return An object of class `lot_cluster_topology_score` containing:
#' \describe{
#'   \item{group_summary}{
#'     Topology diversity and duplication statistics for each clinical group.
#'   }
#'   \item{overall}{
#'     Overall topology-aware summary of the clustering solution.
#'   }
#' }
#'
#' @details
#' For each clinical group, topology coverage is defined as:
#'
#' \deqn{
#'   \text{topology coverage}
#'   =
#'   \frac{\text{number of unique dominant topologies}}
#'        {\text{number of observed clusters}}
#' }
#'
#' A value of 1 indicates that every cluster has a different dominant
#' topology. Lower values indicate that multiple clusters share the same
#' dominant trajectory behaviour.
#'
#' The duplicate topology count is:
#'
#' \deqn{
#'   \text{number of clusters}
#'   -
#'   \text{number of unique dominant topologies}
#' }
#'
#' The function reports purity and coverage separately. It does not yet
#' combine BIC and topology metrics into a single weighted score.
#'
#' @export
#'
#' @examples
#' assignments <- data.frame(
#'   gene_id = paste0("G", 1:12),
#'   cluster = c(
#'     rep(1, 4),
#'     rep(2, 4),
#'     rep(3, 4)
#'   ),
#'   group = "Withdrawal",
#'   topology_label = c(
#'     rep("flat_flat", 4),
#'     rep("flat_flat", 3),
#'     "flat_up",
#'     rep("down_up", 4)
#'   )
#' )
#'
#' purity <- lot_cluster_purity(
#'   assignments,
#'   purity_threshold = 0.90
#' )
#'
#' score <- lot_cluster_topology_score(
#'   purity,
#'   k = 3,
#'   bic = 125.4,
#'   minimum_cluster_size = 4
#' )
#'
#' score$group_summary
#' score$overall
lot_cluster_topology_score <- function(
    purity,
    k = NULL,
    bic = NA_real_,
    aic = NA_real_,
    minimum_cluster_size = NA_integer_
) {
  
  # ---------------------------------------------------------------------------
  # Validate and extract purity results
  # ---------------------------------------------------------------------------
  
  if (inherits(purity, "lot_cluster_purity")) {
    
    purity_data <- purity$cluster_group_purity
    
  } else if (
    is.data.frame(purity) ||
    methods::is(purity, "DataFrame")
  ) {
    
    purity_data <- as.data.frame(
      purity
    )
    
  } else {
    
    stop(
      "`purity` must be a `lot_cluster_purity` object, a data frame, ",
      "or an S4Vectors DataFrame.",
      call. = FALSE
    )
  }
  
  purity_data <- as.data.frame(
    purity_data
  )
  
  required_columns <- c(
    "cluster",
    "group",
    "n_genes",
    "dominant_topology",
    "purity",
    "normalised_entropy",
    "passes_purity"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(purity_data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Purity results are missing: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  purity_data <- purity_data |>
    dplyr::filter(
      !is.na(cluster),
      !is.na(group),
      !is.na(dominant_topology),
      !is.na(purity),
      !is.na(n_genes)
    ) |>
    dplyr::mutate(
      group = as.character(group),
      dominant_topology = as.character(
        dominant_topology
      )
    )
  
  if (nrow(purity_data) == 0L) {
    stop(
      "No complete cluster-purity results were available.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Validate optional model information
  # ---------------------------------------------------------------------------
  
  if (!is.null(k)) {
    
    if (
      length(k) != 1L ||
      !is.numeric(k) ||
      is.na(k) ||
      !is.finite(k) ||
      k < 1L ||
      k != as.integer(k)
    ) {
      stop(
        "`k` must be one positive integer or `NULL`.",
        call. = FALSE
      )
    }
    
    k <- as.integer(k)
  }
  
  validate_optional_numeric <- function(
    value,
    argument_name
  ) {
    
    if (
      length(value) != 1L ||
      !is.numeric(value) ||
      is.na(value)
    ) {
      if (
        length(value) == 1L &&
        is.numeric(value) &&
        is.na(value)
      ) {
        return(invisible(TRUE))
      }
      
      stop(
        "`",
        argument_name,
        "` must be one numeric value or `NA`.",
        call. = FALSE
      )
    }
    
    if (!is.finite(value)) {
      stop(
        "`",
        argument_name,
        "` must be finite or `NA`.",
        call. = FALSE
      )
    }
    
    invisible(TRUE)
  }
  
  validate_optional_numeric(
    bic,
    "bic"
  )
  
  validate_optional_numeric(
    aic,
    "aic"
  )
  
  if (
    length(minimum_cluster_size) != 1L ||
    !is.numeric(minimum_cluster_size) ||
    (
      !is.na(minimum_cluster_size) &&
      (
        !is.finite(minimum_cluster_size) ||
        minimum_cluster_size < 1L ||
        minimum_cluster_size !=
        as.integer(minimum_cluster_size)
      )
    )
  ) {
    stop(
      "`minimum_cluster_size` must be one positive integer or `NA`.",
      call. = FALSE
    )
  }
  
  if (!is.na(minimum_cluster_size)) {
    minimum_cluster_size <- as.integer(
      minimum_cluster_size
    )
  }
  
  # ---------------------------------------------------------------------------
  # Calculate topology diversity within each group
  # ---------------------------------------------------------------------------
  
  group_summary <- purity_data |>
    dplyr::group_by(
      group
    ) |>
    dplyr::summarise(
      observed_clusters = dplyr::n_distinct(
        cluster
      ),
      
      unique_dominant_topologies =
        dplyr::n_distinct(
          dominant_topology
        ),
      
      topology_coverage =
        unique_dominant_topologies /
        observed_clusters,
      
      duplicate_topology_count =
        observed_clusters -
        unique_dominant_topologies,
      
      duplicate_topology_fraction =
        duplicate_topology_count /
        observed_clusters,
      
      minimum_purity = min(
        purity
      ),
      
      mean_purity = mean(
        purity
      ),
      
      weighted_purity =
        stats::weighted.mean(
          purity,
          w = n_genes
        ),
      
      passing_fraction = mean(
        passes_purity
      ),
      
      maximum_normalised_entropy = max(
        normalised_entropy
      ),
      
      mean_normalised_entropy = mean(
        normalised_entropy
      ),
      
      total_features = sum(
        n_genes
      ),
      
      .groups = "drop"
    )
  
  # ---------------------------------------------------------------------------
  # Identify repeated dominant topology labels
  # ---------------------------------------------------------------------------
  
  topology_repetition <- purity_data |>
    dplyr::count(
      group,
      dominant_topology,
      name = "n_clusters"
    ) |>
    dplyr::arrange(
      group,
      dplyr::desc(n_clusters),
      dominant_topology
    ) |>
    dplyr::mutate(
      is_duplicated_topology =
        n_clusters > 1L
    )
  
  # ---------------------------------------------------------------------------
  # Calculate overall metrics
  # ---------------------------------------------------------------------------
  
  observed_cluster_group_panels <- nrow(
    purity_data
  )
  
  unique_group_topology_pairs <- purity_data |>
    dplyr::distinct(
      group,
      dominant_topology
    ) |>
    nrow()
  
  overall_topology_coverage <-
    unique_group_topology_pairs /
    observed_cluster_group_panels
  
  overall_duplicate_count <-
    observed_cluster_group_panels -
    unique_group_topology_pairs
  
  overall <- data.frame(
    requested_k = if (is.null(k)) {
      NA_integer_
    } else {
      k
    },
    
    observed_clusters = dplyr::n_distinct(
      purity_data$cluster
    ),
    
    n_groups = dplyr::n_distinct(
      purity_data$group
    ),
    
    n_cluster_group_panels =
      observed_cluster_group_panels,
    
    unique_dominant_topologies =
      dplyr::n_distinct(
        purity_data$dominant_topology
      ),
    
    unique_group_topology_pairs =
      unique_group_topology_pairs,
    
    topology_coverage =
      overall_topology_coverage,
    
    duplicate_topology_count =
      overall_duplicate_count,
    
    duplicate_topology_fraction =
      overall_duplicate_count /
      observed_cluster_group_panels,
    
    minimum_purity = min(
      purity_data$purity
    ),
    
    mean_purity = mean(
      purity_data$purity
    ),
    
    weighted_purity =
      stats::weighted.mean(
        purity_data$purity,
        w = purity_data$n_genes
      ),
    
    passing_fraction = mean(
      purity_data$passes_purity
    ),
    
    maximum_normalised_entropy = max(
      purity_data$normalised_entropy
    ),
    
    mean_normalised_entropy = mean(
      purity_data$normalised_entropy
    ),
    
    minimum_cluster_size =
      minimum_cluster_size,
    
    BIC = bic,
    
    AIC = aic,
    
    stringsAsFactors = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # Return structured result
  # ---------------------------------------------------------------------------
  
  out <- list(
    group_summary = group_summary,
    topology_repetition =
      topology_repetition,
    overall = overall
  )
  
  class(out) <- "lot_cluster_topology_score"
  
  out
}