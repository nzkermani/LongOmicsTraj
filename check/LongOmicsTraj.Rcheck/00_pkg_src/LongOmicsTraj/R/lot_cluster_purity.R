#' Calculate topology purity within clusters
#'
#' Calculates the dominant Omics Trajectory Signature, topology purity
#' and entropy within each cluster and clinical group.
#'
#' Topology purity is defined as the proportion of features belonging to
#' the most frequent OTS class within a cluster-group combination.
#'
#' @param assignments A data frame containing feature-to-cluster assignments
#'   and feature-level OTS labels. Required columns are `gene_id`, `cluster`,
#'   `group` and `topology_label`.
#' @param purity_threshold Numeric value between 0 and 1 used to classify a
#'   cluster-group combination as sufficiently topology-homogeneous.
#'
#' @return An object of class `lot_cluster_purity` containing:
#' \describe{
#'   \item{cluster_group_purity}{
#'     Purity and entropy statistics for each cluster-group combination.
#'   }
#'   \item{overall}{
#'     Overall minimum, mean and weighted purity, passing fraction and
#'     maximum normalised entropy.
#'   }
#' }
#'
#' @export
#'
#' @examples
#' assignments <- data.frame(
#'   gene_id = paste0("G", 1:10),
#'   cluster = c(rep(1, 5), rep(2, 5)),
#'   group = "Withdrawal",
#'   topology_label = c(
#'     rep("up_down", 4),
#'     "flat_flat",
#'     rep("down_up", 5)
#'   )
#' )
#'
#' result <- lot_cluster_purity(
#'   assignments,
#'   purity_threshold = 0.90
#' )
#'
#' result$cluster_group_purity
#' result$overall
lot_cluster_purity <- function(
    assignments,
    purity_threshold = 0.90
) {
  
  # ---------------------------------------------------------------------------
  # Validate inputs
  # ---------------------------------------------------------------------------
  
  if (
    !is.data.frame(assignments) &&
    !methods::is(assignments, "DataFrame")
  ) {
    stop(
      "`assignments` must be a data frame or S4Vectors DataFrame.",
      call. = FALSE
    )
  }
  
  assignments <- as.data.frame(
    assignments
  )
  
  required_columns <- c(
    "gene_id",
    "cluster",
    "group",
    "topology_label"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(assignments)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Assignments are missing: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
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
      "`purity_threshold` must be one finite number between 0 and 1.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Remove incomplete assignments
  # ---------------------------------------------------------------------------
  
  assignments_complete <- assignments |>
    dplyr::filter(
      !is.na(gene_id),
      !is.na(cluster),
      !is.na(group),
      !is.na(topology_label)
    ) |>
    dplyr::mutate(
      gene_id = as.character(gene_id),
      group = as.character(group),
      topology_label = as.character(topology_label)
    )
  
  if (nrow(assignments_complete) == 0L) {
    stop(
      "No complete cluster and topology assignments were available.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Ensure one topology label per feature, cluster and group
  # ---------------------------------------------------------------------------
  
  conflicting_assignments <- assignments_complete |>
    dplyr::distinct(
      gene_id,
      cluster,
      group,
      topology_label
    ) |>
    dplyr::count(
      gene_id,
      cluster,
      group,
      name = "n_topology_labels"
    ) |>
    dplyr::filter(
      n_topology_labels > 1L
    )
  
  if (nrow(conflicting_assignments) > 0L) {
    stop(
      "Some features have more than one topology label within the same ",
      "cluster and group.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Count topology classes
  # ---------------------------------------------------------------------------
  
  topology_counts <- assignments_complete |>
    dplyr::distinct(
      gene_id,
      cluster,
      group,
      topology_label
    ) |>
    dplyr::count(
      cluster,
      group,
      topology_label,
      name = "topology_count"
    )
  
  # ---------------------------------------------------------------------------
  # Calculate purity and entropy for each cluster-group combination
  # ---------------------------------------------------------------------------
  
  cluster_group_purity <- topology_counts |>
    dplyr::group_by(
      cluster,
      group
    ) |>
    dplyr::mutate(
      total_genes = sum(
        topology_count
      ),
      proportion = topology_count /
        total_genes
    ) |>
    dplyr::arrange(
      dplyr::desc(proportion),
      topology_label,
      .by_group = TRUE
    ) |>
    dplyr::summarise(
      n_genes = dplyr::first(
        total_genes
      ),
      dominant_topology = dplyr::first(
        topology_label
      ),
      dominant_count = dplyr::first(
        topology_count
      ),
      purity = dplyr::first(
        proportion
      ),
      entropy = -sum(
        proportion * log(proportion)
      ),
      n_topologies = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      normalised_entropy = dplyr::if_else(
        n_topologies > 1L,
        entropy / log(n_topologies),
        0
      ),
      passes_purity =
        purity >= purity_threshold
    )
  
  # ---------------------------------------------------------------------------
  # Calculate overall summary statistics
  # ---------------------------------------------------------------------------
  
  overall <- cluster_group_purity |>
    dplyr::summarise(
      minimum_purity = min(
        purity
      ),
      mean_purity = mean(
        purity
      ),
      weighted_purity = stats::weighted.mean(
        purity,
        w = n_genes
      ),
      passing_fraction = mean(
        passes_purity
      ),
      maximum_normalised_entropy = max(
        normalised_entropy
      ),
      n_cluster_group_combinations =
        dplyr::n(),
      .groups = "drop"
    )
  
  # ---------------------------------------------------------------------------
  # Return structured result
  # ---------------------------------------------------------------------------
  
  out <- list(
    cluster_group_purity =
      cluster_group_purity,
    overall = overall,
    purity_threshold =
      purity_threshold
  )
  
  class(out) <- "lot_cluster_purity"
  
  out
}