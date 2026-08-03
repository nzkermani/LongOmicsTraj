#' Build the cluster evaluation table
#'
#' Combines molecular-feature-to-cluster relationships with feature-level
#' Omics Trajectory Signature assignments.
#'
#' @param object A `LongOmicsTraj` object containing cluster relationships
#'   and feature-level OTS results.
#'
#' @return A data frame containing one row per clustered feature and its
#'   associated topology label.
#'
#' @keywords internal
.build_cluster_eval_table <- function(
    object
) {
  
  if (!methods::is(object, "LongOmicsTraj")) {
    stop(
      "`object` must be a LongOmicsTraj object.",
      call. = FALSE
    )
  }
  
  cluster_map <- as.data.frame(
    object@relationships
  )
  
  if (
    is.null(object@results$gene_ots) ||
    nrow(object@results$gene_ots) == 0L
  ) {
    stop(
      "`gene_ots` is missing. Run topology mapping on the base object first.",
      call. = FALSE
    )
  }
  
  gene_ots <- as.data.frame(
    object@results$gene_ots
  )
  
  required_relationship_columns <- c(
    "child_id",
    "parent_id"
  )
  
  missing_relationship_columns <- setdiff(
    required_relationship_columns,
    colnames(cluster_map)
  )
  
  if (length(missing_relationship_columns) > 0L) {
    stop(
      "Cluster relationships are missing: ",
      paste(
        missing_relationship_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
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
  
  cluster_map$child_id <- toupper(
    as.character(
      cluster_map$child_id
    )
  )
  
  cluster_map$parent_id <- as.character(
    cluster_map$parent_id
  )
  
  gene_ots$object_id <- toupper(
    as.character(
      gene_ots$object_id
    )
  )
  
  gene_ots$topology_label <- as.character(
    gene_ots$topology_label
  )
  
  evaluation_data <- merge(
    x = cluster_map,
    y = gene_ots,
    by.x = "child_id",
    by.y = "object_id",
    all = FALSE,
    sort = FALSE
  )
  
  if (nrow(evaluation_data) == 0L) {
    stop(
      "No clustered features could be matched to gene-level OTS results.",
      call. = FALSE
    )
  }
  
  names(evaluation_data)[
    names(evaluation_data) == "parent_id"
  ] <- "cluster_id"
  
  evaluation_data
}


#' Compute OTS metrics for molecular-feature clusters
#'
#' Calculates the topology distribution, purity, entropy, dominant topology,
#' cluster size and agreement for each molecular-feature cluster.
#'
#' @param object A `LongOmicsTraj` object containing cluster relationships
#'   and feature-level OTS results.
#'
#' @return A `LongOmicsTraj` object with cluster OTS metrics stored in
#'   `object@results$ots_metrics`.
#'
#' @keywords internal
.compute_ots_metrics <- function(
    object
) {
  
  if (!methods::is(object, "LongOmicsTraj")) {
    stop(
      "`object` must be a LongOmicsTraj object.",
      call. = FALSE
    )
  }
  
  if (
    is.null(object@results$gene_ots) ||
    nrow(object@results$gene_ots) == 0L
  ) {
    stop(
      "`gene_ots` is missing. Run topology mapping on the base object first.",
      call. = FALSE
    )
  }
  
  evaluation_data <- .build_cluster_eval_table(
    object
  )
  
  distribution_data <- evaluation_data |>
    dplyr::filter(
      !is.na(.data$cluster_id),
      !is.na(.data$topology_label)
    ) |>
    dplyr::count(
      .data$cluster_id,
      .data$topology_label,
      name = "n"
    ) |>
    dplyr::group_by(
      .data$cluster_id
    ) |>
    dplyr::mutate(
      freq = .data$n /
        sum(.data$n)
    ) |>
    dplyr::ungroup()
  
  if (nrow(distribution_data) == 0L) {
    stop(
      "No complete cluster and topology assignments were available.",
      call. = FALSE
    )
  }
  
  metrics <- distribution_data |>
    dplyr::group_by(
      .data$cluster_id
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$freq),
      .data$topology_label,
      .by_group = TRUE
    ) |>
    dplyr::summarise(
      purity = max(
        .data$freq
      ),
      entropy = -sum(
        .data$freq *
          log(.data$freq)
      ),
      dominant_topology = dplyr::first(
        .data$topology_label
      ),
      n_genes = sum(
        .data$n
      ),
      agreement = max(
        .data$freq
      ),
      .groups = "drop"
    )
  
  object@results$ots_metrics <- S4Vectors::DataFrame(
    metrics
  )
  
  object
}