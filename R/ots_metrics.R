.build_cluster_eval_table <- function(object) {

  cluster_map <- as.data.frame(object@relationships)
  gene_ots    <- as.data.frame(object@results$gene_ots)

  # Standardise once
  cluster_map$child_id <- toupper(cluster_map$child_id)
  gene_ots$object_id   <- toupper(gene_ots$object_id)

  # Build canonical table
  df <- merge(
    cluster_map,
    gene_ots,
    by.x = "child_id",
    by.y = "object_id"
  )

  # Rename for clarity (important)
  colnames(df)[colnames(df) == "parent_id"] <- "cluster_id"

  df
}
.compute_ots_metrics <- function(object) {

  if (is.null(object@results$gene_ots)) {
    stop("gene_ots missing — run topology on base object first")
  }

  df <- .build_cluster_eval_table(object)

  library(dplyr)

  # Distribution
  dist_df <- df %>%
    group_by(cluster_id, topology_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(cluster_id) %>%
    mutate(freq = n / sum(n))

  # Metrics
  metrics <- dist_df %>%
    summarise(
      purity = max(freq),
      entropy = -sum(freq * log(freq)),
      dominant_topology = topology_label[which.max(freq)],
      n_genes = sum(n),
      agreement = max(freq),   # ✅ THIS IS THE CORRECT DEFINITION
      .groups = "drop"
    )

  object@results$ots_metrics <- S4Vectors::DataFrame(metrics)

  object
}
