plot_cluster_spaghetti_single <- function(object, cluster_id, assay_name, visits_vector) {

  library(ggplot2)
  library(dplyr)
  library(tidyr)

  # ----------------------------------------------------------------------------
  # 1. Relationships
  # ----------------------------------------------------------------------------
  rel_df <- as.data.frame(object@relationships)
  rel_df <- rel_df[rel_df$assay == assay_name, ]
  rel_df$child_id <- toupper(rel_df$child_id)

  # ----------------------------------------------------------------------------
  # 2. Expression + metadata
  # ----------------------------------------------------------------------------
  se <- MultiAssayExperiment::experiments(object@mae)[[assay_name]]
  raw_mat <- SummarizedExperiment::assay(se)
  meta_df <- as.data.frame(SummarizedExperiment::colData(se))

  rownames(raw_mat) <- toupper(rownames(raw_mat))

  meta_df$sample_id <- rownames(meta_df)
  meta_df <- meta_df[meta_df$visit %in% visits_vector, ]

  groups <- unique(meta_df$group)

  # ----------------------------------------------------------------------------
  # 3. Trajectory matrix
  # ----------------------------------------------------------------------------
  traj_list <- lapply(groups, function(g) {
    meta_g <- meta_df[meta_df$group == g, ]

    sapply(visits_vector, function(v) {
      samples <- meta_g$sample_id[meta_g$visit == v]
      if (length(samples) == 0) return(rep(NA, nrow(raw_mat)))
      rowMeans(raw_mat[, samples, drop = FALSE], na.rm = TRUE)
    })
  })

  traj_mat <- do.call(cbind, traj_list)
  colnames(traj_mat) <- as.vector(outer(groups, visits_vector, paste, sep = "_"))

  # Filter genes in this cluster only
  genes_in_cluster <- rel_df %>%
    filter(parent_id == cluster_id) %>%
    pull(child_id)

  traj_mat <- traj_mat[rownames(traj_mat) %in% genes_in_cluster, , drop = FALSE]

  # Scale
  traj_scaled <- t(scale(t(traj_mat)))

  # ----------------------------------------------------------------------------
  # 4. Long format
  # ----------------------------------------------------------------------------
  gene_long <- as.data.frame(traj_scaled) %>%
    mutate(gene_id = rownames(traj_scaled)) %>%
    pivot_longer(cols = -gene_id, names_to = "combined", values_to = "expression") %>%
    separate(combined, into = c("group", "visit"), sep = "_") %>%
    mutate(cluster = cluster_id)

  gene_long$visit <- factor(gene_long$visit, levels = visits_vector)

  # ----------------------------------------------------------------------------
  # 5. Add topology
  # ----------------------------------------------------------------------------
  gene_ots <- as.data.frame(object@results$gene_ots)
  gene_ots$object_id <- toupper(gene_ots$object_id)

  gene_long <- gene_long %>%
    left_join(gene_ots[, c("object_id", "topology_label")],
              by = c("gene_id" = "object_id"))

  # ----------------------------------------------------------------------------
  # 6. Cluster center
  # ----------------------------------------------------------------------------
  vmm_df <- as.data.frame(object@vmm) %>%
    filter(assay == assay_name, object_id == cluster_id) %>%
    rename(cluster = object_id,
           visit = visit,
           expression = estimated_value)

  vmm_df$visit <- factor(vmm_df$visit, levels = visits_vector)

  # ----------------------------------------------------------------------------
  # 7. Topology % (legend labels)
  # ----------------------------------------------------------------------------
  topo_dist <- gene_long %>%
    count(topology_label) %>%
    mutate(
      pct = round(100 * n / sum(n), 1),
      label = paste0(topology_label, " (", pct, "%)")
    )

  gene_long <- gene_long %>%
    left_join(topo_dist[, c("topology_label", "label")],
              by = "topology_label")

  # ----------------------------------------------------------------------------
  # 8. Metrics
  # ----------------------------------------------------------------------------
  metrics <- as.data.frame(object@results$ots_metrics)
  cl_metrics <- metrics %>% filter(cluster_id == cluster_id)

  purity_val  <- round(cl_metrics$purity[1], 2)
  entropy_val <- round(cl_metrics$entropy[1], 2)

  # ----------------------------------------------------------------------------
  # 9. Plot
  # ----------------------------------------------------------------------------
  p <- ggplot() +

    geom_line(
      data = gene_long,
      aes(x = visit, y = expression, group = gene_id, color = label),
      alpha = 0.4,
      linewidth = 0.5
    ) +

    geom_line(
      data = vmm_df,
      aes(x = visit, y = expression, group = cluster),
      color = "black",
      linewidth = 1.5
    ) +

    geom_point(
      data = vmm_df,
      aes(x = visit, y = expression),
      color = "black",
      size = 3
    ) +

    scale_color_brewer(palette = "Set2", na.value = "grey70") +

    labs(
      title = paste("Cluster:", cluster_id),
      subtitle = paste0(
        "n = ", length(unique(gene_long$gene_id)),
        " | Purity = ", purity_val,
        " | Entropy = ", entropy_val
      ),
      x = "Visit",
      y = "Scaled expression",
      color = "Topology (%)"
    ) +

    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    )

  return(p)
}
