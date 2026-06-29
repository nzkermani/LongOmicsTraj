#' @export
setMethod(
  "lot_from_clusters",
  "LongOmicsTraj",
  function(
    object,
    assay,
    group_col = "group",
    visit_col = "visit",
    visits,
    k = 4,
    estimator = c("kmeans", "pam"),
    overwrite = FALSE,
    verbose = TRUE
  ) {

    estimator <- match.arg(estimator)

    # ------------------------------------------------------------
    # Checks
    # ------------------------------------------------------------
    if (!assay %in% names(MultiAssayExperiment::experiments(object@mae))) {
      stop("Assay not found", call. = FALSE)
    }

    if (nrow(object@vmm) > 0 && !overwrite) {
      stop("VMM exists. Use overwrite=TRUE.", call. = FALSE)
    }

    # ------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------
    se <- MultiAssayExperiment::experiments(object@mae)[[assay]]
    mat <- SummarizedExperiment::assay(se)
    meta <- as.data.frame(SummarizedExperiment::colData(se))

    mat <- mat[, rownames(meta), drop = FALSE]

    # ------------------------------------------------------------
    # Metadata
    # ------------------------------------------------------------
    meta$sample_id <- rownames(meta)
    meta <- meta[meta[[visit_col]] %in% visits, , drop = FALSE]
    meta[[visit_col]] <- factor(meta[[visit_col]], levels = visits, ordered = TRUE)

    groups <- unique(meta[[group_col]])

    # ------------------------------------------------------------
    # Build trajectory matrix
    # ------------------------------------------------------------
    traj_list <- lapply(groups, function(g) {
      meta_g <- meta[meta[[group_col]] == g, , drop = FALSE]
      sapply(visits, function(v) {
        samples <- meta_g$sample_id[meta_g[[visit_col]] == v]
        if (length(samples) == 0) return(rep(NA, nrow(mat)))
        rowMeans(mat[, samples, drop = FALSE], na.rm = TRUE)
      })
    })

    traj_mat <- do.call(cbind, traj_list)

    # Clean IDs
    rownames(traj_mat) <- toupper(rownames(traj_mat))
    center_cols <- as.vector(outer(groups, visits, paste, sep = "_"))
    colnames(traj_mat) <- center_cols

    # ------------------------------------------------------------
    # Remove low variance
    # ------------------------------------------------------------
    keep <- matrixStats::rowSds(traj_mat, na.rm = TRUE) > 0
    traj_mat <- traj_mat[keep, , drop = FALSE]

    if (nrow(traj_mat) == 0) {
      stop("No variable genes available.", call. = FALSE)
    }

    # ------------------------------------------------------------
    # Scale
    # ------------------------------------------------------------
    traj_scaled <- t(scale(t(traj_mat)))

    # ------------------------------------------------------------
    # Clustering
    # ------------------------------------------------------------
    if (estimator == "kmeans") {

      set.seed(1)
      km <- kmeans(traj_scaled, centers = k)

      cluster_assignments <- as.integer(factor(km$cluster))
      names(cluster_assignments) <- rownames(traj_scaled)

    } else if (estimator == "pam") {

      if (!requireNamespace("cluster", quietly = TRUE)) {
        stop("Install 'cluster' package for PAM.", call. = FALSE)
      }

      pam_fit <- cluster::pam(traj_scaled, k)

      cluster_assignments <- as.integer(factor(pam_fit$clustering))
      names(cluster_assignments) <- rownames(traj_scaled)
    }

    # ------------------------------------------------------------
    # Compute centers
    # ------------------------------------------------------------
    clusters <- sort(unique(cluster_assignments))

    centers <- do.call(rbind, lapply(clusters, function(cl) {
      members <- traj_scaled[cluster_assignments == cl, , drop = FALSE]
      if (nrow(members) == 0) return(rep(NA, ncol(traj_scaled)))
      colMeans(members, na.rm = TRUE)
    }))

    rownames(centers) <- paste0("cluster_", clusters)

    # ------------------------------------------------------------
    # Build VMM
    # ------------------------------------------------------------
    parsed <- do.call(rbind, strsplit(center_cols, "_"))

    col_df <- data.frame(
      group = parsed[, 1],
      visit = parsed[, 2],
      stringsAsFactors = FALSE
    )

    vmm_list <- lapply(seq_len(nrow(centers)), function(i) {

      df <- col_df
      df$estimated_value <- as.numeric(centers[i, ])

      df$object_id <- rownames(centers)[i]
      df$object_type <- "cluster"
      df$assay <- assay
      df$estimator <- estimator

      df
    })

    object@vmm <- S4Vectors::DataFrame(do.call(rbind, vmm_list))

    # ------------------------------------------------------------
    # Relationships
    # ------------------------------------------------------------
    cluster_map <- data.frame(
      child_id = rownames(traj_mat),
      parent_id = paste0("cluster_", cluster_assignments[rownames(traj_mat)]),
      relationship = "cluster_membership",
      weight = 1,
      assay = assay,
      method = estimator,
      stringsAsFactors = FALSE
    )

    new_rel <- S4Vectors::DataFrame(cluster_map)

    if (nrow(object@relationships) == 0) {
      object@relationships <- new_rel
    } else {
      object@relationships <- S4Vectors::rbind(object@relationships, new_rel)
    }
    # ------------------------------------------------------------
    # Inject gene-level OTS (CRITICAL FOR METRICS)
    # ------------------------------------------------------------

    if (!is.null(object@ots) && nrow(object@ots) > 0) {

      gene_ots <- as.data.frame(object@ots)

      # Standardise IDs (match clustering)
      gene_ots$object_id <- toupper(gene_ots$object_id)

      object@results$gene_ots <- S4Vectors::DataFrame(gene_ots)

    } else {
      warning("No gene-level OTS found in object. Metrics may fail.")
    }

    # ------------------------------------------------------------
    # History
    # ------------------------------------------------------------
    object@history <- c(
      object@history,
      list(list(
        step = "lot_from_clusters",
        assay = assay,
        estimator = estimator,
        n_clusters = k,
        n_features = nrow(traj_mat),
        time = Sys.time()
      ))
    )

    if (verbose) {
      message(
        "Created cluster-based VMM with ", k,
        " clusters using '", estimator,
        "' on ", nrow(traj_mat), " genes."
      )
    }

    validObject(object)
    object
  }
)
