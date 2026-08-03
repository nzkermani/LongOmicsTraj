#' Cluster molecular features by longitudinal trajectory
#'
#' Clusters molecular features using their combined group-by-visit
#' longitudinal profiles. Cluster centres are stored as visit-level
#' molecular models, while feature-to-cluster memberships and clustering
#' diagnostics are stored in the results and relationships slots.
#'
#' When feature-level OTS results already exist, they are preserved and
#' summarised within each trajectory cluster.
#'
#' @param object A `LongOmicsTraj` object.
#' @param assay Name of the assay to analyse.
#' @param group_col Metadata column defining clinical or treatment groups.
#' @param visit_col Metadata column defining visits.
#' @param visits Character vector giving the chronological visit order.
#' @param k Number of clusters.
#' @param estimator Clustering method: `"kmeans"` or `"pam"`.
#' @param assay_data Optional assay-matrix name. The first available matrix
#'   is used by default.
#' @param scale_features Logical; standardise each feature trajectory before
#'   clustering. The default is `TRUE`.
#' @param standardise_ids Logical; convert feature identifiers to uppercase
#'   before clustering and matching OTS results. The default is `FALSE`.
#' @param seed Random seed used by k-means.
#' @param require_gene_ots Logical; require existing feature-level OTS results.
#' @param overwrite Logical; overwrite an existing visit-level model.
#' @param verbose Logical; print progress messages.
#'
#' @return A `LongOmicsTraj` object containing cluster-centre trajectories,
#'   feature-to-cluster relationships, raw and scaled cluster profiles,
#'   clustering diagnostics and, when available, feature-level OTS
#'   distributions within clusters.
#'
#' @details
#' Each molecular feature is represented by a vector containing its mean
#' value for every observed group-by-visit combination. Clustering is
#' therefore based on the feature's combined longitudinal behaviour across
#' all included groups.
#'
#' When `scale_features = TRUE`, clustering and cluster-level OTS mapping
#' use standardised feature profiles. Both raw and scaled cluster centres
#' are retained in `object@results$cluster_profiles`.
#' @rdname lot_from_clusters
#' @export
methods::setMethod(
  "lot_from_clusters",
  signature = "LongOmicsTraj",
  function(
    object,
    assay,
    group_col = "group",
    visit_col = "visit",
    visits,
    k = 4L,
    estimator = c("kmeans", "pam"),
    assay_data = NULL,
    scale_features = TRUE,
    standardise_ids = FALSE,
    seed = 1L,
    require_gene_ots = FALSE,
    overwrite = FALSE,
    verbose = TRUE
  ) {
    
    estimator <- match.arg(estimator)
    
    # -------------------------------------------------------------------------
    # Validate arguments
    # -------------------------------------------------------------------------
    
    if (
      missing(assay) ||
      !is.character(assay) ||
      length(assay) != 1L ||
      is.na(assay) ||
      !nzchar(assay)
    ) {
      stop(
        "`assay` must be one non-empty character string.",
        call. = FALSE
      )
    }
    
    if (
      missing(visits) ||
      !is.character(visits) ||
      length(visits) < 2L ||
      anyNA(visits) ||
      any(!nzchar(visits)) ||
      anyDuplicated(visits)
    ) {
      stop(
        "`visits` must contain at least two unique, non-missing visit labels.",
        call. = FALSE
      )
    }
    
    if (
      length(k) != 1L ||
      is.na(k) ||
      !is.numeric(k) ||
      !is.finite(k) ||
      k < 2L ||
      k != as.integer(k)
    ) {
      stop(
        "`k` must be one integer greater than or equal to 2.",
        call. = FALSE
      )
    }
    
    if (
      length(seed) != 1L ||
      is.na(seed) ||
      !is.numeric(seed) ||
      !is.finite(seed)
    ) {
      stop(
        "`seed` must be one finite numeric value.",
        call. = FALSE
      )
    }
    
    if (
      !is.logical(scale_features) ||
      length(scale_features) != 1L ||
      is.na(scale_features)
    ) {
      stop(
        "`scale_features` must be `TRUE` or `FALSE`.",
        call. = FALSE
      )
    }
    
    if (
      !is.logical(standardise_ids) ||
      length(standardise_ids) != 1L ||
      is.na(standardise_ids)
    ) {
      stop(
        "`standardise_ids` must be `TRUE` or `FALSE`.",
        call. = FALSE
      )
    }
    
    if (
      !is.logical(require_gene_ots) ||
      length(require_gene_ots) != 1L ||
      is.na(require_gene_ots)
    ) {
      stop(
        "`require_gene_ots` must be `TRUE` or `FALSE`.",
        call. = FALSE
      )
    }
    
    k <- as.integer(k)
    seed <- as.integer(seed)
    
    experiments <- MultiAssayExperiment::experiments(
      object@mae
    )
    
    experiment_names <- names(experiments)
    
    if (!assay %in% experiment_names) {
      stop(
        "Assay `", assay, "` was not found. Available assays: ",
        paste(experiment_names, collapse = ", "),
        call. = FALSE
      )
    }
    
    if (nrow(object@vmm) > 0L && !isTRUE(overwrite)) {
      stop(
        "Visit-level trajectories already exist. ",
        "Use `overwrite = TRUE` to replace them.",
        call. = FALSE
      )
    }
    
    # -------------------------------------------------------------------------
    # Extract molecular data and metadata
    # -------------------------------------------------------------------------
    
    se <- experiments[[assay]]
    
    available_matrices <- SummarizedExperiment::assayNames(
      se
    )
    
    if (length(available_matrices) == 0L) {
      stop(
        "The selected assay contains no data matrices.",
        call. = FALSE
      )
    }
    
    if (is.null(assay_data)) {
      assay_data <- available_matrices[[1L]]
    }
    
    if (
      !is.character(assay_data) ||
      length(assay_data) != 1L ||
      is.na(assay_data) ||
      !assay_data %in% available_matrices
    ) {
      stop(
        "Data matrix `", assay_data, "` was not found. Available matrices: ",
        paste(available_matrices, collapse = ", "),
        call. = FALSE
      )
    }
    
    mat <- SummarizedExperiment::assay(
      se,
      assay_data
    )
    
    if (!is.numeric(mat)) {
      stop(
        "The selected assay matrix must be numeric.",
        call. = FALSE
      )
    }
    
    metadata <- as.data.frame(
      SummarizedExperiment::colData(se)
    )
    
    required_metadata <- unique(
      c(
        "sample_id",
        group_col,
        visit_col
      )
    )
    
    missing_metadata <- setdiff(
      required_metadata,
      colnames(metadata)
    )
    
    if (length(missing_metadata) > 0L) {
      stop(
        "Metadata are missing: ",
        paste(missing_metadata, collapse = ", "),
        call. = FALSE
      )
    }
    
    metadata$sample_id <- as.character(
      metadata$sample_id
    )
    
    if (
      anyNA(metadata$sample_id) ||
      any(!nzchar(metadata$sample_id)) ||
      anyDuplicated(metadata$sample_id)
    ) {
      stop(
        "`sample_id` must contain unique, non-missing identifiers.",
        call. = FALSE
      )
    }
    
    metadata <- metadata[
      match(
        colnames(mat),
        metadata$sample_id
      ),
      ,
      drop = FALSE
    ]
    
    if (
      anyNA(metadata$sample_id) ||
      !identical(
        colnames(mat),
        metadata$sample_id
      )
    ) {
      stop(
        "Expression samples could not be aligned with the metadata.",
        call. = FALSE
      )
    }
    
    metadata <- metadata[
      as.character(metadata[[visit_col]]) %in% visits &
        !is.na(metadata[[group_col]]) &
        nzchar(as.character(metadata[[group_col]])),
      ,
      drop = FALSE
    ]
    
    if (nrow(metadata) == 0L) {
      stop(
        "No samples remained after filtering the requested visits.",
        call. = FALSE
      )
    }
    
    metadata[[visit_col]] <- factor(
      as.character(metadata[[visit_col]]),
      levels = visits,
      ordered = TRUE
    )
    
    mat <- mat[
      ,
      metadata$sample_id,
      drop = FALSE
    ]
    
    groups <- unique(
      as.character(metadata[[group_col]])
    )
    
    groups <- groups[
      !is.na(groups) &
        nzchar(groups)
    ]
    
    if (length(groups) == 0L) {
      stop(
        "No valid groups were found.",
        call. = FALSE
      )
    }
    
    # -------------------------------------------------------------------------
    # Prepare feature identifiers
    # -------------------------------------------------------------------------
    
    original_feature_ids <- rownames(mat)
    
    if (
      is.null(original_feature_ids) ||
      anyNA(original_feature_ids) ||
      any(!nzchar(original_feature_ids))
    ) {
      stop(
        "The assay matrix must have valid feature row names.",
        call. = FALSE
      )
    }
    
    feature_ids <- if (isTRUE(standardise_ids)) {
      toupper(original_feature_ids)
    } else {
      original_feature_ids
    }
    
    if (anyDuplicated(feature_ids)) {
      
      duplicated_ids <- unique(
        feature_ids[
          duplicated(feature_ids)
        ]
      )
      
      stop(
        "Feature identifiers are duplicated",
        if (isTRUE(standardise_ids)) {
          " after converting to uppercase"
        } else {
          ""
        },
        ": ",
        paste(
          utils::head(
            duplicated_ids,
            10L
          ),
          collapse = ", "
        ),
        if (length(duplicated_ids) > 10L) {
          " ..."
        } else {
          ""
        },
        call. = FALSE
      )
    }
    
    rownames(mat) <- feature_ids
    
    # -------------------------------------------------------------------------
    # Construct the complete group-by-visit design
    # -------------------------------------------------------------------------
    
    design_grid <- expand.grid(
      group = groups,
      visit = visits,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    
    design_counts <- vapply(
      seq_len(nrow(design_grid)),
      function(i) {
        
        sum(
          as.character(metadata[[group_col]]) ==
            design_grid$group[[i]] &
            as.character(metadata[[visit_col]]) ==
            design_grid$visit[[i]],
          na.rm = TRUE
        )
      },
      integer(1)
    )
    
    design_grid$n_samples <- design_counts
    
    missing_profiles <- design_grid[
      design_grid$n_samples == 0L,
      ,
      drop = FALSE
    ]
    
    if (
      nrow(missing_profiles) > 0L &&
      isTRUE(verbose)
    ) {
      message(
        "Omitting ",
        nrow(missing_profiles),
        " group-by-visit combination",
        if (nrow(missing_profiles) == 1L) "" else "s",
        " with no samples."
      )
    }
    
    design_grid <- design_grid[
      design_grid$n_samples > 0L,
      ,
      drop = FALSE
    ]
    
    if (nrow(design_grid) < 2L) {
      stop(
        "Fewer than two observed group-by-visit profiles were available.",
        call. = FALSE
      )
    }
    
    # -------------------------------------------------------------------------
    # Build the feature trajectory matrix
    # -------------------------------------------------------------------------
    
    trajectory_columns <- lapply(
      seq_len(nrow(design_grid)),
      function(i) {
        
        current_group <- design_grid$group[[i]]
        current_visit <- design_grid$visit[[i]]
        
        sample_ids <- metadata$sample_id[
          as.character(metadata[[group_col]]) == current_group &
            as.character(metadata[[visit_col]]) == current_visit
        ]
        
        rowMeans(
          mat[
            ,
            sample_ids,
            drop = FALSE
          ],
          na.rm = TRUE
        )
      }
    )
    
    traj_mat <- do.call(
      cbind,
      trajectory_columns
    )
    
    rownames(traj_mat) <- feature_ids
    
    colnames(traj_mat) <- paste0(
      "profile_",
      seq_len(ncol(traj_mat))
    )
    
    # -------------------------------------------------------------------------
    # Remove incomplete or invariant features
    # -------------------------------------------------------------------------
    
    finite_counts <- rowSums(
      is.finite(traj_mat)
    )
    
    feature_sd <- matrixStats::rowSds(
      traj_mat,
      na.rm = TRUE
    )
    
    keep <- finite_counts == ncol(traj_mat) &
      is.finite(feature_sd) &
      feature_sd > 0
    
    removed_features <- sum(!keep)
    
    traj_mat <- traj_mat[
      keep,
      ,
      drop = FALSE
    ]
    
    if (
      removed_features > 0L &&
      isTRUE(verbose)
    ) {
      message(
        "Removed ",
        removed_features,
        " incomplete or invariant feature",
        if (removed_features == 1L) "" else "s",
        "."
      )
    }
    
    if (nrow(traj_mat) == 0L) {
      stop(
        "No complete, variable features were available for clustering.",
        call. = FALSE
      )
    }
    
    if (k > nrow(traj_mat)) {
      stop(
        "`k` cannot exceed the number of retained features (",
        nrow(traj_mat),
        ").",
        call. = FALSE
      )
    }
    
    # -------------------------------------------------------------------------
    # Scale feature trajectories
    # -------------------------------------------------------------------------
    
    if (isTRUE(scale_features)) {
      
      traj_for_clustering <- t(
        base::scale(
          t(traj_mat),
          center = TRUE,
          scale = TRUE
        )
      )
      
    } else {
      
      traj_for_clustering <- traj_mat
    }
    
    if (any(!is.finite(traj_for_clustering))) {
      stop(
        "The clustering matrix contains non-finite values after scaling.",
        call. = FALSE
      )
    }
    
    # -------------------------------------------------------------------------
    # Cluster molecular features
    # -------------------------------------------------------------------------
    
    if (estimator == "kmeans") {
      
      fit <- withr::with_seed(
        seed,
        stats::kmeans(
          x = traj_for_clustering,
          centers = k,
          nstart = 25
        )
      )
      
      cluster_assignments <- as.integer(
        fit$cluster
      )
      
    } else {
      
      if (!requireNamespace(
        "cluster",
        quietly = TRUE
      )) {
        stop(
          "Package `cluster` is required for PAM.",
          call. = FALSE
        )
      }
      
      fit <- cluster::pam(
        x = traj_for_clustering,
        k = k
      )
      
      cluster_assignments <- as.integer(
        fit$clustering
      )
    }
    
    names(cluster_assignments) <- rownames(
      traj_for_clustering
    )
    
    cluster_levels <- sort(
      unique(cluster_assignments)
    )
    
    cluster_ids <- paste0(
      "cluster_",
      cluster_levels
    )
    
    # -------------------------------------------------------------------------
    # Calculate scaled and raw cluster centres
    # -------------------------------------------------------------------------
    
    scaled_centers <- do.call(
      rbind,
      lapply(
        cluster_levels,
        function(cluster_number) {
          
          members <- traj_for_clustering[
            cluster_assignments == cluster_number,
            ,
            drop = FALSE
          ]
          
          colMeans(
            members,
            na.rm = TRUE
          )
        }
      )
    )
    
    raw_centers <- do.call(
      rbind,
      lapply(
        cluster_levels,
        function(cluster_number) {
          
          members <- traj_mat[
            cluster_assignments == cluster_number,
            ,
            drop = FALSE
          ]
          
          colMeans(
            members,
            na.rm = TRUE
          )
        }
      )
    )
    
    rownames(scaled_centers) <- cluster_ids
    rownames(raw_centers) <- cluster_ids
    
    colnames(scaled_centers) <- colnames(traj_mat)
    colnames(raw_centers) <- colnames(traj_mat)
    
    # -------------------------------------------------------------------------
    # Store cluster-centre trajectories as VMM
    # -------------------------------------------------------------------------
    
    centres_for_vmm <- if (isTRUE(scale_features)) {
      scaled_centers
    } else {
      raw_centers
    }
    
    vmm_value_type <- if (isTRUE(scale_features)) {
      "scaled"
    } else {
      "raw"
    }
    
    vmm_list <- lapply(
      seq_len(nrow(centres_for_vmm)),
      function(i) {
        
        data.frame(
          object_id = rownames(centres_for_vmm)[[i]],
          object_type = "cluster",
          assay = assay,
          group = design_grid$group,
          visit = factor(
            design_grid$visit,
            levels = visits,
            ordered = TRUE
          ),
          estimator = estimator,
          estimated_value = as.numeric(
            centres_for_vmm[i, ]
          ),
          value_type = vmm_value_type,
          stringsAsFactors = FALSE
        )
      }
    )
    
    object@vmm <- S4Vectors::DataFrame(
      do.call(
        rbind,
        vmm_list
      )
    )
    
    # -------------------------------------------------------------------------
    # Store feature-to-cluster relationships
    # -------------------------------------------------------------------------
    
    cluster_map <- data.frame(
      child_id = names(cluster_assignments),
      parent_id = paste0(
        "cluster_",
        cluster_assignments
      ),
      relationship = "cluster_membership",
      weight = 1,
      assay = assay,
      method = estimator,
      stringsAsFactors = FALSE
    )
    
    existing_relationships <- as.data.frame(
      object@relationships
    )
    
    if (nrow(existing_relationships) > 0L) {
      
      required_relationship_columns <- c(
        "child_id",
        "parent_id",
        "relationship",
        "weight",
        "assay",
        "method"
      )
      
      if (
        !all(
          required_relationship_columns %in%
          colnames(existing_relationships)
        )
      ) {
        stop(
          "Existing relationships do not use the expected schema: ",
          paste(
            required_relationship_columns,
            collapse = ", "
          ),
          call. = FALSE
        )
      }
      
      existing_relationships <- existing_relationships[
        !(
          existing_relationships$relationship ==
            "cluster_membership" &
            existing_relationships$assay == assay &
            existing_relationships$method == estimator
        ),
        required_relationship_columns,
        drop = FALSE
      ]
      
      cluster_map <- cluster_map[
        ,
        required_relationship_columns,
        drop = FALSE
      ]
      
      object@relationships <- S4Vectors::DataFrame(
        rbind(
          existing_relationships,
          cluster_map
        )
      )
      
    } else {
      
      object@relationships <- S4Vectors::DataFrame(
        cluster_map
      )
    }
    
    # -------------------------------------------------------------------------
    # Preserve and summarise existing feature-level OTS
    # -------------------------------------------------------------------------
    
    gene_ots_available <- nrow(
      as.data.frame(object@ots)
    ) > 0L
    
    if (
      isTRUE(require_gene_ots) &&
      !gene_ots_available
    ) {
      stop(
        "Feature-level OTS results are required. ",
        "Run the feature-level trajectory pipeline before clustering.",
        call. = FALSE
      )
    }
    
    if (gene_ots_available) {
      
      gene_ots <- as.data.frame(
        object@ots
      )
      
      required_ots_columns <- c(
        "object_id",
        "group",
        "topology_label"
      )
      
      missing_ots_columns <- setdiff(
        required_ots_columns,
        colnames(gene_ots)
      )
      
      if (length(missing_ots_columns) > 0L) {
        stop(
          "Feature-level OTS results are missing: ",
          paste(
            missing_ots_columns,
            collapse = ", "
          ),
          call. = FALSE
        )
      }
      
      gene_ots$object_id <- as.character(
        gene_ots$object_id
      )
      
      if (isTRUE(standardise_ids)) {
        gene_ots$object_id <- toupper(
          gene_ots$object_id
        )
      }
      
      gene_ots <- gene_ots[
        gene_ots$object_id %in%
          cluster_map$child_id,
        ,
        drop = FALSE
      ]
      
      gene_ots <- merge(
        gene_ots,
        cluster_map[
          ,
          c(
            "child_id",
            "parent_id"
          ),
          drop = FALSE
        ],
        by.x = "object_id",
        by.y = "child_id",
        all = FALSE,
        sort = FALSE
      )
      
      names(gene_ots)[
        names(gene_ots) == "parent_id"
      ] <- "cluster_id"
      
      cluster_ots_distribution <- gene_ots |>
        dplyr::filter(
          !is.na(.data$group),
          !is.na(.data$topology_label)
        ) |>
        dplyr::count(
          .data$cluster_id,
          .data$group,
          .data$topology_label,
          name = "count"
        ) |>
        dplyr::group_by(
          .data$cluster_id,
          .data$group
        ) |>
        dplyr::mutate(
          proportion = .data$count /
            sum(.data$count)
        ) |>
        dplyr::ungroup() |>
        dplyr::arrange(
          .data$cluster_id,
          .data$group,
          dplyr::desc(.data$count),
          .data$topology_label
        )
      
      object@results$gene_ots <-
        S4Vectors::DataFrame(
          gene_ots
        )
      
      object@results$cluster_ots_distribution <-
        S4Vectors::DataFrame(
          cluster_ots_distribution
        )
      
    } else if (isTRUE(verbose)) {
      
      message(
        "No feature-level OTS results were present. ",
        "Cluster centres and memberships were created, but ",
        "within-cluster OTS distributions were not calculated."
      )
    }
    
    # -------------------------------------------------------------------------
    # Store cluster memberships and sizes
    # -------------------------------------------------------------------------
    
    cluster_table <- table(
      cluster_assignments
    )
    
    cluster_sizes <- data.frame(
      cluster_id = paste0(
        "cluster_",
        names(cluster_table)
      ),
      n_features = as.integer(
        cluster_table
      ),
      stringsAsFactors = FALSE
    )
    
    object@results$cluster_membership <-
      S4Vectors::DataFrame(
        cluster_map
      )
    
    object@results$cluster_sizes <-
      S4Vectors::DataFrame(
        cluster_sizes
      )
    
    # -------------------------------------------------------------------------
    # Store raw and scaled cluster profiles
    # -------------------------------------------------------------------------
    
    scaled_profiles <- data.frame(
      cluster_id = rep(
        rownames(scaled_centers),
        each = ncol(scaled_centers)
      ),
      group = rep(
        design_grid$group,
        times = nrow(scaled_centers)
      ),
      visit = factor(
        rep(
          design_grid$visit,
          times = nrow(scaled_centers)
        ),
        levels = visits,
        ordered = TRUE
      ),
      value = as.vector(
        t(scaled_centers)
      ),
      value_type = "scaled",
      estimator = estimator,
      stringsAsFactors = FALSE
    )
    
    raw_profiles <- data.frame(
      cluster_id = rep(
        rownames(raw_centers),
        each = ncol(raw_centers)
      ),
      group = rep(
        design_grid$group,
        times = nrow(raw_centers)
      ),
      visit = factor(
        rep(
          design_grid$visit,
          times = nrow(raw_centers)
        ),
        levels = visits,
        ordered = TRUE
      ),
      value = as.vector(
        t(raw_centers)
      ),
      value_type = "raw",
      estimator = estimator,
      stringsAsFactors = FALSE
    )
    
    object@results$cluster_profiles <-
      S4Vectors::DataFrame(
        rbind(
          scaled_profiles,
          raw_profiles
        )
      )
    
    # -------------------------------------------------------------------------
    # Store clustering diagnostics
    # -------------------------------------------------------------------------
    
    total_withinss <- NA_real_
    average_silhouette <- NA_real_
    
    if (estimator == "kmeans") {
      total_withinss <- fit$tot.withinss
    }
    
    if (
      estimator == "pam" &&
      !is.null(fit$silinfo$avg.width)
    ) {
      average_silhouette <-
        as.numeric(
          fit$silinfo$avg.width
        )
    }
    
    if (
      estimator == "kmeans" &&
      requireNamespace(
        "cluster",
        quietly = TRUE
      ) &&
      nrow(traj_for_clustering) > k
    ) {
      
      silhouette_result <- cluster::silhouette(
        cluster_assignments,
        stats::dist(
          traj_for_clustering
        )
      )
      
      average_silhouette <- mean(
        silhouette_result[
          ,
          "sil_width"
        ],
        na.rm = TRUE
      )
    }
    
    clustering_diagnostics <- data.frame(
      estimator = estimator,
      k = k,
      n_features = nrow(traj_for_clustering),
      n_profiles = ncol(traj_for_clustering),
      total_withinss = total_withinss,
      average_silhouette = average_silhouette,
      scale_features = scale_features,
      seed = seed,
      stringsAsFactors = FALSE
    )
    
    object@results$cluster_diagnostics <-
      S4Vectors::DataFrame(
        clustering_diagnostics
      )
    
    # -------------------------------------------------------------------------
    # Store the retained trajectory matrix
    # -------------------------------------------------------------------------
    
    trajectory_matrix_long <- data.frame(
      feature_id = rep(
        rownames(traj_mat),
        each = ncol(traj_mat)
      ),
      group = rep(
        design_grid$group,
        times = nrow(traj_mat)
      ),
      visit = factor(
        rep(
          design_grid$visit,
          times = nrow(traj_mat)
        ),
        levels = visits,
        ordered = TRUE
      ),
      raw_value = as.vector(
        t(traj_mat)
      ),
      clustering_value = as.vector(
        t(traj_for_clustering)
      ),
      cluster_id = rep(
        paste0(
          "cluster_",
          cluster_assignments[
            rownames(traj_mat)
          ]
        ),
        each = ncol(traj_mat)
      ),
      stringsAsFactors = FALSE
    )
    
    object@results$feature_trajectory_matrix <-
      S4Vectors::DataFrame(
        trajectory_matrix_long
      )
    
    # -------------------------------------------------------------------------
    # Record analysis history
    # -------------------------------------------------------------------------
    
    object@history <- c(
      object@history,
      list(
        list(
          step = "lot_from_clusters",
          assay = assay,
          assay_data = assay_data,
          estimator = estimator,
          n_clusters = k,
          n_features = nrow(traj_mat),
          n_profiles = ncol(traj_mat),
          scale_features = scale_features,
          standardise_ids = standardise_ids,
          seed = seed,
          visits = visits,
          groups = groups,
          omitted_group_visit_profiles =
            nrow(missing_profiles),
          time = Sys.time()
        )
      )
    )
    
    if (isTRUE(verbose)) {
      message(
        "Created ",
        k,
        " trajectory clusters using `",
        estimator,
        "` from ",
        nrow(traj_mat),
        " features and ",
        ncol(traj_mat),
        " group-by-visit profiles."
      )
    }
    
    methods::validObject(
      object
    )
    
    object
  }
)