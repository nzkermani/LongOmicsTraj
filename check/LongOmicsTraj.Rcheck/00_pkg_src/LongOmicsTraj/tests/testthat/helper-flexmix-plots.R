make_flexmix_plot_object <- function() {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  object <- lot_glucold_ots
  
  group_name <- "Example treatment"
  
  assignments <- data.frame(
    gene_id = c(
      "G1", "G2", "G3",
      "G4", "G5", "G6"
    ),
    cluster = factor(
      c(
        1, 1, 1,
        2, 2, 2
      )
    ),
    topology_label = c(
      "flat_flat",
      "flat_flat",
      "flat_up",
      "down_up",
      "down_up",
      "flat_flat"
    ),
    group = group_name,
    k = 2L,
    stringsAsFactors = FALSE
  )
  
  cluster_trajectories <- data.frame(
    cluster = factor(
      rep(
        c(1, 2),
        each = 3
      )
    ),
    visit = factor(
      rep(
        c(
          "Baseline",
          "Visit 2",
          "Visit 3"
        ),
        times = 2
      ),
      levels = c(
        "Baseline",
        "Visit 2",
        "Visit 3"
      ),
      ordered = TRUE
    ),
    time = rep(
      0:2,
      times = 2
    ),
    mean_clustering_value = c(
      -0.2, 0.0, 0.2,
      0.8, 0.0, -0.8
    ),
    mean_trajectory_value = c(
      5.0, 5.1, 5.2,
      7.0, 6.0, 5.0
    ),
    n_genes = rep(
      3L,
      6
    ),
    group = group_name,
    k = 2L,
    stringsAsFactors = FALSE
  )
  
  ots_distribution <- data.frame(
    cluster = factor(
      c(
        1, 1,
        2, 2
      )
    ),
    group = group_name,
    topology_label = c(
      "flat_flat",
      "flat_up",
      "down_up",
      "flat_flat"
    ),
    n_genes = c(
      2L, 1L,
      2L, 1L
    ),
    proportion = c(
      2 / 3,
      1 / 3,
      2 / 3,
      1 / 3
    ),
    stringsAsFactors = FALSE
  )
  
  cluster_purity <- data.frame(
    cluster = factor(
      c(1, 2)
    ),
    group = group_name,
    n_genes = c(
      3L,
      3L
    ),
    dominant_topology = c(
      "flat_flat",
      "down_up"
    ),
    dominant_count = c(
      2L,
      2L
    ),
    purity = c(
      2 / 3,
      2 / 3
    ),
    entropy = c(
      0.6365142,
      0.6365142
    ),
    n_topologies = c(
      2L,
      2L
    ),
    normalised_entropy = c(
      0.9182958,
      0.9182958
    ),
    passes_purity = c(
      FALSE,
      FALSE
    ),
    stringsAsFactors = FALSE
  )
  
  trajectory_data <- expand.grid(
    gene_id = paste0(
      "G",
      1:6
    ),
    visit = factor(
      c(
        "Baseline",
        "Visit 2",
        "Visit 3"
      ),
      levels = c(
        "Baseline",
        "Visit 2",
        "Visit 3"
      ),
      ordered = TRUE
    ),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  trajectory_data <- trajectory_data[
    order(
      trajectory_data$gene_id,
      trajectory_data$visit
    ),
    ,
    drop = FALSE
  ]
  
  trajectory_data$time <- rep(
    0:2,
    times = 6
  )
  
  trajectory_data$mean_expression <- c(
    5.0, 5.0, 5.1,
    4.9, 5.0, 5.1,
    5.0, 5.1, 5.4,
    7.0, 6.0, 5.0,
    7.2, 6.1, 5.1,
    5.1, 5.0, 4.9
  )
  
  trajectory_data$clustering_value <- ave(
    trajectory_data$mean_expression,
    trajectory_data$gene_id,
    FUN = function(x) {
      as.numeric(
        scale(x)
      )
    }
  )
  
  trajectory_data$time_scaled <- rep(
    as.numeric(
      scale(
        0:2
      )
    ),
    times = 6
  )
  
  trajectory_data$n_samples <- NA_integer_
  
  trajectory_data <- merge(
    trajectory_data,
    assignments[
      ,
      c(
        "gene_id",
        "topology_label"
      ),
      drop = FALSE
    ],
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  model_comparison <- data.frame(
    k = 2:5,
    trajectory_source = "estimated",
    trajectory_estimator = "masigpro",
    BIC = c(
      470,
      425,
      390,
      405
    ),
    AIC = c(
      440,
      380,
      330,
      325
    ),
    logLik = c(
      -211,
      -176,
      -145,
      -138
    ),
    n_parameters = c(
      9,
      14,
      19,
      24
    ),
    minimum_purity = c(
      0.50,
      0.60,
      0.72,
      0.70
    ),
    mean_purity = c(
      0.58,
      0.66,
      0.82,
      0.80
    ),
    weighted_purity = c(
      0.60,
      0.68,
      0.85,
      0.83
    ),
    passing_fraction = c(
      0.00,
      0.33,
      0.75,
      0.60
    ),
    minimum_cluster_size = c(
      18,
      12,
      7,
      4
    ),
    passes_topology_rule = c(
      FALSE,
      FALSE,
      FALSE,
      FALSE
    ),
    passes_size_rule = c(
      TRUE,
      TRUE,
      TRUE,
      FALSE
    ),
    acceptable = FALSE,
    stringsAsFactors = FALSE
  )
  
  object@results$flexmix_gene_assignments <-
    S4Vectors::DataFrame(
      assignments
    )
  
  object@results$flexmix_cluster_trajectories <-
    S4Vectors::DataFrame(
      cluster_trajectories
    )
  
  object@results$flexmix_ots_distribution <-
    S4Vectors::DataFrame(
      ots_distribution
    )
  
  object@results$flexmix_cluster_purity <-
    S4Vectors::DataFrame(
      cluster_purity
    )
  
  object@results$flexmix_trajectory_data <-
    S4Vectors::DataFrame(
      trajectory_data
    )
  
  object@results$flexmix_model_comparison <-
    S4Vectors::DataFrame(
      model_comparison
    )
  
  object@results$flexmix_selected_k <- 4L
  
  object
}