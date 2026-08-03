test_that("cluster trajectory-arm plot returns a ggplot", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  visits <- c(
    "Baseline",
    "6 Months",
    "30 Months"
  )
  
  clustered <- lot_from_clusters(
    object = lot_glucold_ots,
    assay = "transcriptomics",
    visits = visits,
    k = 3,
    estimator = "kmeans",
    require_gene_ots = TRUE,
    overwrite = TRUE,
    verbose = FALSE
  )
  
  p <- lot_plot_cluster_trajectory_arms(
    clustered,
    top_n = 4
  )
  
  expect_s3_class(
    p,
    "ggplot"
  )
})


test_that("cluster trajectory-arm plot filters groups", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  visits <- c(
    "Baseline",
    "6 Months",
    "30 Months"
  )
  
  clustered <- lot_from_clusters(
    object = lot_glucold_ots,
    assay = "transcriptomics",
    visits = visits,
    k = 3,
    estimator = "kmeans",
    require_gene_ots = TRUE,
    overwrite = TRUE,
    verbose = FALSE
  )
  
  selected_group <- "Placebo"
  
  p <- lot_plot_cluster_trajectory_arms(
    clustered,
    groups = selected_group
  )
  
  expect_true(
    all(
      grepl(
        selected_group,
        p$layers[[1]]$data$panel_id,
        fixed = TRUE
      )
    )
  )
})