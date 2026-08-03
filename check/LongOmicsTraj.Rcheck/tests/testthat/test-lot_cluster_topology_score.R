test_that("lot_cluster_topology_score detects duplicated topologies", {
  
  assignments <- data.frame(
    gene_id = paste0("G", 1:12),
    cluster = c(
      rep(1, 4),
      rep(2, 4),
      rep(3, 4)
    ),
    group = "Withdrawal",
    topology_label = c(
      rep("flat_flat", 4),
      rep("flat_flat", 3),
      "flat_up",
      rep("down_up", 4)
    ),
    stringsAsFactors = FALSE
  )
  
  purity <- lot_cluster_purity(
    assignments,
    purity_threshold = 0.90
  )
  
  score <- lot_cluster_topology_score(
    purity = purity,
    k = 3,
    bic = 125.4,
    minimum_cluster_size = 4
  )
  
  expect_s3_class(
    score,
    "lot_cluster_topology_score"
  )
  
  expect_equal(
    score$overall$requested_k,
    3L
  )
  
  expect_equal(
    score$overall$observed_clusters,
    3L
  )
  
  expect_equal(
    score$overall$unique_dominant_topologies,
    2L
  )
  
  expect_equal(
    score$overall$duplicate_topology_count,
    1L
  )
  
  expect_equal(
    score$overall$topology_coverage,
    2 / 3
  )
  
  expect_equal(
    score$overall$BIC,
    125.4
  )
})


test_that("topology coverage is one when all dominant labels differ", {
  
  assignments <- data.frame(
    gene_id = paste0("G", 1:9),
    cluster = rep(
      1:3,
      each = 3
    ),
    group = "Withdrawal",
    topology_label = c(
      rep("flat_flat", 3),
      rep("up_down", 3),
      rep("down_up", 3)
    ),
    stringsAsFactors = FALSE
  )
  
  purity <- lot_cluster_purity(
    assignments
  )
  
  score <- lot_cluster_topology_score(
    purity,
    k = 3
  )
  
  expect_equal(
    score$overall$topology_coverage,
    1
  )
  
  expect_equal(
    score$overall$duplicate_topology_count,
    0L
  )
})


test_that("lot_cluster_topology_score validates required columns", {
  
  invalid_data <- data.frame(
    cluster = 1,
    group = "Withdrawal"
  )
  
  expect_error(
    lot_cluster_topology_score(
      invalid_data
    ),
    "Purity results are missing"
  )
})