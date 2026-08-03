test_that("lot_cluster_purity calculates expected purity", {
  
  assignments <- data.frame(
    gene_id = paste0("G", 1:10),
    cluster = c(
      rep(1, 5),
      rep(2, 5)
    ),
    group = "Withdrawal",
    topology_label = c(
      rep("up_down", 4),
      "flat_flat",
      rep("down_up", 5)
    )
  )
  
  result <- lot_cluster_purity(
    assignments,
    purity_threshold = 0.90
  )
  
  expect_equal(
    nrow(result$cluster_group_purity),
    2
  )
  
  expect_equal(
    sort(result$cluster_group_purity$purity),
    c(0.8, 1.0)
  )
  
  expect_equal(
    result$overall$passing_fraction,
    0.5
  )
})