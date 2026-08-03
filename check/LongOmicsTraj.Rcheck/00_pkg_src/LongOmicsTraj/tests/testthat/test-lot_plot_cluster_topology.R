test_that("cluster topology plot returns a patchwork object", {
  
  testthat::skip_if_not_installed(
    "patchwork"
  )
  
  object <- make_flexmix_plot_object()
  
  p <- lot_plot_cluster_topology(
    object = object,
    group = "Example treatment",
    value_type = "scaled",
    show_features = TRUE,
    max_feature_lines = 2L,
    purity_threshold = 0.90
  )
  
  expect_s3_class(
    p,
    "patchwork"
  )
})


test_that("cluster topology plot works without feature trajectories", {
  
  testthat::skip_if_not_installed(
    "patchwork"
  )
  
  object <- make_flexmix_plot_object()
  
  p <- lot_plot_cluster_topology(
    object = object,
    group = "Example treatment",
    value_type = "raw",
    show_features = FALSE
  )
  
  expect_s3_class(
    p,
    "patchwork"
  )
})


test_that("cluster topology plot supports cluster filtering", {
  
  testthat::skip_if_not_installed(
    "patchwork"
  )
  
  object <- make_flexmix_plot_object()
  
  p <- lot_plot_cluster_topology(
    object = object,
    clusters = "1",
    group = "Example treatment",
    show_features = FALSE
  )
  
  expect_s3_class(
    p,
    "patchwork"
  )
})


test_that("cluster topology plot supports purity sorting", {
  
  testthat::skip_if_not_installed(
    "patchwork"
  )
  
  object <- make_flexmix_plot_object()
  
  p <- lot_plot_cluster_topology(
    object = object,
    group = "Example treatment",
    sort_by = "purity",
    show_features = FALSE
  )
  
  expect_s3_class(
    p,
    "patchwork"
  )
})


test_that("cluster topology plot rejects an invalid purity threshold", {
  
  testthat::skip_if_not_installed(
    "patchwork"
  )
  
  object <- make_flexmix_plot_object()
  
  expect_error(
    lot_plot_cluster_topology(
      object = object,
      group = "Example treatment",
      purity_threshold = 1.5
    ),
    "must be between 0 and 1"
  )
})


test_that("cluster topology plot rejects missing FlexMix results", {
  
  testthat::skip_if_not_installed(
    "patchwork"
  )
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  object <- lot_glucold_ots
  
  flexmix_result_names <- grep(
    "^flexmix_",
    names(object@results),
    value = TRUE
  )
  
  for (result_name in flexmix_result_names) {
    object@results[[result_name]] <- NULL
  }
  
  expect_error(
    lot_plot_cluster_topology(
      object = object
    ),
    "FlexMix results are missing"
  )
})