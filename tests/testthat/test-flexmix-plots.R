# -------------------------------------------------------------------------
# lot_plot_cluster_topology()
# -------------------------------------------------------------------------

test_that("lot_plot_cluster_topology returns a patchwork plot", {
  
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


test_that("lot_plot_cluster_topology works without feature lines", {
  
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


test_that("lot_plot_cluster_topology rejects missing results", {
  
  testthat::skip_if_not_installed(
    "patchwork"
  )
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  object <- lot_glucold_ots
  
  flexmix_names <- grep(
    "^flexmix_",
    names(object@results),
    value = TRUE
  )
  
  for (result_name in flexmix_names) {
    object@results[[result_name]] <- NULL
  }
  
  expect_error(
    lot_plot_cluster_topology(
      object
    ),
    "FlexMix results are missing"
  )
})


# -------------------------------------------------------------------------
# lot_plot_cluster_selection()
# -------------------------------------------------------------------------

test_that("lot_plot_cluster_selection returns a ggplot", {
  
  object <- make_flexmix_plot_object()
  
  p <- lot_plot_cluster_selection(
    object = object,
    purity_threshold = 0.90,
    required_passing_fraction = 0.80,
    min_cluster_size = 5
  )
  
  expect_s3_class(
    p,
    "ggplot"
  )
})


test_that("lot_plot_cluster_selection accepts selected metrics", {
  
  object <- make_flexmix_plot_object()
  
  p <- lot_plot_cluster_selection(
    object = object,
    metrics = c(
      "BIC",
      "weighted_purity"
    )
  )
  
  expect_s3_class(
    p,
    "ggplot"
  )
})


test_that("lot_plot_cluster_selection rejects unsupported metrics", {
  
  object <- make_flexmix_plot_object()
  
  expect_error(
    lot_plot_cluster_selection(
      object = object,
      metrics = "unsupported_metric"
    ),
    "Unsupported metrics"
  )
})


# -------------------------------------------------------------------------
# lot_plot_topology_deviation()
# -------------------------------------------------------------------------

test_that("lot_plot_topology_deviation returns a ggplot", {
  
  object <- make_flexmix_plot_object()
  
  p <- lot_plot_topology_deviation(
    object = object,
    group = "Example treatment",
    min_abs_deviation = 0
  )
  
  expect_s3_class(
    p,
    "ggplot"
  )
})


test_that("lot_plot_topology_deviation supports top-n filtering", {
  
  object <- make_flexmix_plot_object()
  
  p <- lot_plot_topology_deviation(
    object = object,
    group = "Example treatment",
    top_n = 1L
  )
  
  expect_s3_class(
    p,
    "ggplot"
  )
})


test_that("lot_plot_topology_deviation rejects invalid deviation", {
  
  object <- make_flexmix_plot_object()
  
  expect_error(
    lot_plot_topology_deviation(
      object = object,
      min_abs_deviation = 1.5
    ),
    "must be between 0 and 1"
  )
})