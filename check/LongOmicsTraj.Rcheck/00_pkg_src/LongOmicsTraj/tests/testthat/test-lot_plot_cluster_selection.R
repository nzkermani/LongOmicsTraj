test_that("cluster selection plot returns a ggplot", {
  
  object <- make_flexmix_plot_object()
  
  p <- lot_plot_cluster_selection(
    object = object
  )
  
  expect_s3_class(
    p,
    "ggplot"
  )
})


test_that("cluster selection plot accepts a subset of metrics", {
  
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
  
  expect_equal(
    length(
      unique(
        p$data$metric
      )
    ),
    2L
  )
})


test_that("cluster selection plot rejects unsupported metrics", {
  
  object <- make_flexmix_plot_object()
  
  expect_error(
    lot_plot_cluster_selection(
      object = object,
      metrics = "banana"
    ),
    "Unsupported metrics"
  )
})


test_that("cluster selection plot rejects missing FlexMix results", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  object <- lot_glucold_ots
  
  object@results$flexmix_model_comparison <- NULL
  
  expect_error(
    lot_plot_cluster_selection(
      object = object
    ),
    "No FlexMix model-comparison results were found"
  )
})