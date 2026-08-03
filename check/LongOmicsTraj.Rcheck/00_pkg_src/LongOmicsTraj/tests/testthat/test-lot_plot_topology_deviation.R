test_that("topology deviation plot returns a ggplot", {
  
  object <- make_flexmix_plot_object()
  
  p <- lot_plot_topology_deviation(
    object = object,
    group = "Example treatment"
  )
  
  expect_s3_class(
    p,
    "ggplot"
  )
})


test_that("topology deviation plot supports top-n filtering", {
  
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
  
  expect_lte(
    max(
      table(
        p$data$cluster
      )
    ),
    1L
  )
})


test_that("topology deviation plot supports deviation filtering", {
  
  object <- make_flexmix_plot_object()
  
  p <- lot_plot_topology_deviation(
    object = object,
    group = "Example treatment",
    min_abs_deviation = 0.05
  )
  
  expect_s3_class(
    p,
    "ggplot"
  )
  
  expect_true(
    all(
      abs(
        p$data$deviation
      ) >= 0.05
    )
  )
})


test_that("topology deviation plot rejects invalid deviation thresholds", {
  
  object <- make_flexmix_plot_object()
  
  expect_error(
    lot_plot_topology_deviation(
      object = object,
      min_abs_deviation = 1.5
    ),
    "must be between 0 and 1"
  )
})


test_that("topology deviation plot rejects missing FlexMix results", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  object <- lot_glucold_ots
  
  object@results$flexmix_ots_distribution <- NULL
  
  expect_error(
    lot_plot_topology_deviation(
      object = object
    ),
    "No FlexMix OTS distribution was found"
  )
})