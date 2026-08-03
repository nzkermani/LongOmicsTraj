test_that("lot_plot_ots_distribution returns a ggplot object", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  p <- lot_plot_ots_distribution(
    lot_glucold_ots,
    estimator = "masigpro"
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("proportions sum to one within each group", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  p <- lot_plot_ots_distribution(
    lot_glucold_ots,
    estimator = "masigpro",
    type = "proportion"
  )
  
  group_sums <- aggregate(
    value ~ group,
    data = p$data,
    FUN = sum
  )
  
  expect_equal(
    group_sums$value,
    rep(1, nrow(group_sums)),
    tolerance = 1e-8
  )
})


test_that("count output uses raw counts", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  p <- lot_plot_ots_distribution(
    lot_glucold_ots,
    estimator = "masigpro",
    type = "count"
  )
  
  expect_equal(
    p$data$value,
    p$data$count
  )
})


test_that("unknown estimator gives an informative error", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  expect_error(
    lot_plot_ots_distribution(
      lot_glucold_ots,
      estimator = "unknown_estimator"
    ),
    "Available estimators"
  )
})