test_that("lot_plot_trajectory returns a ggplot", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  p <- lot_plot_trajectory(
    lot_glucold_ots,
    features = "ATP2A3"
  )
  
  expect_s3_class(
    p,
    "ggplot"
  )
})


test_that("lot_plot_trajectory supports selected groups", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  selected_groups <- c(
    "ICS 30 months",
    "Placebo"
  )
  
  p <- lot_plot_trajectory(
    lot_glucold_ots,
    features = "ATP2A3",
    groups = selected_groups
  )
  
  plotted_groups <- unique(
    as.character(
      p$layers[[1]]$data$group
    )
  )
  
  expect_setequal(
    plotted_groups,
    selected_groups
  )
})


test_that("lot_plot_trajectory supports multiple features", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  p <- lot_plot_trajectory(
    lot_glucold_ots,
    features = c(
      "ATP2A3",
      "CAT"
    )
  )
  
  expect_setequal(
    unique(p$layers[[1]]$data$feature),
    c("ATP2A3", "CAT")
  )
})


test_that("lot_plot_trajectory rejects missing features", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  expect_error(
    suppressWarnings(
      lot_plot_trajectory(
        lot_glucold_ots,
        features = "not_a_gene"
      )
    ),
    "None of the requested features"
  )
})