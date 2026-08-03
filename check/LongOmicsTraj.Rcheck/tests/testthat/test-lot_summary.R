test_that("lot_summary returns a structured summary object", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  s <- lot_summary(
    lot_glucold_ots
  )
  
  expect_s3_class(
    s,
    "lot_summary"
  )
  
  expect_equal(
    s$n_features,
    71
  )
  
  expect_equal(
    s$n_samples,
    221
  )
  
  expect_true(
    nrow(s$ots_distribution) > 0L
  )
})


test_that("lot_summary reports estimators and OTS classes", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  s <- lot_summary(
    lot_glucold_ots
  )
  
  expect_true(
    "masigpro" %in% s$estimators
  )
  
  expect_gt(
    s$n_ots_classes,
    0L
  )
})


test_that("lot_summary filters by estimator", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  s <- lot_summary(
    lot_glucold_ots,
    estimator = "masigpro"
  )
  
  expect_equal(
    s$selected_estimator,
    "masigpro"
  )
  
  expect_true(
    nrow(s$ots_distribution) > 0L
  )
})


test_that("lot_summary rejects an unknown estimator", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  expect_error(
    lot_summary(
      lot_glucold_ots,
      estimator = "unknown_estimator"
    ),
    "Available estimators"
  )
})


test_that("lot_summary respects top_n", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  s <- lot_summary(
    lot_glucold_ots,
    top_n = 2
  )
  
  counts_per_group <- table(
    s$top_ots$group
  )
  
  expect_true(
    all(counts_per_group <= 2L)
  )
})


test_that("printed lot_summary contains key sections", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  s <- lot_summary(
    lot_glucold_ots
  )
  
  output <- capture.output(
    print(s)
  )
  
  output_text <- paste(
    output,
    collapse = "\n"
  )
  
  expect_match(
    output_text,
    "LongOmicsTraj summary"
  )
  
  expect_match(
    output_text,
    "Omics Trajectory Signatures"
  )
  
  expect_match(
    output_text,
    "masigpro"
  )
  
  expect_match(
    output_text,
    "Dataset"
  )
  
  expect_match(
    output_text,
    "Analysis status"
  )
})