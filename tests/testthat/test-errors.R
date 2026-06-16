test_that("functions fail correctly", {

  data("lot_glucold", package = "LongOmicsTraj")

  expect_error(
    lot_from_clusters(lot_glucold, assay = "fake"),
    "Assay not found"
  )

  expect_error(
    lot_compute_deltas(lot_glucold),
    "VMM is empty"
  )
})
