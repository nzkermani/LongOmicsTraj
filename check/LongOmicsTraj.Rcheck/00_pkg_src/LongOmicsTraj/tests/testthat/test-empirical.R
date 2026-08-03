
test_that("empirical VMM builds correctly", {

  data("lot_glucold", package = "LongOmicsTraj")

  visits <- c("Baseline", "6 Months", "30 Months")

  res <- lot_from_empirical(
    lot_glucold,
    assay = "transcriptomics",
    visits = visits,
    overwrite = TRUE
  )

  expect_true(nrow(res@vmm) > 0)
})
