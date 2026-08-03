test_that('OTS mapping works', {

  data("lot_glucold", package = "LongOmicsTraj")

  visits <- c("Baseline", "6 Months", "30 Months")

  res <- lot_from_empirical(
    lot_glucold,
    assay = "transcriptomics",
    visits = visits,
    overwrite = TRUE
  )

    expect_s4_class(res, 'LongOmicsTraj')
})
