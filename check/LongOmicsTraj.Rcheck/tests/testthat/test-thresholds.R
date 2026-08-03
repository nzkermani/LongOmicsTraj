test_that('thresholds are computed and valid', {

  data("lot_glucold", package = "LongOmicsTraj")

  visits <- c("Baseline", "6 Months", "30 Months")

  res <- lot_from_empirical(
    lot_glucold,
    assay = "transcriptomics",
    visits = visits,
    overwrite = TRUE
  )


  res <- lot_compute_deltas(
    res,
    assay = 'transcriptomics',
    visits = c("Baseline", "6 Months", "30 Months"),
    overwrite = TRUE
  )

  res <- lot_compute_thresholds(
    res,
    assay = 'transcriptomics',
    visits = c("Baseline", "6 Months", "30 Months"),
    overwrite = TRUE
  )

    expect_s4_class(res, 'LongOmicsTraj')
})
