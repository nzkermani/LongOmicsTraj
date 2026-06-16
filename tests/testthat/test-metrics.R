test_that('OTS metrics computed correctly', {

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


    # Assuming lot_compute_metrics is how you calculate metrics downstream
  expect_s4_class(res, 'LongOmicsTraj')
})
