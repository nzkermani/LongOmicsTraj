test_that("default centring leaves deltas unchanged", {

  data("lot_glucold", package = "LongOmicsTraj")
  visits <- c("Baseline", "6 Months", "30 Months")

  res <- lot_from_empirical(lot_glucold, assay = "transcriptomics",
                            visits = visits, overwrite = TRUE)

  d0 <- lot_compute_deltas(res, assay = "transcriptomics",
                           visits = visits, overwrite = TRUE)
  d1 <- lot_compute_deltas(res, assay = "transcriptomics",
                           visits = visits, overwrite = TRUE, center = "none")

  expect_equal(as.data.frame(d0@deltas), as.data.frame(d1@deltas))
  expect_null(d1@results$visit_centring)
})

test_that("median centring gives zero median change per group and transition", {

  data("lot_glucold", package = "LongOmicsTraj")
  visits <- c("Baseline", "6 Months", "30 Months")

  res <- lot_from_empirical(lot_glucold, assay = "transcriptomics",
                            visits = visits, overwrite = TRUE)

  # add a visit-wide shift of +2 at the last visit in every group
  vmm <- as.data.frame(res@vmm)
  vmm$estimated_value[vmm$visit == "30 Months"] <-
    vmm$estimated_value[vmm$visit == "30 Months"] + 2
  res@vmm <- S4Vectors::DataFrame(vmm)

  raw <- as.data.frame(lot_compute_deltas(res, assay = "transcriptomics",
                                          visits = visits, overwrite = TRUE)@deltas)
  cen <- lot_compute_deltas(res, assay = "transcriptomics", visits = visits,
                            overwrite = TRUE, center = "median")
  dc <- as.data.frame(cen@deltas)

  # the shift dominates uncentred late changes
  late_raw <- raw$delta_hat[raw$to_visit == "30 Months"]
  expect_gt(stats::median(late_raw), 1)

  # after centring, the shift is gone
  late_cen <- dc$delta_hat[dc$to_visit == "30 Months"]
  expect_lt(abs(stats::median(late_cen)), 0.5)

  # offsets are recorded for every group x visit
  off <- as.data.frame(cen@results$visit_centring)
  expect_equal(nrow(off), length(unique(vmm$group)) * length(visits))
  expect_true(all(off$method == "median"))
})

test_that("centring argument is validated", {

  data("lot_glucold", package = "LongOmicsTraj")
  visits <- c("Baseline", "6 Months", "30 Months")
  res <- lot_from_empirical(lot_glucold, assay = "transcriptomics",
                            visits = visits, overwrite = TRUE)

  expect_error(
    lot_compute_deltas(res, assay = "transcriptomics", visits = visits,
                       overwrite = TRUE, center = "zscore")
  )
})
