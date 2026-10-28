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


test_that("relative centring removes group-specific but keeps shared shifts", {

  data("lot_glucold", package = "LongOmicsTraj")
  visits <- c("Baseline", "6 Months", "30 Months")
  res <- lot_from_empirical(lot_glucold, assay = "transcriptomics",
                            visits = visits, overwrite = TRUE)
  vmm <- as.data.frame(res@vmm)
  groups <- unique(vmm$group)
  # shared shift of +1 at the last visit in every group, plus +2 more in one group
  late <- vmm$visit == "30 Months"
  vmm$estimated_value[late] <- vmm$estimated_value[late] + 1
  one <- late & vmm$group == groups[1]
  vmm$estimated_value[one] <- vmm$estimated_value[one] + 2
  res@vmm <- S4Vectors::DataFrame(vmm)

  raw <- as.data.frame(lot_compute_deltas(res, assay = "transcriptomics",
                                          visits = visits, overwrite = TRUE)@deltas)
  rel <- suppressWarnings(lot_compute_deltas(res, assay = "transcriptomics",
                                             visits = visits, overwrite = TRUE,
                                             center = "relative"))
  dr <- as.data.frame(rel@deltas)
  med <- function(d, g) stats::median(d$delta_hat[d$to_visit == "30 Months" & d$group == g])

  # the group-specific extra shift is removed: late medians now agree across groups
  meds <- vapply(groups, function(g) med(dr, g), numeric(1))
  expect_lt(diff(range(meds)), 0.1)
  # the shared shift is kept when groups agree (with two groups it is their average)
  expect_gt(min(meds), 0.5)
  expect_true(all(as.data.frame(rel@results$visit_centring)$method == "relative"))
  # transitions not affected by any shift keep their raw values up to a small offset
  early <- raw$to_visit == "6 Months"
  expect_lt(max(abs(dr$delta_hat[early] - raw$delta_hat[early])), 0.5)
})
