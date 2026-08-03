test_that("full pipeline is stable", {

  data("lot_glucold", package = "LongOmicsTraj")

  visits <- c("Baseline", "6 Months", "30 Months")

  # ------------------------------------------------------------
  # 1. Empirical VMM
  # ------------------------------------------------------------
  res <- lot_from_empirical(
    object = lot_glucold,
    assay = "transcriptomics",
    visits = visits,
    overwrite = TRUE
  )

  expect_s4_class(res, "LongOmicsTraj")
  expect_true(nrow(res@vmm) > 0)

  # ------------------------------------------------------------
  # 2. Deltas
  # ------------------------------------------------------------
  res <- lot_compute_deltas(
    object = res,
    assay = "transcriptomics",
    visits = visits,
    overwrite = TRUE
  )

  expect_true(nrow(res@deltas) > 0)

  # ------------------------------------------------------------
  # 3. Thresholds
  # ------------------------------------------------------------
  res <- lot_compute_thresholds(
    object = res,
    assay = "transcriptomics",
    visits = visits,
    overwrite = TRUE
  )

  expect_true(nrow(res@thresholds) > 0)

  # ------------------------------------------------------------
  # 4. OTS
  # ------------------------------------------------------------
  res <- lot_map_to_ots(
    object = res,
    visits = visits,
    overwrite = TRUE
  )

  expect_true(nrow(res@ots) > 0)

  ots_df <- as.data.frame(res@ots)

  # ------------------------------------------------------------
  # 5. VALIDATION (THIS IS WHAT MATTERS)
  # ------------------------------------------------------------

  # Structural checks
  expect_true("topology_label" %in% colnames(ots_df))
  expect_false(any(is.na(ots_df$topology_label)))

  # Expected topology space
  valid_topologies <- c(
    "up_up","up_down","up_flat",
    "down_up","down_down","down_flat",
    "flat_up","flat_down","flat_flat"
  )

  expect_true(all(ots_df$topology_label %in% valid_topologies))

  # Transition consistency
  expect_true(all(ots_df$n_transitions == length(visits) - 1))

  # Size sanity (not fragile, but meaningful)
  expect_gt(nrow(ots_df), 50)

})
