library(testthat)
library(LongOmicsTraj)

test_that("clustering works (kmeans + pam)", {

  data("lot_glucold", package = "LongOmicsTraj")
  visits <- c("Baseline", "6 Months", "30 Months")

  for (method in c("kmeans", "pam")) {

    res <- lot_from_clusters(
      lot_glucold,
      assay = "transcriptomics",
      visits = visits,
      k = 3,
      estimator = method,
      overwrite = TRUE,
      verbose = FALSE
    )

    expect_s4_class(res, "LongOmicsTraj")

    vmm <- as.data.frame(res@vmm)
    rel <- as.data.frame(res@relationships)

    expect_equal(length(unique(vmm$object_id)), 3)
    expect_true(nrow(rel) > 0)
  }
})
