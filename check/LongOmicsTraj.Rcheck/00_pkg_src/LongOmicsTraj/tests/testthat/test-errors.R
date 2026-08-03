test_that("functions fail correctly", {
  
  data(
    "lot_glucold",
    package = "LongOmicsTraj"
  )
  
  experiments <- MultiAssayExperiment::experiments(
    lot_glucold@mae
  )
  
  first_experiment <- experiments[[1L]]
  
  metadata <- as.data.frame(
    SummarizedExperiment::colData(
      first_experiment
    )
  )
  
  available_visits <- unique(
    as.character(
      metadata$visit[
        !is.na(metadata$visit)
      ]
    )
  )
  
  expect_gte(
    length(available_visits),
    2L
  )
  
  expect_error(
    lot_from_clusters(
      object = lot_glucold,
      assay = "fake",
      visits = available_visits
    ),
    "was not found"
  )
  
  expect_error(
    lot_compute_deltas(
      lot_glucold
    ),
    "VMM is empty"
  )
})