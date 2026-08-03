test_that("lot_create returns a valid LongOmicsTraj object", {
  
  data(
    "lot_glucold_expression",
    package = "LongOmicsTraj"
  )
  
  data(
    "lot_glucold_clinical",
    package = "LongOmicsTraj"
  )
  
  object <- lot_create(
    expression = lot_glucold_expression,
    metadata = lot_glucold_clinical,
    assay = "transcriptomics",
    visit_levels = c(
      "Baseline",
      "6 Months",
      "30 Months"
    ),
    time_values = c(
      "Baseline" = 0,
      "6 Months" = 6,
      "30 Months" = 30
    ),
    feature_type = "gene"
  )
  
  expect_s4_class(
    object,
    "LongOmicsTraj"
  )
  
  expect_true(
    methods::validObject(object)
  )
})


test_that("lot_create aligns shuffled metadata", {
  
  data(
    "lot_glucold_expression",
    package = "LongOmicsTraj"
  )
  
  data(
    "lot_glucold_clinical",
    package = "LongOmicsTraj"
  )
  
  shuffled_metadata <- lot_glucold_clinical[
    rev(seq_len(nrow(lot_glucold_clinical))),
    ,
    drop = FALSE
  ]
  
  object <- lot_create(
    expression = lot_glucold_expression,
    metadata = shuffled_metadata
  )
  
  se <- MultiAssayExperiment::experiments(
    object@mae
  )[["transcriptomics"]]
  
  stored_metadata <- as.data.frame(
    SummarizedExperiment::colData(se)
  )
  
  expect_identical(
    colnames(lot_glucold_expression),
    stored_metadata$sample_id
  )
})


test_that("lot_create rejects missing metadata columns", {
  
  data(
    "lot_glucold_expression",
    package = "LongOmicsTraj"
  )
  
  data(
    "lot_glucold_clinical",
    package = "LongOmicsTraj"
  )
  
  incomplete_metadata <- lot_glucold_clinical[
    ,
    setdiff(
      colnames(lot_glucold_clinical),
      "subject_id"
    ),
    drop = FALSE
  ]
  
  expect_error(
    lot_create(
      expression = lot_glucold_expression,
      metadata = incomplete_metadata
    ),
    "subject_id"
  )
})


test_that("lot_create rejects unmatched sample IDs", {
  
  data(
    "lot_glucold_expression",
    package = "LongOmicsTraj"
  )
  
  data(
    "lot_glucold_clinical",
    package = "LongOmicsTraj"
  )
  
  altered_metadata <- lot_glucold_clinical
  altered_metadata$sample_id[1] <- "unknown_sample"
  
  expect_error(
    lot_create(
      expression = lot_glucold_expression,
      metadata = altered_metadata
    ),
    "sample IDs do not match"
  )
})