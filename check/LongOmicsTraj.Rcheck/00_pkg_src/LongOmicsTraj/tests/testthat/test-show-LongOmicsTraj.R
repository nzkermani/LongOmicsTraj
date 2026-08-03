test_that("show displays a useful object summary", {
  
  data(
    "lot_glucold_ots",
    package = "LongOmicsTraj"
  )
  
  output <- capture.output(
    show(
      lot_glucold_ots
    )
  )
  
  output_text <- paste(
    output,
    collapse = "\n"
  )
  
  expect_match(
    output_text,
    "LongOmicsTraj object"
  )
  
  expect_match(
    output_text,
    "transcriptomics"
  )
  
  expect_match(
    output_text,
    "71 features"
  )
  
  expect_match(
    output_text,
    "221 samples"
  )
  
  expect_match(
    output_text,
    "Baseline"
  )
  
  expect_match(
    output_text,
    "masigpro"
  )
})