test_that("QC_ReferenceLibrary returns a dataframe with more than 1 row", {
  # Prepare and Execute the experiment

  TheData <- QC_ReferenceLibrary(FluorNameContains = "FITC", NumberDetectors=64)

  # Is it more than 1 row
  expect_gt(nrow(TheData), 1)
})
