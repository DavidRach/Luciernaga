test_that("QC_SimilarFluorophores returns a dataframe with more than 1 row", {
  # Prepare and Execute the experiment
  Results <- QC_SimilarFluorophores(TheFluorophore="Spark Blue 550",
                                    NumberDetectors=64, NumberHits = 10,
                                    returnPlots=FALSE)

  # Is it more than 1 row
  expect_gt(nrow(Results), 1)
})
