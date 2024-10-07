library(testthat)
library(Luciernaga)

test_that("QC_FilePrep returns a dataframe with more than 1 row", {
  # Prepare the experiment
  File_Location <- system.file("extdata", package = "Luciernaga")
  CSV_Pattern <- ".CSV$"
  CSV_Files <- list.files(path = File_Location, pattern = CSV_Pattern, full.names = TRUE)

  # Execute the experiment
  TidyData <- QC_FilePrep(CSV_Files, TrackChange = FALSE)

  # Did it return a data.frame
  expect_s3_class(TidyData, "data.frame")

  # Is it more than 1 row
  expect_gt(nrow(TidyData), 1)
})
