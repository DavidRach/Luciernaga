test_that("QC_ProspectiveAdditions returns a dataframe with more than 1 row", {
  # Prepare the test

  Folder_Location <- system.file("extdata", package = "Luciernaga")
  ThePanelLocation <- list.files(Folder_Location, pattern="^Panel.csv",
                                 full.names=TRUE)
  OutPath <- file.path("C:", "Users", "JohnDoe", "Desktop")

  # Execute the test
  ProspectiveAdditions <- QC_ProspectiveAdditions(path=ThePanelLocation,
                                                  NumberDetectors=64,
                                                  TheCutoff=0.9,
                                                  returnAll=FALSE,
                                                  returnCSV=FALSE,
                                                  filename="ProspectiveAdditions",
                                                  outpath=OutPath)
  # Did it return a data.frame
  expect_s3_class(ProspectiveAdditions, "data.frame")

  # Is it more than 1 row
  expect_gt(nrow(ProspectiveAdditions), 1)
})
