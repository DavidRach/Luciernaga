test_that("QC_LibraryParse returns a ggplot object", {
  library(flowCore)
  library(flowWorkspace)
  library(purrr)

  # Prepare the experiment
  Folder_Location <- system.file("extdata", package = "Luciernaga")
  XML_Pattern <- ".XML$"
  XML_Files <- list.files(path = Folder_Location, pattern = XML_Pattern,
                          full.names = TRUE, recursive = FALSE)

  # Execute the experiment
  SinglePlot <- QC_LibraryParse(XML_Files[2], returntype="plots", references=FALSE)

  # Did it return a ggplot?
  expect_true(inherits(SinglePlot, "gg"))
})
