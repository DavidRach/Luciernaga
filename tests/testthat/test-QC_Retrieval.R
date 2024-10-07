test_that("QC_Retrieval returns a dataframe with more than 1 row", {
  library(flowCore)
  library(flowWorkspace)

  # Prepare the experiment
  File_Location <- system.file("extdata", package = "Luciernaga")
  FCS_Pattern <- ".fcs$"
  FCS_Files <- list.files(path = File_Location, pattern = FCS_Pattern,
                          full.names = TRUE, recursive = FALSE)
  QCBeads <- FCS_Files[grep("Before", FCS_Files)]
  MyCytoSet <- load_cytoset_from_fcs(QCBeads[1], truncate_max_range = FALSE,
                                     transform = FALSE)
  MyGatingSet <- GatingSet(MyCytoSet)

  # Execute the experiment
  SingleSpecimen <- QC_Retrieval(x=MyGatingSet[[1]], sample.name="TUBENAME")

  # Did it return a data.frame
  expect_s3_class(SingleSpecimen, "data.frame")

  # Is it more than 1 row
  expect_gt(nrow(SingleSpecimen), 0)
})
