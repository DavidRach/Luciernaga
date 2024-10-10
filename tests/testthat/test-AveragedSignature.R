test_that("AveragedSignature returns a dataframe with more than 1 row", {
  # Prepare the test

  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(data.table)

  File_Location <- system.file("extdata", package = "Luciernaga")
  FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
                          full.names = TRUE)
  UnstainedFCSFiles <- FCS_Files[grep("Unstained", FCS_Files)]
  UnstainedCells <- UnstainedFCSFiles[-grep("Beads", UnstainedFCSFiles)]
  MyCytoSet <- load_cytoset_from_fcs(UnstainedCells[1],
                                     truncate_max_range = FALSE,
                                     transformation = FALSE)
  MyGatingSet <- GatingSet(MyCytoSet)
  MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
  MyGatingTemplate <- gatingTemplate(MyGates)
  gt_gating(MyGatingTemplate, MyGatingSet)
  removestrings <-  c("DR_", "Cells", ".fcs", "-", " ")
  StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")

  PopulationInterest <- gs_pop_get_data(MyGatingSet[1], subset="lymphocytes")
  TheDataValues <- exprs(PopulationInterest[[1]])
  TheDataValues <- data.frame(TheDataValues, check.names=FALSE)

  # Execute the test

  Signature <- AveragedSignature(TheDataValues, stats="median")

  # Is it more than 1 row
  expect_equal(nrow(Signature), 1)
})
