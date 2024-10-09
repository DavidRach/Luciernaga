test_that("Utility_IterativeGating returns a ggplot2 object", {
  # Prepare the experiment
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

  # Execute the experiment
  SingleSpecimen <- Utility_IterativeGating(x=MyGatingSet[1],
                                            sample.name = "GUID",
                                            removestrings = removestrings,
                                            subset = "nonDebris",
                                            gate = "lymphocytes",
                                            xValue = "FSC-A",
                                            yValue = "SSC-A",
                                            bins = 270)

  # Did it return a ggplot2 object?
  expect_true(inherits(SingleSpecimen[[1]], "gg"))
})
