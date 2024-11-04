test_that("Luciernaga_Plots returns a ggplot2 object", {

  # Prepare the test

  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(data.table)
  library(dplyr)
  library(purrr)

  File_Location <- system.file("extdata", package = "Luciernaga")
  FCS_Files <- list.files(path = File_Location, pattern = ".fcs", full.names = TRUE)
  UnstainedFCSFiles <- FCS_Files[grep("Unstained", FCS_Files)]
  UnstainedCells <- UnstainedFCSFiles[-grep(
    "Beads", UnstainedFCSFiles)]
  MyCytoSet <- load_cytoset_from_fcs(UnstainedCells[c(1,3,5)],
                                     truncate_max_range = FALSE,transformation = FALSE)
  MyGatingSet <- GatingSet(MyCytoSet)
  MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
  MyGatingTemplate <- gatingTemplate(MyGates)
  gt_gating(MyGatingTemplate, MyGatingSet)
  removestrings <-  c("DR_", "Cells", ".fcs", "-", " ")
  StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
  FileLocation <- system.file("extdata", package = "Luciernaga")
  pattern = "AutofluorescentOverlaps.csv"
  AFOverlap <- list.files(path=FileLocation, pattern=pattern, full.names = TRUE)

  reports <- map(.x=MyGatingSet[1:3], .f=Luciernaga_QC, subsets="lymphocytes",
                          removestrings=removestrings, sample.name="GROUPNAME",
                          unmixingcontroltype = "cells", Unstained = TRUE,
                          ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
                          stats = "median", ExportType = "data", SignatureReturnNow = FALSE,
                          outpath = TemporaryFolder, Increments=0.1, SecondaryPeaks=2,
                          experiment = "Lymphocytes", condition = "Ctrl",
                          Subtraction = "Internal", SCData="subtracted",
                          NegativeType="default")

  # Execute the test

  plot <- Luciernaga_GroupHeatmap(reports=reports, nameColumn="Sample", cutoff=0.02, returntype = "plot")


  # Did it return a data.frame
  expect_true(inherits(plot, "gg"))
})
