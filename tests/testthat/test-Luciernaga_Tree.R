test_that("Luciernaga_Tree returns a dataframe with at least 1 row", {

  # Prepare the test
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(data.table)
  library(dplyr)
  library(purrr)
  library(stringr)

  StorageLocation <- file.path(tempdir(), "LuciernagaTreeExample")
  if (!dir.exists(StorageLocation)) {dir.create(StorageLocation)}
  File_Location <- system.file("extdata", package = "Luciernaga")
  FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
                          full.names = TRUE)
  CellSingleColorFiles <- FCS_Files[grep("Cells", FCS_Files)]
  CellSingleColors <- CellSingleColorFiles[!str_detect("Unstained", CellSingleColorFiles)]
  MyCytoSet <- load_cytoset_from_fcs(CellSingleColors[1:2],
                                     truncate_max_range = FALSE,
                                     transformation = FALSE)
  MyGatingSet <- GatingSet(MyCytoSet)
  MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
  MyGatingTemplate <- gatingTemplate(MyGates)
  gt_gating(MyGatingTemplate, MyGatingSet)

  UnstainedFCSFiles <- FCS_Files[grep("Unstained", FCS_Files)]
  UnstainedCells <- UnstainedFCSFiles[-grep("Beads", UnstainedFCSFiles)]
  MyUnstainedCytoSet <- load_cytoset_from_fcs(UnstainedCells[1],
                                              truncate_max_range = FALSE,
                                              transformation = FALSE)
  MyUnstainedGatingSet <- GatingSet(MyUnstainedCytoSet)
  gt_gating(MyGatingTemplate, MyUnstainedGatingSet)
  removestrings <-  c(".fcs", "(", ")", "Cells")
  pattern = "AutofluorescentOverlaps.csv"
  AFOverlap <- list.files(path=File_Location, pattern=pattern, full.names = TRUE)

  SingleColor_Data <- map(.x=MyGatingSet[1:2], .f=Luciernaga_QC,
    subsets="lymphocytes", removestrings=removestrings, sample.name="GUID",
    unmixingcontroltype = "cells", Unstained = FALSE, ratiopopcutoff = 0.001,
    Verbose = FALSE, AFOverlap = AFOverlap, stats = "median", ExportType = "fcs",
    Brightness=TRUE, SignatureReturnNow = FALSE,outpath = StorageLocation,
    Increments=0.1, SecondaryPeaks=2, experiment = "FirstExperiment",
    condition = "ILTPanel", Subtraction = "Internal", CellAF=TheCellAF,
    SCData="subtracted",NegativeType="default")

  Unstained_Data  <- map(.x=MyUnstainedGatingSet[1], .f=Luciernaga_QC,
                         subsets="nonDebris", removestrings=removestrings,
                         sample.name="GUID", unmixingcontroltype = "cells",
                         Unstained = TRUE, ratiopopcutoff = 0.001,
                         Verbose = FALSE, AFOverlap = AFOverlap,
                         stats = "median", ExportType = "fcs", Brightness=TRUE,
                         SignatureReturnNow = FALSE, outpath = StorageLocation,
                         Increments=0.1, SecondaryPeaks=2, experiment = "FirstExperiment",
                         condition = "ILTPanel", Subtraction = "Internal", CellAF=TheCellAF,
                         SCData="subtracted",NegativeType="default", minimalfcscutoff=0.01)

  PanelPath <- file.path(File_Location, "UnmixingPanel.csv")

  # Execute the test

  MoveThese <- Luciernaga_Tree(BrightnessFilePath = StorageLocation, PanelPath = PanelPath)

  # Did it return a data.frame
  expect_s3_class(MoveThese, "data.frame")

  # Is it more than 1 row
  expect_gt(nrow(MoveThese), 0)
})
