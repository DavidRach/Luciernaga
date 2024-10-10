test_that("Lucierna_Move relocated the .fcs files", {
  # Prepare the test
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(data.table)
  library(dplyr)
  library(purrr)
  library(stringr)

  StorageLocation <- file.path(tempdir(), "LuciernagaMoveExample1")
  if (!dir.exists(StorageLocation)) {dir.create(StorageLocation)}
  File_Location <- system.file("extdata", package = "Luciernaga")
  FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
                          full.names = TRUE)
  CellSingleColorFiles <- FCS_Files[grep("Cells", FCS_Files)]
  CellSingleColors <- CellSingleColorFiles[!str_detect("Unstained", CellSingleColorFiles)]
  MyCytoSet <- load_cytoset_from_fcs(CellSingleColors,
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

  TheCellAF  <- Luciernaga_QC(x=MyUnstainedGatingSet[[1]], subsets="nonDebris",
                              removestrings=removestrings, sample.name="GUID",
                              unmixingcontroltype = "cells", Unstained = TRUE,
                              ratiopopcutoff = 0.001, Verbose = FALSE,
                              AFOverlap = AFOverlap, stats = "median",
                              ExportType = "fcs", Brightness=TRUE,
                              SignatureReturnNow = TRUE, outpath = StorageLocation,
                              Increments=0.1, SecondaryPeaks=2,
                              experiment = "FirstExperiment", condition = "ILTPanel",
                              Subtraction = "Internal", CellAF=TheCellAF,
                              SCData="subtracted",NegativeType="default",
                              minimalfcscutoff=0.01)

  MyGatingSet <- subset(MyGatingSet, !str_detect(name, "CD14"))

  SingleColor_Data <- map(.x=MyGatingSet, .f=Luciernaga_QC,
                          subsets="lymphocytes", removestrings=removestrings,
                          sample.name="GUID", unmixingcontroltype = "cells",
                          Unstained = FALSE, ratiopopcutoff = 0.001,
                          Verbose = FALSE, AFOverlap = AFOverlap,
                          stats = "median", ExportType = "fcs", Brightness=TRUE,
                          SignatureReturnNow = FALSE,outpath = StorageLocation,
                          Increments=0.1, SecondaryPeaks=2,
                          experiment = "FirstExperiment", condition = "ILTPanel",
                          Subtraction = "Internal", CellAF=TheCellAF,
                          SCData="subtracted",NegativeType="default",
                          minimalfcscutoff=0.01)

  Unstained_Data  <- map(.x=MyUnstainedGatingSet[1], .f=Luciernaga_QC,
                         subsets="nonDebris", removestrings=removestrings,
                         sample.name="GUID", unmixingcontroltype = "cells",
                         Unstained = TRUE, ratiopopcutoff = 0.001,
                         Verbose = FALSE, AFOverlap = AFOverlap,
                         stats = "median", ExportType = "fcs", Brightness=TRUE,
                         SignatureReturnNow = FALSE, outpath = StorageLocation,
                         Increments=0.1, SecondaryPeaks=2,
                         experiment = "FirstExperiment", condition = "ILTPanel",
                         Subtraction = "Internal", CellAF=TheCellAF,
                         SCData="subtracted",NegativeType="default",
                         minimalfcscutoff=0.01)

  PanelPath <- file.path(File_Location, "UnmixingPanel.csv")

  # Execute the test

  MoveThese <- Luciernaga_Tree(BrightnessFilePath = StorageLocation, PanelPath = PanelPath)
  SortedStorageLocation <- file.path(tempdir(), "LuciernagaMoveSelected1")
  if (!dir.exists(SortedStorageLocation)) {dir.create(SortedStorageLocation)}
  UnmixingPanel <- read.csv(PanelPath, check.names=FALSE)
  TheseFluorophores <- UnmixingPanel %>% pull(Fluorophore)

  walk(.x=TheseFluorophores, .f=Luciernaga_Move, data=MoveThese,
       input=StorageLocation, output=SortedStorageLocation)
  MovedFiles <- list.files(SortedStorageLocation, pattern="fcs", full.names=TRUE)
  Moved_Length <- length(MovedFiles)

  # Did it return a pdf?
  expect_true(nrow(MoveThese) == Moved_Length)
})
