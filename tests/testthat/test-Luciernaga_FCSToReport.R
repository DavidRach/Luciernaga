test_that("Luciernaga_FCSToReport returns a dataframe with at least 1 row", {

  # Prepare the test
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(data.table)
  library(dplyr)
  library(purrr)
  library(stringr)

  StorageLocation <- file.path(tempdir(), "LuciernagaFCSToReportExample1")
  if (!dir.exists(StorageLocation)) {dir.create(StorageLocation)}
  File_Location <- system.file("extdata", package = "Luciernaga")
  FCS_Files <- list.files(path = File_Location, pattern = ".fcs", full.names = TRUE)
  CellSingleColorFiles <- FCS_Files[grep("Cells", FCS_Files)]
  CellSingleColors <- CellSingleColorFiles[!str_detect("Unstained", CellSingleColorFiles)]
  MyCytoSet <- load_cytoset_from_fcs(CellSingleColors[1:2],
                                     truncate_max_range = FALSE,
                                     transformation = FALSE)
  MyGatingSet <- GatingSet(MyCytoSet)
  MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
  MyGatingTemplate <- gatingTemplate(MyGates)
  gt_gating(MyGatingTemplate, MyGatingSet)
  removestrings <-  c(".fcs", "(", ")", "Cells")
  FileLocation <- system.file("extdata", package = "Luciernaga")
  pattern = "AutofluorescentOverlaps.csv"
  AFOverlap <- list.files(path=FileLocation, pattern=pattern, full.names = TRUE)

  SingleColor_Data <- map(.x=MyGatingSet[1:2], .f=Luciernaga_QC, subsets="lymphocytes",
                          removestrings=removestrings, sample.name="GUID",
                          unmixingcontroltype = "cells", Unstained = FALSE,
                          ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
                          stats = "median", ExportType = "fcs", Brightness=TRUE,
                          SignatureReturnNow = FALSE,outpath = StorageLocation,
                          Increments=0.1, SecondaryPeaks=2, experiment = "FirstExperiment",
                          condition = "ILTPanel", Subtraction = "Internal",
                          CellAF=TheCellAF, SCData="subtracted",NegativeType="default")

  TheLuciernagaOutputs_FCS <- list.files(StorageLocation, pattern="fcs", full.names = TRUE)
  TheLuciernagaOutputs_CSV <- list.files(StorageLocation, pattern="csv", full.names = TRUE)
  PanelPath <- file.path(File_Location, "Panel.csv")

  # Execute the test
  ReportOutput <- Luciernaga_FCSToReport(path=StorageLocation, reference=PanelPath,
                                         stats="median", RetainedType = "normalized",
                                         experiment="FirstExperiment",
                                         condition="ILTExperiment", TheSummary = TRUE)

  # Did it return a data.frame
  expect_s3_class(ReportOutput, "data.frame")

  # Is it more than 1 row
  expect_gt(nrow(ReportOutput), 0)
})
