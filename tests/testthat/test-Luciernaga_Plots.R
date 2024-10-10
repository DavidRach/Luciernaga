test_that("Luciernaga_Plots returns a ggplot2 object", {

  # Prepare the test

  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(data.table)
  library(dplyr)
  library(purrr)
  #'
  File_Location <- system.file("extdata", package = "Luciernaga")
  FCS_Files <- list.files(path = File_Location, pattern = ".fcs", full.names = TRUE)
  CellSingleColorFiles <- FCS_Files[grep("Cells", FCS_Files)]
  CellSingleColors <- CellSingleColorFiles[!str_detect("Unstained", CellSingleColorFiles)]
  MyCytoSet <- load_cytoset_from_fcs(CellSingleColors[1:2],
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

  SingleColor_Data <- map(.x=MyGatingSet[1:2], .f=Luciernaga_QC, subsets="lymphocytes",
                          removestrings=removestrings, sample.name="GUID",
                          unmixingcontroltype = "cells", Unstained = FALSE,
                          ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
                          stats = "median", ExportType = "data", SignatureReturnNow = FALSE,
                          outpath = TemporaryFolder, Increments=0.1, SecondaryPeaks=2,
                          experiment = "FirstExperiment", condition = "ILTPanel",
                          Subtraction = "Internal", CellAF=TheCellAF, SCData="subtracted",
                          NegativeType="default") %>% bind_rows()

  pattern = "^Panel.csv"
  CSV <- list.files(path=FileLocation, pattern=pattern, full.names=TRUE)
  TheFluorophoreOrder <- read.csv(CSV, check.names = FALSE)

  # Execute the test

  ThePlots <- Luciernaga_Plots(data=SingleColor_Data, RetainedType="normalized",
                               CellPopRatio=0.05, outfolder=NULL,
                               filename="LuciernagaReport", returntype="plots",
                               LinePlots=FALSE, CosinePlots=FALSE,
                               StackedBarPlots = FALSE, HeatmapPlots = TRUE,
                               reference = TheFluorophoreOrder)

  ThePlots <- flatten(ThePlots)

  # Did it return a data.frame
  expect_true(inherits(ThePlots[[1]], "gg"))
})
