test_that("Luciernaga_QC returns a dataframe with at least 1 row", {

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
  MyCytoSet <- load_cytoset_from_fcs(CellSingleColors[1:3],
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

  # Execute the test
  SingleColor_Data <- map(.x=MyGatingSet[1:3], .f=Luciernaga_QC, subsets="lymphocytes",
                          removestrings=removestrings, sample.name="GUID",
                          unmixingcontroltype = "cells", Unstained = FALSE,
                          ratiopopcutoff = 0.001, Verbose = FALSE, AFOverlap = AFOverlap,
                          stats = "median", ExportType = "data", SignatureReturnNow = FALSE,
                          outpath = TemporaryFolder, Increments=0.1, SecondaryPeaks=2,
                          experiment = "FirstExperiment", condition = "ILTPanel",
                          Subtraction = "Internal", CellAF=TheCellAF, SCData="subtracted",
                          NegativeType="default", RetainedType="normalized") %>% bind_rows()

  # Did it return a data.frame
  expect_s3_class(SingleColor_Data, "data.frame")

  # Is it more than 1 row
  expect_gt(nrow(SingleColor_Data), 0)
})
