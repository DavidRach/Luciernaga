test_that("Luciernaga_FCSToReport returns a dataframe with at least 1 row", {

  tmp <- withr::local_tempdir(pattern = "Luciernaga")
  withr::local_dir(tmp)
  StorageLocation <- tmp

  expect_true(length(MySCsGatingSet) > 0)

  FileLocation <- system.file("extdata", package = "Luciernaga")
  pattern = "AutofluorescentOverlaps.csv"
  AFOverlap <- list.files(path=FileLocation, pattern=pattern, full.names = TRUE)

  SingleColor_Data <- map(.x=MySCsGatingSet[1:2], .f=Luciernaga_QC, subsets="lymphocytes",
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
