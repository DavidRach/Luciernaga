test_that("Luciernaga_GroupHeatmap returns a ggplot2 object", {

  # NonFunctional #TODOLIST #MultiSpecimen

  expect_true(length(MyGatingSet) > 0)

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
