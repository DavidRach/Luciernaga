test_that("Luciernaga_LinearSlices returns a ggplot2 object", {

  # Prepare the test
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(data.table)
  library(dplyr)
  library(purrr)
  library(stringr)

  File_Location <- system.file("extdata", package = "Luciernaga")
  FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
                          full.names = TRUE)
  CellSingleColorFiles <- FCS_Files[grep("Cells", FCS_Files)]
  CellSingleColors <- CellSingleColorFiles[!str_detect("Unstained", CellSingleColorFiles)]

  APC_Example <- CellSingleColors[grep("CD16_", CellSingleColors)]

  MyCytoSet <- load_cytoset_from_fcs(APC_Example,
                                     truncate_max_range = FALSE,
                                     transformation = FALSE)
  MyGatingSet <- GatingSet(MyCytoSet)
  MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
  MyGatingTemplate <- gatingTemplate(MyGates)
  gt_gating(MyGatingTemplate, MyGatingSet)
  removestrings <- c(".fcs")

  # Execute the test

  NormalizedSlices <- Luciernaga_LinearSlices(x=MyGatingSet[1],
                                              subset="lymphocytes",
                                              sample.name="GUID",
                                              removestrings=removestrings,
                                              stats="median",
                                              returntype="normalized",
                                              probsratio=0.1,
                                              output="plot",
                                              desiredAF="R1-A")


  # Did it return a data.frame
  expect_true(inherits(NormalizedSlices, "gg"))
})
