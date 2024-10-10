test_that("Utility_ColAppends returns a flowframe", {
  # Prepare the experiment
  library(BiocGenerics)
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(data.table)

  File_Location <- system.file("extdata", package = "Luciernaga")
  FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
                          full.names = TRUE)
  Unmixed_FullStained <- FCS_Files[grep("Unmixed", FCS_Files)]
  UnmixedFCSFiles <- Unmixed_FullStained[1]
  UnmixedCytoSet <- load_cytoset_from_fcs(UnmixedFCSFiles[1],
                                          truncate_max_range = FALSE,
                                          transformation = FALSE)
  UnmixedGatingSet <- GatingSet(UnmixedCytoSet)
  Markers <- colnames(UnmixedCytoSet)
  KeptMarkers <- Markers[-grep("Time|FS|SC|SS|Original|-W$|-H$|AF", Markers)]
  MyBiexponentialTransform <- flowjo_biexp_trans(channelRange = 256,
                                                 maxValue = 1000000,
                                                 pos = 4.5, neg = 0,
                                                 widthBasis = -1000)
  TransformList <- transformerList(KeptMarkers, MyBiexponentialTransform)
  UnmixedGatingSet <- flowWorkspace::transform(UnmixedGatingSet, TransformList)
  FileLocation <- system.file("extdata", package = "Luciernaga")
  UnmixedGates <- fread(file.path(path = FileLocation, pattern = 'GatesUnmixed.csv'))
  UnmixedGating <- gatingTemplate(UnmixedGates)
  gt_gating(UnmixedGating, UnmixedGatingSet)
  removestrings <-  c("DTR_", ".fcs")
  StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")

  ff <- gs_pop_get_data(UnmixedGatingSet[1], subsets="live", inverse.transform = FALSE)
  BeforeParameters <- ff[[1, returnType = "flowFrame"]]
  MainDataFrame <- as.data.frame(exprs(ff[[1]]), check.names = FALSE)
  NewData <- MainDataFrame %>% mutate(ExposureStatus = sample(1:3, n(), replace = TRUE))
  NewData <- NewData %>% select(ExposureStatus)

  # Execute the experiment

  AfterParameters <- Utility_ColAppend(ff=ff, DF=MainDataFrame, columnframe = NewData)

  # Did it return a patchwork object?
  expect_true(inherits(AfterParameters, "flowFrame"))
})
