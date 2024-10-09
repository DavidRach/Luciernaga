test_that("Utility_NbyNPlots returns a patchwork object", {
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
  flowWorkspace::transform(UnmixedGatingSet, TransformList)
  FileLocation <- system.file("extdata", package = "Luciernaga")
  UnmixedGates <- fread(file.path(path = FileLocation, pattern = 'GatesUnmixed.csv'))
  UnmixedGating <- gatingTemplate(UnmixedGates)
  gt_gating(UnmixedGating, UnmixedGatingSet)
  removestrings <-  c("DTR_", ".fcs")
  StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")

  # Execute the experiment
  IndividualNxN <- Utility_NbyNPlots(x=UnmixedGatingSet[[1]],
                                     sample.name = "GROUPNAME",
                                     removestrings = removestrings,
                                     marginsubset = "lymphocytes",
                                     gatesubset = "live",
                                     ycolumn = "Spark Blue 550-A",
                                     bins = 70, clearance = 0.2,
                                     gatelines = FALSE,
                                     reference = NULL,
                                     outpath = StorageLocation,
                                     returntype="patchwork")

  # Did it return a patchwork object?
  expect_true(inherits(IndividualNxN[[1]], "patchwork"))
})
