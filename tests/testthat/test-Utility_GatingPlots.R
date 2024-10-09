test_that("Utility_GatingPlots returns a patchwork object", {
  # Prepare the experiment
  library(flowCore)
  library(flowWorkspace)
  library(openCyto)
  library(data.table)

  File_Location <- system.file("extdata", package = "Luciernaga")
  FCS_Files <- list.files(path = File_Location, pattern = ".fcs",
                          full.names = TRUE)
  UnstainedFCSFiles <- FCS_Files[grep("Unstained", FCS_Files)]
  UnstainedCells <- UnstainedFCSFiles[-grep("Beads", UnstainedFCSFiles)]
  MyCytoSet <- load_cytoset_from_fcs(UnstainedCells[1],
                                     truncate_max_range = FALSE,
                                     transformation = FALSE)
  MyGatingSet <- GatingSet(MyCytoSet)
  MyGates <- fread(file.path(path = File_Location, pattern = 'Gates.csv'))
  MyGatingTemplate <- gatingTemplate(MyGates)
  gt_gating(MyGatingTemplate, MyGatingSet)
  removestrings <-  c("DR_", "Cells", ".fcs", "-", " ")
  StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")

  # Execute the experiment
  IndividualPlot <- Utility_GatingPlots(x=MyGatingSet[[1]],
                                        sample.name = "GUID",
                                        removestrings = removestrings,
                                        gtFile = MyGates,
                                        DesiredGates = NULL,
                                        outpath = StorageLocation,
                                        export = FALSE)

  # Did it return a patchwork object?
  expect_true(inherits(IndividualPlot[[1]], "patchwork"))
})
