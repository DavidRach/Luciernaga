test_that("Utility_GatingPlots returns a patchwork object", {

  expect_true(length(MyUnstainedGatingSet) > 0)

  StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")

  IndividualPlot <- Utility_GatingPlots(x=MyUnstainedGatingSet[[1]],
                                        sample.name = "GUID",
                                        removestrings = removestrings,
                                        gtFile = MyGates,
                                        DesiredGates = NULL,
                                        outpath = StorageLocation,
                                        returnType = "patchwork")
  
  # Debugging
  x <- MyUnstainedGatingSet[[1]]
  sample.name <- "GUID"
  gtFile <- MyGates
  outpath <- StorageLocation
  returnType <- "patchwork"
  plotname <- FALSE
  DesiredGates <- NULL
  optionalName <- NULL
  subset<-"root"
  bins <- 70

  # Did it return a patchwork object?
  expect_true(inherits(IndividualPlot[[1]], "patchwork"))
})
