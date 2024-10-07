test_that("QC_Plots returns a ggplot2 object", {
  library(ggplot2)

  # Prepare the experiment
  File_Location <- system.file("extdata", package = "Luciernaga")
  CSV_Pattern <- ".CSV$"
  CSV_Files <- list.files(path = File_Location, pattern = CSV_Pattern, full.names = TRUE)

  TidyData <- QC_FilePrep(CSV_Files, TrackChange = FALSE)

  # Execute the experiment
  StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")
  SinglePlot <- "UV7-Gain"

  TheSinglePlot <- QC_Plots(
    x = TidyData,
    MeasurementType = SinglePlot,
    FailedFlag = TRUE,
    returntype = "patchwork",
    path = StorageLocation,
    filename = "CytekAurora5L_QC"
  )

  # Did it return a ggplot?
  expect_true(inherits(TheSinglePlot[[1]], "gg"))
})
