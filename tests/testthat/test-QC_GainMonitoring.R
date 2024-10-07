test_that("QC_GainMonitoring returns a dataframe with more than 1 row", {
  library(flowCore)
  library(flowWorkspace)
  library(purrr)

  # Prepare the experiment
  File_Location <- system.file("extdata", package = "Luciernaga")
  FCS_Pattern <- ".fcs$"
  FCS_Files <- list.files(path = File_Location, pattern = FCS_Pattern,
                          full.names = TRUE, recursive = FALSE)
  QCBeads <- FCS_Files[grep("Before|After", FCS_Files)]
  BeforeAfter_CS <- load_cytoset_from_fcs(files=QCBeads,
                                          transform=FALSE,
                                          truncate_max_range = FALSE)

  # Execute the experiment
  BeforeAfter <- map(.x=BeforeAfter_CS[1:2], .f=QC_GainMonitoring,
                     sample.name = "TUBENAME", stats="median") %>% bind_rows()

  # Did it return a data.frame
  expect_s3_class(BeforeAfter, "data.frame")

  # Is it more than 1 row
  expect_gt(nrow(BeforeAfter), 1)
})
