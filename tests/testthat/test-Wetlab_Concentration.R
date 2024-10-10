test_that("Wetlab_Concentration returns a dataframe >= 1 row", {
  # Prepare the test
  library(flowCore)
  library(flowWorkspace)
  library(CytoML)
  library(dplyr)
  library(purrr)

  File_Location <- system.file("extdata", package = "Luciernaga")
  WSP_File <- list.files(File_Location, pattern=".wsp", full.names = TRUE)
  ws <- open_flowjo_xml(WSP_File[1])
  gs <- flowjo_to_gatingset(ws, name=1, path = File_Location)
  nameKeyword <- c("GROUPNAME", "TUBENAME")

  # Execute the experiment
  TheData <- map(.x=gs, Wetlab_Concentration, subset = "CD45+",
                 nameKeyword=nameKeyword, DilutionMultiplier=100,
                 TotalVolume=1) %>% bind_rows()

  # Did it return a data.frame
  expect_s3_class(TheData, "data.frame")

  # Is it more than 1 row
  expect_gt(nrow(TheData), 0)
})
