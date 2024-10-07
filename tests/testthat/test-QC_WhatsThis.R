test_that("QC_WhatsThis returns a dataframe with more than 1 row", {

  # Prepare the experiment
  Folder_Location <- system.file("extdata", package = "Luciernaga")
  XML_Pattern <- ".XML$"
  XML_Files <- list.files(path = Folder_Location, pattern = XML_Pattern,
                          full.names = TRUE, recursive = FALSE)
  LibraryData <- QC_LibraryParse(XML_Files[2], returntype="dataframe", references=FALSE)
  Data <- LibraryData %>% select(-Sample, -Creator, -Date) %>%
    rename(Sample=Fluorochrome)
  TheFluorophore <- Data %>% pull(Sample)

  # Execute the experiment
  Results <- QC_WhatsThis(x=TheFluorophore, data=Data, NumberHits = 10, returnPlots=FALSE)

  # Did it return a data.frame
  expect_s3_class(Results, "data.frame")

  # Is it more than 1 row
  expect_gt(nrow(Results), 1)
})
