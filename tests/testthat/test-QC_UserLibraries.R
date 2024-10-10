test_that("QC_UserLibraries returns a pdf object", {
  # Prepare the experiment
  library(purrr)
  library(dplyr)
  StorageLocation <- file.path(tempdir(), "LuciernagaTemporaryExamples")
  if (!dir.exists(StorageLocation)) {dir.create(StorageLocation)}

  Folder_Location <- system.file("extdata", package = "Luciernaga")
  XML_Pattern <- ".XML$"
  XML_Files <- list.files(path = Folder_Location, pattern = XML_Pattern,
                          full.names = TRUE, recursive = FALSE)
  Data <- map(.x=XML_Files[1:4], .f=QC_LibraryParse,
              returntype="dataframe", references=FALSE) %>% bind_rows()
  TheIndividuals <- Data %>% pull(Creator) %>% unique()

  # Execute the experiment
  JohnDoesLibrary <- QC_UserLibraries(x=TheIndividuals[[1]], Data=Data,
                                      NameAppend="_LibraryQC",
                                      outpath=StorageLocation,
                                      references = TRUE, thecolumns = 3,
                                      therows=4, width=7, height=9,
                                      saveCSV=FALSE)

  ThePDF <- list.files(StorageLocation, pattern="_LibraryQC.pdf")

  # Did it return a pdf?
  expect_true(length(ThePDF) >= 1)
})
