test_that("QC_ViewSignature returns a ggplot2 object", {

  expect_true(length(MyGatingSet) > 0)

  StorageLocation <- file.path("C:", "Users", "JohnDoe", "Desktop")

  PopulationInterest <- gs_pop_get_data(MyGatingSet[1], subset="lymphocytes")
  TheDataValues <- exprs(PopulationInterest[[1]])
  TheDataValues <- data.frame(TheDataValues, check.names=FALSE)
  Signature <- AveragedSignature(TheDataValues, stats="median")
  TheData1 <- Signature[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Signature))]
  TheData1 <- TheData1 |> mutate(Sample="lymphocytes") |>
    relocate(Sample, .before=1)

  PopulationInterest <- gs_pop_get_data(MyGatingSet[1], subset="nonDebris")
  TheDataValues <- exprs(PopulationInterest[[1]])
  TheDataValues <- data.frame(TheDataValues, check.names=FALSE)
  Signature <- AveragedSignature(TheDataValues, stats="median")
  TheData2 <- Signature[,-grep("Time|FS|SC|SS|Original|W$|H$", names(Signature))]
  TheData2 <- TheData2 |> mutate(Sample="nonDebris") |>
    relocate(Sample, .before=1)

  FinalData <- rbind(TheData1, TheData2)

  # Execute the test

  Plot <- Luciernaga_Cosine(data=FinalData, returntype="plot")

  # Is it more than 1 row
  expect_true(inherits(Plot, "gg"))
})
