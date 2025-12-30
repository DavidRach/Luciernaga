test_that("AveragedSignature returns a dataframe with more than 1 row", {

  expect_true(length(MyGatingSet) > 0)

  PopulationInterest <- gs_pop_get_data(MyGatingSet[1], subset="lymphocytes")
  TheDataValues <- exprs(PopulationInterest[[1]])
  TheDataValues <- data.frame(TheDataValues, check.names=FALSE)

  # Execute the test

  Signature <- AveragedSignature(TheDataValues, stats="median")

  # Is it more than 1 row
  expect_equal(nrow(Signature), 1)
})
