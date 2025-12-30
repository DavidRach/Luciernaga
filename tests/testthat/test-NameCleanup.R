test_that("Name gets cleaned up", {
  
  name <- "DR BUV496 CD8 (Cells).fcs"
  removestrings <- c("DR", "(Cells)", ".fcs", " ")

  expected_outcome <- "BUV496CD8"

  Cleaned_Name <- NameCleanUp(name, removestrings)

  expect_equal(Cleaned_Name, expected_outcome)
})
