SingleColorPlots <- function(x, subsets, column, sample.name, experiment = NULL, experiment.name = NULL, condition = NULL, condition.name = NULL, bins, clearance, outpath, sourcelocation){
  #library(flowCore); #library(dplyr)
  x <- x
  name <- keyword(x, sample.name)
  name <- gsub(".fcs$", "", name)
  #name <- gsub(" ", "", name)

  experiment <- keyword(x, experiment.name)

  AggregateName <- paste(experiment, name, sep = "_")
  StorageLocation <- paste(outpath, AggregateName, sep = "/", collapse = NULL)

  ff <- gs_pop_get_data(x, subsets)
  df <- flowCore::exprs(ff[[1]]) #Is the one necessary in this case? Unclear how works with lapply...
  TheDF <- data.frame(df, check.names = FALSE)
  source(sourcelocation, local = TRUE)
}
