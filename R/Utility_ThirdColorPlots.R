Utility_ThirdColorPlots <- function(x, subset, xaxis, yaxis, zaxis, sample.name, removestrings){

  AggregateName <- NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, ...)
  #StorageLocation <- file.path(outpath, AggregateName)

  # Retrieving margin info for the x specimen
  Margin <- gs_pop_get_data(x, subset)
  df <- exprs(Margin[[1]])
  TheDF <- data.frame(df, check.names = FALSE)



}
