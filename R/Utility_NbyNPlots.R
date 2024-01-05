Utility_NbyNPlots <- function(x, sample.name, removestrings, experiment = NULL, experiment.name = NULL,  outpath, marginsubset, gatesubset, ycolumn, bins, clearance, sourcelocation){
  #library(flowCore); #library(dplyr)
  x <- x
  name <- keyword(x, sample.name)

  NameCleanUp <- function(name, removestrings){
    for(i in removestrings){
      name <- str_replace_all(name, i, "")
    }
    return(name)
  }

  name <- NameCleanUp(name = name, removestrings = remove_strings)

  #Retrieving Experiment Info #Switched to an exist statement.
  if(!is.null(experiment)){experiment <- experiment
  } else {experiment <- keyword(x, experiment.name)}
  #suppressWarnings(rm(experiment)) #Being Used Somewhere

  AggregateName <- paste(name, experiment, sep = "_")
  StorageLocation <- paste(outpath, AggregateName, sep = "/", collapse = NULL)

  mff <- gs_pop_get_data(x, marginsubset)
  df <- flowCore::exprs(mff[[1]])
  TheDF <- data.frame(df, check.names = FALSE)
  DFNames <- colnames(TheDF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(TheDF))])

  ff <- gs_pop_get_data(x, gatesubset)

  source(sourcelocation, local = TRUE)
}
