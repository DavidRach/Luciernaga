#' Generate NxN plots for one or all fluorophores.
#'
#' @param x A GatingSet object (ex. gs or gs[[1]])
#' @param sample.name Keyword for which samples names are stored (ex. "GUID")
#' @param removestrings A string of characters to remove from name (ex. c("DR_ILT_2023_", "Cord"))
#' @param experiment Provide directly an experiment name (ex. "Jan2024")
#' @param experiment.name Keyword for which experiment name is stored (ex. "GROUPNAME")
#' @param outpath Location for which to store the generated .pdf
#' @param marginsubset The Gating Hierarchy level from which values will be used to estimate the plot margins (ex. "nonDebris")
#' @param gatesubset The Gating Hierarchy level of the cells that you want to see plotted (ex. "lymph")
#' @param ycolumn The ycolumn that you want to see everything plotted by (ex. "APC-A") or ALL to see all comparisons
#' @param bins Bins for which the plotted cells will be divided into providing granularity
#' @param clearance The additional ratio added to the margins to avoid clipping main population but exclude outliers.
#' @param sourcelocation Location of the file for individual or all NxN plotting.
#'
#' @return NULL
#' @export
#'
#' @examples NULL
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
