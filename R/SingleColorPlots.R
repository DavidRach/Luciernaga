#' Visualize Gating Hierarchy for each .fcs file to check the gating.
#'
#' @param x  A Gating Set object (example, gs or gs[[1]])
#' @param subsets The Gating Hierarchy node of cells to be plotted (ex. "root")
#' @param column Unclear
#' @param sample.name Keyword under which samples names are stored (ex. "GUID")
#' @param experiment Provide an experiment name directly. (ex. "Jan2024")
#' @param experiment.name Keyword under which experiment names are stored (ex. "GROUPNAME")
#' @param condition Provide a condition name directly (ex. "PMA")
#' @param condition.name Keyword under which condition names are stored (ex. "TUBENAME")
#' @param bins Bins for which the ggplot object (ex. 270)
#' @param clearance Margin by which the quantile MFI is multiplied by to avoid clipping main pop, but exclude outliers
#' @param outpath Location where the generated .pdf should be stored.
#' @param sourcelocation Location where the .csv for openCyto gating is stored.
#'
#' @return A pdf of the visualized ggplots.
#' @export
#'
#' @examples NULL
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
