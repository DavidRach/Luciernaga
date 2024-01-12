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
#' @importFrom flowCore keyword
#' @importFrom flowCore exprs
#'
#' @return NULL
#' @export
#'
#' @examples NULL

Utility_NbyNPlots <- function(x, sample.name, removestrings, experiment = NULL, experiment.name = NULL, condition = NULL, condition.name = NULL,
                              marginsubset, gatesubset, ycolumn, bins, clearance, gatelines, reference = NULL, outpath){
  ycolumn <- ycolumn
  x <- x
  name <- flowCore::keyword(x, sample.name)

  name <- NameCleanUp(name = name, removestrings)

  if(!is.null(experiment)){experiment <- experiment
  } else {experiment <- flowCore::keyword(x, experiment.name)}

  if(!is.null(condition)){condition <- condition
  } else {condition <- flowCore::keyword(x, condition.name)}

  AggregateName <- paste(name, experiment, sep = "_") #Additional for condition (we need to think this through)
  StorageLocation <- paste(outpath, AggregateName, sep = "/", collapse = NULL)

  mff <- gs_pop_get_data(x, marginsubset)
  df <- flowCore::exprs(mff[[1]])
  TheDF <- data.frame(df, check.names = FALSE)
  DFNames <- colnames(TheDF[,-grep("Time|FS|SC|SS|Original|W$|H$", names(TheDF))])

  ff <- gs_pop_get_data(x, gatesubset)

  if (ycolumn == "ALL"){print("run ComprehensiveNxNPlot fctn")
  } else {
    columnlist <- DFNames[DFNames != ycolumn]
    Plots <- map(.x = columnlist, .f = Utility_GeneralGating, name = name, ff = ff, yValue = ycolumn, columnlist = DFNames,
                 TheDF = TheDF, gatelines = gatelines, reference = reference, clearance, bins)
  }
  return(Plots)
}









