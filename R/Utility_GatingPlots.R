#' Check gate placement for individual .fcs files in a GatingSet
#'
#' @param x A GatingSet object
#' @param sample.name The .fcs keyword that contains an unique name for that sample
#' @param removestrings A string of character values to remove from the sample.name
#' @param experiment.name The .fcs keyword that contains the name of the experiment
#' @param experiment Directly provide the name of the experiment (alternative to experiment.name)
#' @param condition.name The .fcs keyword that contains the name of the condition.
#' @param condition Directly provide the name of the condition (alternative to condition.name)
#' @param subsets The GatingSet subset that you want to visualize for data plotting, "root" is the default.
#' @param bins Argument to geom_hex for number of bins to visualize the plotted data density.
#' @param gtFile The data.table imported .csv file containing the gating template.
#' @param DesiredGates A vector of gates that you want plotted, for example Desired <- c("nonDebris, "lymphocytes")
#' @param export Whether to return as a .pdf file or as individual patchwork ggplot objects.
#' @param thecolumns Number of desired columns for the .pdf file
#' @param therows Number of desired rows for the .pdf file
#' @param outpath Location to store the generated .pdf file
#'
#' @param clearance Area around edge
#'
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr pull
#' @importFrom purrr map
#'

#' @importFrom patchwork wrap_plots
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom stats quantile
#'
#' @return Additional information to be added
#' @export
#'
#' @examples NULL

Utility_GatingPlots <- function(x, sample.name, removestrings, subsets="root", gtFile,
    DesiredGates = NULL, thecolumn = 2, therows = 2, bins, clearance, outpath = NULL, export = NULL, ...){

  # Setting up individual file name
  if (is.null(outpath)){outpath <- getwd()}

  AggregateName <- NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, ...)
  StorageLocation <- file.path(outpath, AggregateName)

  # Pulling Gating Information
  TheXYZgates <- gtFile %>% pull(alias)

  # Desired Gate
  TheXYZgates <- intersect(DesiredGates, TheXYZgates)

  # Pulling Gating Set Data
  ff <- gs_pop_get_data(x, subsets)
  df <- exprs(ff[[1]])
  TheDF <- data.frame(df, check.names = FALSE)
  x2 <- x

  #Plot Generation
  CompiledPlots <- map(.x = TheXYZgates, data=x2, .f = GatePlot, gtFile = gtFile)



  # Send to Utility_Patchwork
  theList <- CompiledPlots
  theListLength <- length(theList)

  thecolumns
  therows <- 3
  theoreticalitems <- therows*thecolumns

  split_list <- function(input_list, chunk_size) {
    split(input_list, ceiling(seq_along(input_list) / chunk_size))
  }

  sublists <- split_list(theList, theoreticalitems)

  AssembledPlots <- map(.x = sublists, .f = wrap_plots, ncol = thecolumns,
                        nrow = therows, widths = 0.8, heights = 0.8)

  if (export == TRUE){
  MergedName <- StorageLocation

  pdf(file = paste(MergedName, ".pdf", sep = "", collapse = NULL),
      width = 9, height = 7) #Optional Adjustments for Second

  print(AssembledPlots)

  dev.off()

  } else {return(AssembledPlots)}
}
