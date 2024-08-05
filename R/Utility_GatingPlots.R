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
#' @param clearance A buffer area around the plot edge
#' @param gtFile The data.table imported .csv file containing the gating template.
#' @param DesiredGates A vector of gates that you want plotted, for example Desired <- c("nonDebris, "lymphocytes")
#' @param export Whether to return as a .pdf file or as individual patchwork ggplot objects.
#' @param thecolumns Number of desired columns for the .pdf file
#' @param therows Number of desired rows for the .pdf file
#' @param width Desired page width
#' @param height Desired page height
#'
#' @param outpath Location to store the generated .pdf file
#'
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr pull
#' @importFrom purrr map
#'
#' @return Additional information to be added
#' @export
#'
#' @examples NULL

Utility_GatingPlots <- function(x, sample.name, removestrings, subsets="root", gtFile, DesiredGates = NULL,
                                outpath = NULL, export = TRUE, ...){

  # Setting up individual file name
  if (is.null(outpath)){outpath <- getwd()}

  AggregateName <- NameForSample(x=x, sample.name=sample.name, removestrings=removestrings, ...)
  # AggregateName <- Luciernaga:::NameForSample(x=x, sample.name=sample.name, removestrings=removestrings)

  # Pulling Gating Information
  TheXYZgates <- gtFile %>% pull(alias)

  # Desired Gate
  if(!is.null(DesiredGates)){TheXYZgates <- intersect(DesiredGates, TheXYZgates)}

  # Pulling Gating Set Data
  ff <- gs_pop_get_data(x, subsets)
  df <- exprs(ff[[1]])
  TheDF <- data.frame(df, check.names = FALSE)
  x2 <- x

  #Plot Generation
  CompiledPlots <- map(.x = TheXYZgates, .f = GatePlot, data=x2, gtFile = gtFile)
  # CompiledPlots <- map(.x = TheXYZgates, .f = Luciernaga:::GatePlot, data=x2, gtFile = gtFile)

  if (export == TRUE){
    AssembledPlots <- Utility_Patchwork(x=CompiledPlots, filename=AggregateName, outfolder=outpath, returntype = "pdf")
  } else if (export == FALSE){
    AssembledPlots <- Utility_Patchwork(x=CompiledPlots, filename=AggregateName, outfolder=outpath, returntype = "patchwork")}

  return(AssembledPlots)
}
