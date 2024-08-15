#' Check gate placement for individual .fcs files in a GatingSet
#'
#' @param x A GatingSet object
#' @param sample.name The .fcs keyword that contains an unique name for that sample
#' @param removestrings A string of character values to remove from the sample.name
#' @param experiment.name The .fcs keyword that contains the name of the experiment
#' @param experiment Directly provide the name of the experiment (alternative to experiment.name)
#' @param condition.name The .fcs keyword that contains the name of the condition.
#' @param condition Directly provide the name of the condition (alternative to condition.name)
#' @param subset The GatingSet subset that you want to visualize for data plotting, "root" is the default.
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

Utility_GatingPlots <- function(x, sample.name, removestrings, subset="root", gtFile, DesiredGates = NULL,
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
  ff <- gs_pop_get_data(x, subset)
  df <- exprs(ff[[1]])
  TheDF <- data.frame(df, check.names = FALSE)
  x2 <- x

  #Plot Generation
  CompiledPlots <- map(.x = TheXYZgates, .f = GatePlot, data=x2, TheDF = TheDF, gtFile = gtFile)
  # CompiledPlots <- map(.x = TheXYZgates, .f = Luciernaga:::GatePlot, data=x2, gtFile = gtFile)

  if (export == TRUE){
    AssembledPlots <- Utility_Patchwork(x=CompiledPlots, filename=AggregateName, outfolder=outpath, returntype = "pdf")
  } else if (export == FALSE){
    AssembledPlots <- Utility_Patchwork(x=CompiledPlots, filename=AggregateName, outfolder=outpath, returntype = "patchwork")}

  return(AssembledPlots)
}

#' Generates called plots from Utility_GatingPlots
#'
#' @param x A specific gate, ex. "nonDebris"
#' @param data A GatingSet object
#' @param TheDF A data.frame object of the flow file's expr data
#' @param gtFile The data.table imported .csv file containing the gating template.
#' @param bins Argument to geom_hex for number of bins to visualize the plotted data density.
#' @param clearance A buffer area around the plot edge
#'
#' @importFrom dplyr filter
#' @importFrom stringr str_split
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom ggcyto ggcyto
#' @importFrom ggcyto as.ggplot
#' @import ggplot2
#'
#' @return A ggplot corresponding to the given inputs
#' @noRd
GatePlot <- function(x, data, TheDF, gtFile, bins=270, clearance = 0.2){
  i <- x
  gtFile <- data.frame(gtFile, check.names = FALSE)
  RowData <- gtFile %>% filter(alias %in% i)
  theSubset <- RowData %>% pull(parent)
  theGate <- RowData %>% pull(alias)
  theParameters <- RowData %>% pull(dims) %>% str_split(",", simplify = TRUE)

  theParameters <- gsub("^\\s+|\\s+$", "", theParameters)

  if(length(theParameters) == 2){xValue <- theParameters[[1]]
  yValue <- theParameters[[2]]
  } else if (length(theParameters) == 1){xValue <- theParameters[[1]]
  yValue <- "SSC-A" #or an alternate variable specify
  } else {message(
    "Plotting Parameters for Axis were not 1 or 2, please check the .csv file")}


  #Please Note, All the Below Are Raw Values With No Transforms Yet Applied.

  if (!grepl("FSC|SSC", xValue)) {ExprsData <- TheDF %>%
    select(all_of(xValue)) %>% pull()
  theXmin <- ExprsData %>% quantile(., 0.001)
  theXmax <- ExprsData %>% quantile(., 0.999)
  theXmin <- theXmin - abs((clearance*theXmin))
  theXmax <- theXmax + (clearance*theXmax)}
  if (!grepl("FSC|SSC", yValue)) {ExprsData <- TheDF %>%
    select(all_of(yValue)) %>% pull()
  theYmin <- ExprsData %>% quantile(., 0.001)
  theYmax <- ExprsData %>% quantile(., 0.999)
  theYmin <- theYmin - abs((clearance*theYmin))
  theYmax <- theYmax + (clearance*theYmax)}

  if (!exists("theYmax") || !exists("theXmax")){
    Plot <- as.ggplot(ggcyto(data, aes(x = .data[[xValue]], y = .data[[yValue]]),
                             subset = theSubset) + geom_hex(bins=bins) + geom_gate(theGate) +
                        theme_bw() + labs(title = NULL) + theme(
                          strip.background = element_blank(), strip.text.x = element_blank(),
                          panel.grid.major = element_line(linetype = "blank"),
                          panel.grid.minor = element_line(linetype = "blank"),
                          axis.title = element_text(size = 10, face = "bold"),
                          legend.position = "none"))

  } else {Plot <- as.ggplot(ggcyto(data, aes(
    x = .data[[xValue]], y = .data[[yValue]]), subset = theSubset) +
      geom_hex(bins=bins) + coord_cartesian(xlim = c(theXmin, theXmax),
                                            ylim = c(theYmin, theYmax), default = TRUE) + geom_gate(theGate) +
      theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(),
                                              strip.text.x = element_blank(), panel.grid.major = element_line(
                                                linetype = "blank"), panel.grid.minor = element_line(
                                                  linetype = "blank"), axis.title = element_text(size = 10,
                                                                                                 face = "bold"),
                                              legend.position = "none"))

  }
}


