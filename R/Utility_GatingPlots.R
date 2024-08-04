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
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom stringr str_split
#' @importFrom ggcyto as.ggplot
#' @importFrom ggcyto ggcyto
#' @import ggplot2
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


  # Pulling Gating Set Data
  ff <- gs_pop_get_data(x, subsets)
  df <- exprs(ff[[1]])
  TheDF <- data.frame(df, check.names = FALSE)

  x2 <- x
  CompiledPlots <- map(.x = TheXYZgates, .f = GatePlot, gtFile = gtFile)

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










GatePlot <- function(x, gtFile){
  i <- x
  gtFile <- data.frame(gtFile, check.names = FALSE)
  RowData <- gtFile %>% dplyr::filter(alias %in% i)
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
    Plot <- as.ggplot(ggcyto(x2, aes(x = .data[[xValue]], y = .data[[yValue]]),
                             subset = theSubset) + geom_hex(bins=bins) + geom_gate(theGate) +
                        theme_bw() + labs(title = NULL) + theme(
                          strip.background = element_blank(), strip.text.x = element_blank(),
                          panel.grid.major = element_line(linetype = "blank"),
                          panel.grid.minor = element_line(linetype = "blank"),
                          axis.title = element_text(size = 10, face = "bold"),
                          legend.position = "none"))

  } else {Plot <- as.ggplot(ggcyto(x2, aes(
    x = .data[[xValue]], y = .data[[yValue]]), subset = theSubset) +
      geom_hex(bins=bins) + coord_cartesian(xlim = c(theXmin, theXmax),
                                            ylim = c(theYmin, theYmax), default = TRUE) + geom_gate(theGate) +
      theme_bw() + labs(title = NULL) + theme(strip.background = element_blank(),
                                              strip.text.x = element_blank(), panel.grid.major = element_line(
                                                linetype = "blank"), panel.grid.minor = element_line(
                                                  linetype = "blank"), axis.title = element_text(size = 10,
                                                                                                 face = "bold"),
                                              legend.position = "none"))

  tryCatch({rm("theXmin", "theXmax", "theYmin", "theYmax")})}
}
