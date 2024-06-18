#' Check gate placement for individual .fcs files
#'
#' @param x A Gating Set object
#' @param sample.name The .fcs keyword containing the samples name
#' @param removestrings The values you want removed from sample.name
#' @param experiment Provide an alternative experiment name
#' @param experiment.name The .fcs keyword containing the experiment value
#' @param condition Provide an alternative condition name
#' @param condition.name The .fcs keyword containing the condition name
#' @param subsets The subset level you want the data plotted at. "root" as default
#' @param gtFile The imported data.table containing the gate information.
#' @param column Residual superceeded
#' @param bins Number of geom hexbins to bin the plotted data for visualization.
#' @param clearance Area around edge
#' @param outpath The location to send the generated .pdf file
#' @param export Whether to export as a .pdf file to the outpath, TRUE, FALSE
#'
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom stringr str_split
#' @importFrom ggcyto as.ggplot
#' @importFrom ggcyto ggcyto
#' @import ggplot2
#' @importFrom purrr map
#' @importFrom patchwork wrap_plots
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom stats quantile
#'
#' @return Additional information to be added
#' @export
#'
#' @examples NULL

Utility_GatingPlots <- function(x, sample.name, removestrings,
    experiment = NULL, experiment.name = NULL, condition = NULL,
    condition.name = NULL, subsets, gtFile, column, bins, clearance,
    outpath, export = NULL){

  name <- keyword(x, sample.name)
  name <- NameCleanUp(name = name, removestrings)

  if(exists("experiment")) {experiment <- experiment
  } else if (exists("experiment.name")) {experiment <- keyword(x, experiment.name)
  } else {experiment <- NULL}

  if(exists("condition")){condition <- condition
  } else if (exists("condition.name")) {condition <- keyword(x, condition.name)
  } else {condition <- NULL}

  if (!is.null(experiment)){
    AggregateName <- paste0(name, experiment)
  } else {AggregateName <- name}

  StorageLocation <- file.path(outpath, AggregateName)

  ff <- gs_pop_get_data(x, subsets)
  df <- exprs(ff[[1]]) #Is the one necessary in this case? Unclear
  #how works with lapply...
  TheDF <- data.frame(df, check.names = FALSE)

  TheXYZgates <- gtFile %>% pull(alias)

  x2 <- x

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

  CompiledPlots <- map(.x = TheXYZgates, .f = GatePlot, gtFile = gtFile)

  theList <- CompiledPlots
  theListLength <- length(theList)

  thecolumns <- 2
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
