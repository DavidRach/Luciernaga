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
#' @keywords internal
#'
#' @examples NULL

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
