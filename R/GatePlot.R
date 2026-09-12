#' Generates called plots from Utility_GatingPlots
#'
#' @param x A specific gate, ex. "nonDebris"
#' @param data A GatingSet object
#' @param TheDF A data.frame object of the flow file's expr data
#' @param gtFile The data.table imported .csv file containing the gating template.
#' @param bins Argument to geom_hex for number of bins to visualize the plotted
#' data density.
#' @param clearance A buffer area around the plot edge
#' @param name Sets the title for the plot, default is NULL
#'
#' @importFrom dplyr filter pull select
#' @importFrom stringr str_split
#' @importFrom tidyselect all_of
#' @importFrom ggcyto ggcyto as.ggplot geom_gate
#' @importFrom ggplot2 geom_hex theme_bw aes labs theme element_blank
#'  element_line element_text coord_cartesian
#' @importFrom flowWorkspace gs_pop_get_data cytoframe_to_flowFrame
#' @importFrom stats quantile

#'
#' @return A ggplot corresponding to the given inputs
#'
#' @noRd
GatePlot <- function(x, data, TheDF, gtFile, bins=270, clearance = 0.2,
  name){
    i <- x
    gtFile <- data.frame(gtFile, check.names = FALSE)
    RowData <- gtFile |> filter(alias %in% i)
    theSubset <- RowData |> pull(parent)
    theGate <- RowData |> pull(alias)
    theParameters <- RowData |> pull(dims) |>
      str_split(",", simplify = TRUE)

    theParameters <- gsub("^\\s+|\\s+$", "", theParameters)

    if(length(theParameters) == 2){xValue <- theParameters[[1]]
    yValue <- theParameters[[2]]
    } else if (length(theParameters) == 1){
      xValue <- theParameters[[1]]
      yValue <- "SSC-A" #or an alternate variable specify
    } else {message(
    "Plotting Parameters for Axis were not 1 or 2, please check the .csv file")
    }


  #Please Note, All the Below Are Raw Values With No Transforms Yet Applied.

  if (!grepl("FSC|SSC", xValue)) {

  if (!xValue %in% colnames(TheDF)){
    internal_cs <- gs_pop_get_data(data)
    ff <- cytoframe_to_flowFrame(internal_cs[[1]])
    Workaround <- ff@parameters@data
    xValue <- Workaround |> filter(desc %in% xValue) |> pull(name)
  }
    
  ExprsData <- TheDF |> select(all_of(xValue)) |> pull()
  theXmin <- ExprsData %>% quantile(., 0.001)
  theXmax <- ExprsData %>% quantile(., 0.999)
  theXmin <- theXmin - abs((clearance*theXmin))
  theXmax <- theXmax + (clearance*theXmax)}

  if (!grepl("FSC|SSC", yValue)) {
  
  if (!yValue %in% colnames(TheDF)){
    internal_cs <- gs_pop_get_data(data)
    ff <- cytoframe_to_flowFrame(internal_cs[[1]])
    Workaround <- ff@parameters@data
    yValue <- Workaround |> filter(desc %in% yValue) |> pull(name)
  }

  ExprsData <- TheDF |>select(all_of(yValue)) |> pull()
  theYmin <- ExprsData %>% quantile(., 0.001)
  theYmax <- ExprsData %>% quantile(., 0.999)
  theYmin <- theYmin - abs((clearance*theYmin))
    theYmax <- theYmax + (clearance*theYmax)}

  if (!exists("theYmax") || !exists("theXmax")){
    Plot <- ggcyto(data, aes(x = .data[[xValue]], y = .data[[yValue]]),
       subset = theSubset) + geom_hex(bins=bins) + geom_gate(theGate) + theme_bw() +
       labs(title = name) + theme(strip.background = element_blank(),
       strip.text.x = element_blank(), panel.grid.major = element_line(
       linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"),
       axis.title = element_text(size = 10, face = "bold"), legend.position = "none")
    Plot <- as.ggplot(Plot)
  } else {
    Plot <- as.ggplot(ggcyto(data, aes(x = .data[[xValue]], y = .data[[yValue]]), subset = theSubset)) +
      geom_hex(bins=bins) +
      coord_cartesian(xlim = c(theXmin, theXmax), ylim = c(theYmin, theYmax), default = TRUE) +
      geom_gate(theGate) + theme_bw() + labs(title = name) +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank(),
            panel.grid.major = element_line(linetype = "blank"),
            panel.grid.minor = element_line( linetype = "blank"),
            axis.title = element_text(size = 10),
            legend.position = "none")
  }
}
