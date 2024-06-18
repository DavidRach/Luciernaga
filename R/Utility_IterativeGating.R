#' Take a Gating Set, and return individual plots for specified node
#'
#' @param x A Gating Set Object
#' @param subset Desired node of cells to plot
#' @param gate Desired gate to show on the same plot as the subset, else NULL
#' @param xValue Desired x axis
#' @param yValue Desired y axis
#' @param sample.name Keyword under which the sample name is stored
#' @param remove.strings Character values to remove from the name
#' @param bins Geo_hex bins for the dots
#'
#' @importFrom Biobase pData
#' @importFrom flowCore keyword
#' @importFrom ggcyto ggcyto
#' @importFrom ggcyto as.ggplot
#' @importFrom ggplot2 ggplot
#' @importFrom purrr map
#'
#' @return A list of ggplot objects
#' @export
#'
#' @examples NULL

Utility_IterativeGating <- function(x, subset, gate, xValue, yValue, sample.name, remove.strings, bins){
  gs <- x
  TheSpecimens <- pData(gs)$name

  ThePlots <- map(.x=TheSpecimens, .InternalIterator, subset=subset, gate=gate, xValue=xValue, yValue=yValue, sample.name=sample.name,
                  remove.strings=remove.strings, bins=bins)

  return(ThePlots)
}


.InternalIterator <- function(x, subset, gate, xValue, yValue, sample.name, remove.strings, bins){
  theGSsubset <- subset(gs, name == x)
  name <- keyword(theGSsubset, sample.name)
  name <- NameCleanUp(name, remove.strings)

  if(!is.null(gate)){

  Plot <- as.ggplot(ggcyto(theGSsubset, aes(x=.data[[xValue]], y=.data[[yValue]]), subset=subset) + geom_gate(gate) + geom_hex(bins=bins) +
    labs(title = name) + theme_bw() + theme(strip.background = element_blank(), strip.text.x = element_blank(),
                                            panel.grid.major = element_line(linetype = "blank"),
                                            panel.grid.minor = element_line(linetype = "blank"),
                                            axis.title = element_text(size = 10, face = "bold"),
                                            legend.position = "none"))
  } else {
    Plot <- as.ggplot(ggcyto(theGSsubset, aes(x=.data[[xValue]], y=.data[[yValue]]), subset=subset) + geom_hex(bins=bins) +
    labs(title = name) + theme_bw() + theme(strip.background = element_blank(), strip.text.x = element_blank(),
                                            panel.grid.major = element_line(linetype = "blank"),
                                            panel.grid.minor = element_line(linetype = "blank"),
                                            axis.title = element_text(size = 10, face = "bold"),
                                            legend.position = "none"))
  }


  return(Plot)
}
