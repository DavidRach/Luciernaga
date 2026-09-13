#' Generates called plots from Utility_IterativeGating
#'
#' @param x The name of a individual specimen in a GatingSet
#' @param gs The GatingSet Object
#' @param subset Desired node of cells to plot
#' @param gate Desired gate to show on the same plot as the subset, else NULL
#' @param xValue Desired x axis
#' @param yValue Desired y axis
#' @param sample.name Keyword under which the sample name is stored
#' @param removestrings Character values to remove from the name
#' @param bins Geo_hex bins for the dots
#' @param plotname Default FALSE, adds name to title
#'
#' @importFrom Biobase pData
#' @importFrom flowCore keyword
#' @importFrom BiocGenerics subset
#' @importFrom ggcyto ggcyto as.ggplot geom_gate
#' @importFrom ggplot2 ggplot aes geom_hex labs theme_bw theme element_blank
#'  element_line element_text
#' @importFrom purrr map
#' @importFrom rlang .data
#'
#' @return A ggplot corresponding to the given inputs
#'
#' @noRd
InternalIterator <- function(x,
                              gs,
                              subset,
                              gate,
                              xValue,
                              yValue,
                              sample.name,
                              removestrings,
                              bins,
                              plotname) {

  theGSsubset <- subset(gs, name == x)

  if (length(sample.name) == 2) {
    first <- sample.name[[1]]
    second <- sample.name[[2]]
    first <- keyword(theGSsubset, first)
    second <- keyword(theGSsubset, second)
    name <- paste(first, second, sep = "_")
  } else {
    name <- keyword(theGSsubset, sample.name)
  }

  name <- NameCleanUp(name, removestrings)

  if (!is.null(gate)) {

    Plot <- as.ggplot(ggcyto(theGSsubset,
                              aes(x = .data[[xValue]], y = .data[[yValue]]),
                              subset = subset) +
      geom_gate(gate) +
      geom_hex(bins = bins) +
      labs(title = name) +
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank(),
            panel.grid.major = element_line(linetype = "blank"),
            panel.grid.minor = element_line(linetype = "blank"),
            axis.title = element_text(size = 10, face = "bold"),
            legend.position = "none"))
  } else {
    Plot <- as.ggplot(ggcyto(theGSsubset,
                              aes(x = .data[[xValue]], y = .data[[yValue]]),
                              subset = subset) +
      geom_hex(bins = bins) +
      labs(title = name) +
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank(),
            panel.grid.major = element_line(linetype = "blank"),
            panel.grid.minor = element_line(linetype = "blank"),
            axis.title = element_text(size = 10, face = "bold"),
            legend.position = "none"))
  }

  return(Plot)
}