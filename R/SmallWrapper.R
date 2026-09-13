#' Internal for QC_ReferenceLibrary
#' 
#' @param x TBD
#' @param data TBD
#' @param TheseFluorophores TBD
#' @param unstained TBD
#' 
#' @importFrom dplyr filter
#' 
#' @return A ggplot object
#' 
#' @noRd
SmallWrapper <- function(x, data, TheseFluorophores, unstained) {
  Internal <- data |> filter(Instrument %in% x)
  ThePlot <- SimilarFluorPlots(TheseFluorophores = TheseFluorophores,
                                TheFluorophore = NULL, data = Internal,
                                plotlinecolor = plotlinecolor,
                                legend = legend, plotname = plotname,
                                unstained = unstained)
}