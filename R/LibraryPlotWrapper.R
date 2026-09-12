#' Internal, thin wrapper that filters for Fluorophore and then
#'  group plots from Library Data
#'
#' @param x A iterated Fluorophore Name
#' @param data The data output from QC_Library
#' @param myfactor The desired factor for group
#' @param animate Whether to convert to ggplotly output,
#'  default FALSE
#' 
#' @importFrom dplyr filter
#' @importFrom plotly ggplotly
#'
#' @return A ggplot2 or a ggplotly object
#' @noRd
LibraryPlotWrapper <- function(x, data, myfactor, animate=FALSE){
  Subset <- data %>% filter(Fluorophore %in% x)
  plot <- Luciernaga:::LibraryPlot(x=Subset, myfactor=myfactor)
  if (animate == TRUE){
    plot <- plotly::ggplotly(plot)
  }
  return(plot)
}
