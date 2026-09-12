
#' Internal for Utility_NbyNPlots
#' @importFrom flowWorkspace keyword gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom patchwork wrap_plots  plot_spacer
#' @importFrom purrr map
#'
#' @return An internal value
#'
#' @noRd
UniversalIterator <- function(x, x_ff,
                              TheDF, yValue, columnlist, gatelines,
                              reference, clearance, bins, AltNameX,
                              AltNameY, colorX, colorY){
  ff <- x_ff
  DFNames <- columnlist
  columnlist <- columnlist[columnlist != x] # Remove the universal Y value

  Plots <- map(.x = columnlist, .f = GeneralGating, name = name, ff = ff,
    yValue = x, columnlist = DFNames, TheDF = TheDF, gatelines = gatelines,
    reference = reference, clearance=clearance, bins=bins)

  #Plots <- flatten(Plots)
  #Plots1 <- Plots
}