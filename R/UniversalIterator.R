#' Internal for Utility_NbyNPlots
#'
#' @param x Parameter x
#' @param x_ff Parameter x_ff
#' @param TheDF Parameter TheDF
#' @param yValue Parameter yValue
#' @param columnlist Parameter columnlist
#' @param gatelines Parameter gatelines
#' @param reference Parameter reference
#' @param clearance Parameter clearance
#' @param bins Parameter bins
#' @param AltNameX Parameter AltNameX
#' @param AltNameY Parameter AltNameY
#' @param colorX Parameter colorX
#' @param colorY Parameter colorY
#'
#' @importFrom flowWorkspace keyword gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom patchwork wrap_plots plot_spacer
#' @importFrom purrr map
#'
#' @return An internal value
#'
#' @noRd
UniversalIterator <- function(x, x_ff,
                              TheDF, yValue, columnlist, gatelines,
                              reference, clearance, bins, AltNameX,
                              AltNameY, colorX, colorY) {

  ff <- x_ff
  DFNames <- columnlist
  columnlist <- columnlist[columnlist != x] # Remove the universal Y value

  Plots <- map(
    .x = columnlist, .f = GeneralGating, name = name, ff = ff,
    yValue = x, columnlist = DFNames, TheDF = TheDF, gatelines = gatelines,
    reference = reference, clearance = clearance, bins = bins
  )

  # Plots <- flatten(Plots)
  # Plots1 <- Plots
}