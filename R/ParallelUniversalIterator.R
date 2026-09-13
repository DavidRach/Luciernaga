#' Internal for Utility_ParallelNbyNPlots
#'
#' @param x TBD
#' @param x_ff TBD
#' @param y_ff TBD
#' @param TheDF TBD
#' @param yValue TBD
#' @param columnlist TBD
#' @param gatelines TBD
#' @param reference TBD
#' @param clearance TBD
#' @param bins TBD
#' @param AltNameX TBD
#' @param AltNameY TBD
#' @param colorX TBD
#' @param colorY TBD
#'
#' @importFrom purrr map flatten
#'
#' @return An internal value
#'
#' @noRd
ParallelUniversalIterator <- function(x, x_ff, y_ff,
                                      TheDF, yValue, columnlist, gatelines,
                                      reference, clearance, bins, AltNameX,
                                      AltNameY, colorX, colorY) {
  DFNames <- columnlist
  columnlist <- columnlist[columnlist != x] # Remove the universal Y value

  Plots <- map(.x = columnlist, .f = ParallelGating, x_ff = x_ff,
    y_ff = y_ff, TheDF = TheDF, yValue = x, columnlist = DFNames,
    gatelines = gatelines, reference = reference, clearance = clearance,
    bins = bins, AltNameX = AltNameX, AltNameY = AltNameY, colorX = colorX,
    colorY = colorY) #Name
  Plots <- flatten(Plots)

  #Plots1 <- Plots
  #Plots <- flatten(Plots)
  return(Plots)
}