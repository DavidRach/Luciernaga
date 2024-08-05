#' Internal for Utility_NbyNPlots
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom patchwork wrap_plots
#' @importFrom patchwork plot_spacer
#' @importFrom purrr map
#'
#' @keywords internal

UniversalIterator <- function(x, x_ff,
                               TheDF, yValue, columnlist, gatelines,
                               reference, clearance, bins, AltNameX,
                               AltNameY, colorX, colorY){
  ff <- x_ff
  DFNames <- columnlist
  columnlist <- columnlist[columnlist != x] # Remove the universal Y value

  Plots <- map(.x = columnlist, .f = Utility_GeneralGating, name = name, ff = ff, yValue = x, columnlist = DFNames,
               TheDF = TheDF, gatelines = gatelines, reference = reference, clearance=clearance, bins=bins)

  #Plots <- flatten(Plots)
  #Plots1 <- Plots
}

#' Internal for Utility_ParallelNbyNPlots
#' @importFrom flowWorkspace keyword
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowCore exprs
#' @importFrom patchwork wrap_plots
#' @importFrom patchwork plot_spacer
#' @importFrom purrr map
#'
#' @keywords internal

ParallelUniversalIterator <- function(x, x_ff, y_ff,
                               TheDF, yValue, columnlist, gatelines,
                               reference, clearance, bins, AltNameX,
                               AltNameY, colorX, colorY){
  DFNames <- columnlist
  columnlist <- columnlist[columnlist != x] # Remove the universal Y value

  Plots <- map(.x = columnlist, .f = Internal_ParallelGating, x_ff=x_ff, y_ff=y_ff,
               TheDF=TheDF, yValue=x, columnlist=DFNames, gatelines=gatelines,
               reference=reference, clearance=clearance, bins=bins, AltNameX=AltNameX,
               AltNameY=AltNameY, colorX=colorX, colorY=colorY) #Name
  Plots <- flatten(Plots)

  #Plots1 <- Plots
  #Plots <- flatten(Plots)
  return(Plots)
}
