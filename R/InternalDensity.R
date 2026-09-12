#' Internal for Utility_DensityOverlay
#'
#' @param x A passed x-axis name
#' @param cs The passed CytoSet
#' @param TheFill The passed factor to fill plots by
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_density
#'
#' @return An internal value
#'
#' @keywords internal
InternalDensity <- function(x, cs, TheFill){
  Plot <- ggplot(cs, aes(x = .data[[x]], fill = .data[[TheFill]])) +
    geom_density(alpha = 0.2) + coord_cartesian(xlim = c(0, NA)) + theme_bw() + theme(
      legend.position = "none")
  return(Plot)
}