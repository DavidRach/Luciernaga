

#' Internal to Utility_RidgePlots
#'
#' @param x Passed Argument for X axis
#' @param cs The CytoSet with the data
#' @param TheY The passed argument for the Y axis
#' @param TheFill The passed argument for the Fill color by factor.
#'
#' @importFrom ggcyto ggcyto
#' @importFrom ggcyto axis_x_inverse_trans
#' @importFrom ggcyto as.ggplot
#' @importFrom ggridges geom_density_ridges
#' @importFrom ggplot2 facet_null
#'
#' @return An internal value
#'
#' @noRd
RidgePlots <- function(x, cs, TheY, TheFill){
  p <- ggcyto(cs, aes(x = .data[[x]]))
  p1 <- p + geom_density_ridges(aes(y = .data[[TheY]], fill = .data[[TheFill]], alpha = 0.2)) +
    facet_null() + theme_bw() + coord_cartesian(xlim = c(0, NA)) + theme(legend.position = "none") +
    axis_x_inverse_trans()
  p1 <- as.ggplot(p1)
  return(p1)
}