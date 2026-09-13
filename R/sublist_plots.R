#' Wraps the sublist to hand off to patchwork
#'
#' @param x TBD
#' @param thecolumns TBD
#' @param therows TBD
#'
#' @importFrom patchwork wrap_plots
#'
#' @return An internal value
#'
#' @noRd
sublist_plots <- function(x, thecolumns, therows) {
  p <- wrap_plots(x, ncol = thecolumns, nrow = therows, widths = 0.8,
                   heights = 0.8)
}