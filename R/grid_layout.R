#' Internal, used as basic html framework InteractiveLuciernaga
#' 
#' @param plot_list A list of plot objects
#' @param ncol TBD
#' 
#' @importFrom htmltools div
#' 
#' @return A list of HTML div elements arranged in a grid layout
#' 
#' @noRd
grid_layout <- function(plot_list, ncol = 2) {
  rows <- split(plot_list, ceiling(seq_along(plot_list) / ncol))
  html <- lapply(rows, function(row_plots) {
    div(style = "display: flex; justify-content: center;",
        lapply(row_plots, function(p) {
          div(style = "width: 50%; padding: 10px;", p)
        })
    )
  })
  return(html)
}