
#' Internal for Luciernaga_Lists
#'
#' @param i Passed Indicies
#' @param data The Transposed List of Lists
#'
#' @importFrom purrr flatten
#' @importFrom plotly ggplotly
#' @importFrom htmltools tagList div
#'
#' @return An internal value
#'
#' @keywords internal
Subplots <- function(i, data) {
  Components <- flatten(data[i])

  plotlyobjs <- lapply(Components, ggplotly)

  subplot <- tagList(
    lapply(plotlyobjs, function(x) {
      div(
        x,
        style = "float:left; width:50%;",
        tags$br()
      )
    }),
    tags$br(style = "clear:both;"),
    tags$br(), tags$br(), tags$br()
  )

  return(subplot)
}
