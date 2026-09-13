#' Dashboard Internal, rearranges vector by scatter parameter
#'
#' @param colors A vector of scatter parameters to be rearranged
#'
#' @return Rearranged vector of scatter parameters
#' @noRd
ScalePriority <- function(colors) {
  Ordered <- colors[order(grepl("^FSC", colors) * -1,
                           grepl("^SSC", colors) * -1,
                           grepl("-A", colors) * -1,
                           grepl("-H", colors) * -1,
                           grepl("-W", colors) * -1)]
  return(Ordered)
}