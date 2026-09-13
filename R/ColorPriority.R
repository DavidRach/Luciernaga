#' Dashboard Internal, rearranges vector by color
#'
#' @param colors A vector of cytometer parameters to be rearranged
#'
#' @return Reordered vector according to light wavelength
#' @noRd
ColorPriority <- function(colors) {

  Ordered <- colors[order(grepl("^Ultra", colors) * -1,
                          grepl("^UV", colors) * -1,
                          grepl("^Violet", colors) * -1,
                          grepl("^Blue", colors) * -1,
                          grepl("^Yellow", colors) * -1,
                          grepl("^yellow", colors) * -1,
                          grepl("^Red", colors) * -1)]
  return(Ordered)
}