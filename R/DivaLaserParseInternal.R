#' Internal for TubeIterate, returns the iterated laser information
#'
#' @param x The iterrated xml_node with the laser information
#' @param y The name of the laser being iterrated on
#'
#' @importFrom xml2 xml_children xml_name xml_text
#'
#' @return A data.frame row containing the parsed data
#'
#' @noRd
DivaLaserParseInternal <- function(x, y) {
  Laser <- y
  InternalLanding <- xml_children(x)

  Delay <- InternalLanding[xml_name(InternalLanding) == "delay"]
  Delay <- xml_text(Delay)

  AreaScaling <- InternalLanding[xml_name(InternalLanding) == "area_scaling"]
  AreaScaling <- xml_text(AreaScaling)

  Data <- data.frame(cbind(Laser, Delay, AreaScaling))
  return(Data)
}