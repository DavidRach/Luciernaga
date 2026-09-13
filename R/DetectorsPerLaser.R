#' Internal FrankensteinsConfig, determines number of detectors per laser
#'
#' @param x The iterated laser
#' @param data The vector of detectors
#'
#' @importFrom stringr str_detect
#'
#' @return TBD
#'
#' @noRd
DetectorsPerLaser <- function(x, data) {
  Lasers <- x
  if (Lasers == "V") {
    Lasers <- paste0("^", Lasers)
    Conditional <- TRUE
  } else {
    Conditional <- FALSE
  }

  DetectorCount <- length(data[str_detect(data, Lasers)])

  if (Conditional == TRUE) {
    Lasers <- "V"
  }

  Stats <- data.frame(Lasers, DetectorCount, check.names = FALSE)
  return(Stats)
}