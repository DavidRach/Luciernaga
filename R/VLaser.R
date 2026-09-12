#' Generates a V laser template
#' 
#' @return An V laser data.frame
#' 
#' @importFrom tibble tibble
#' 
#' @noRd
VLaser <- function() {
  Laser <- tibble(
    Detector = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8",
                 "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16"),
    Wavelength = c(428, 443, 458, 473, 508, 525, 542, 581, 598, 615, 662,
                   695, 718, 750, 781, 812)
  )
  return(Laser)
}