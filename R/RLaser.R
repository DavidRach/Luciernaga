#' Generates a R laser template
#' 
#' @return An R laser data.frame
#' 
#' @importFrom tibble tibble
#' 
#' @noRd
RLaser <- function() {
  Laser <- tibble(
    Detector = c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"),
    Wavelength = c(662, 679, 695, 718, 738, 760, 781, 812)
  )
  return(Laser)
}