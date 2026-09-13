#' Generates a B laser template
#'
#' @importFrom tibble tibble
#'
#' @return An B laser data.frame
#'
#' @noRd
BLaser <- function() {

  Laser <- tibble(
   Detector=c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8",
   "B9", "B10", "B11", "B12","B13","B14"),
   Wavelength=c(508, 525, 542, 581, 598, 615, 662,
   679, 695, 718, 738, 760, 781, 812)
  )
  return(Laser)
}