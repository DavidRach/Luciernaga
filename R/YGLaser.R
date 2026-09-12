#' Generates a YG laser template
#' 
#' @return An YG laser data.frame
#' 
#' @noRd
YGLaser <- function(){
  Laser <- tibble(
   Detector=c("YG1", "YG2", "YG3", "YG4", "YG5", "YG6", "YG7", "YG8",
   "YG9", "YG10"),
   Wavelength=c(577, 598, 615, 662, 679, 695, 718, 750, 781, 812)
  )
  return(Laser)
 }