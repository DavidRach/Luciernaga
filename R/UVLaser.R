#' Generates a UV laser template
#' 
#' @return An UV laser data.frame
#' 
#' @noRd
UVLaser <- function(){
 Laser <- tibble(
  Detector=c("UV1", "UV2", "UV3", "UV4", "UV5", "UV6", "UV7", "UV8",
  "UV9", "UV10", "UV11", "UV12","UV13","UV14","UV15","UV16"),
  Wavelength=c(373, 388, 428, 443, 458, 473, 514, 
  542, 582, 613, 662, 695, 718, 750, 781, 812)
 )
 return(Laser)
}