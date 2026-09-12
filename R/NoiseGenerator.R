#' Internal for Simulated Data, returns vector of noise
#'
#' @param x The data.frame to get dimensions from
#'
#' @importFrom stats rnorm
#'
#' @return A vector of noise values
#' @noRd
NoiseGenerator <- function(x){
  rnorm(nrow(x), mean=0, sd=6)
}