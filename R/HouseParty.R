#' Internal for Simulated Data, generates noise data.frame
#'
#' @param x The data.frame we wish to generate noise for
#'
#' @importFrom dplyr bind_cols
#'
#' @return An equivalent data.frame of noise values
#' @noRd
HouseParty <- function(x){
  TheNames <- colnames(x)
  TheParticipants <- list()
  for (i in seq_along(x)){
    TheParticipants[[i]] <- NoiseGenerator(x)
  }
  names(TheParticipants) <- TheNames
  TheParticipants <- bind_cols(TheParticipants)
  return(TheParticipants)
}