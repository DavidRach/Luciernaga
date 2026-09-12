#' Internal for LuciernagaQC
#'
#' @importFrom dplyr filter pull
#'
#' @return An internal value
#'
#' @noRd
TroubleChannelExclusion <- function(x, TheSCData, MainDetector, AFChannels){
  Internal <- TheSCData %>% filter(Fluorophore %in% x) %>% pull(MainDetector)
  Internal <- gsub("-A", "", Internal)
  Exclusion <- setdiff(AFChannels, Internal)
  return(Exclusion)
}