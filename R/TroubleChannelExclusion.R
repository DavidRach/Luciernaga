#' Internal for LuciernagaQC
#'
#' @param x Vector of fluorophore names to exclude
#' @param TheSCData Data frame containing spectral unmixing matrix information
#' @param MainDetector Column name in `TheSCData` containing main detector IDs
#' @param AFChannels Vector of autofluorescence channels
#'
#' @importFrom dplyr filter pull
#'
#' @return A character vector of excluded channels.
#'
#' @noRd
TroubleChannelExclusion <- function(x, TheSCData, MainDetector, AFChannels) {

  Internal <- TheSCData |>
    filter(Fluorophore %in% x) |>
    pull(MainDetector)
  Internal <- gsub("-A", "", Internal)
  Exclusion <- setdiff(AFChannels, Internal)

  return(Exclusion)
}