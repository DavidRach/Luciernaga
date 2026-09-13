#' Internal for InstrumentComparison, filters for given fluorophores
#' info on a per instrument basis
#'
#' @param x The number detectors corresponding to the instrument
#' @param fluorophores The vector of fluorophores to filter for
#'
#' @importFrom dplyr filter
#'
#' @return Data for the fluorophores
#'
#' @noRd
InstrumentData <- function(x, fluorophores) {
  Data <- InstrumentReferences(NumberDetectors = x)
  Data <- Data |> filter(Fluorophore %in% fluorophores)
  return(Data)
}