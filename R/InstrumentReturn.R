#' Internal for InstrumentComparison, pulls fluorophore list for
#'  respective instrument
#'
#' @param x Number detectors corresponding desired instrument
#'
#' @importFrom dplyr pull
#'
#' @return The reference fluorophores for that instrument
#'
#' @noRd
InstrumentReturn <- function(x) {
  Data <- InstrumentReferences(NumberDetectors = x)
  Fluorophores <- Data |> pull(Fluorophore) |> unique()
  return(Fluorophores)
}