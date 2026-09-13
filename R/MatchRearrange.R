#' Internal for InstrumentComparison
#'
#' @param x The instrument being filtered for
#' @param thematch The established laser ordered fluorophores
#' @param data The passed data for all instruments filtering from
#'
#' @importFrom dplyr filter
#'
#' @return Rearranged data.frame with laser ordered rows
#'
#' @noRd
MatchRearrange <- function(x, thematch, data) {

  Subset <- data |> filter(Instrument %in% x)
  MatchedReferences <- match(thematch, Subset$Fluorophore)
  MatchedReferences <- Subset[MatchedReferences, ]
  return(MatchedReferences)
}