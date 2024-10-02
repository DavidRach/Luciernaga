#' Queries the available reference library for available fluorophores (and their naming conventions)
#'
#' @param FluorNameContains A character string pattern to match, example "APC"
#' @param NumberDetectors The Number of Detectors for your instrument
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#'
#' @return A dataframe column containing matching Fluorophores from your querry
#' @export
#'
#' @examples NULL
QC_ReferenceLibrary <- function(FluorNameContains, NumberDetectors){
  ReferenceData <- InstrumentReferences(NumberDetectors=NumberDetectors)
  TheList <- ReferenceData %>% select(Fluorophore) %>% unique()
  rownames(TheList) <- NULL

  Subset <- TheList %>% filter(str_detect(Fluorophore, FluorNameContains))
  return(Subset)
}
