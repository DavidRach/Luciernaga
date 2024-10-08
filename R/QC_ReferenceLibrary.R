#' Queries the available reference library for available fluorophores (and their naming conventions)
#'
#' @param FluorNameContains A character string pattern to match, example "APC"
#' @param NumberDetectors The Number of Detectors for your instrument
#' @param returnPlots Whether to return signature plot as well. Default FALSE.
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#'
#' @return A dataframe column containing matching Fluorophores from your querry
#' @export
#'
#' @examples
#'
#' QC_ReferenceLibrary(FluorNameContains = "FITC", NumberDetectors=64)
QC_ReferenceLibrary <- function(FluorNameContains, NumberDetectors, returnPlots=FALSE){
  ReferenceData <- InstrumentReferences(NumberDetectors=NumberDetectors)
  if (returnPlots == TRUE){ReferenceData1 <- ReferenceData}

  TheList <- ReferenceData %>% select(Fluorophore) %>% unique()
  rownames(TheList) <- NULL

  Subset <- TheList %>% filter(str_detect(Fluorophore, FluorNameContains))

  if (returnPlots==TRUE){
    TheseFluorophores <- Subset %>% pull(Fluorophore)

    ThePlot <- SimilarFluorPlots(TheseFluorophores=TheseFluorophores,
                                 TheFluorophore=NULL, data=ReferenceData1)
    ReturnThese <- list(Subset, ThePlot)
    return(ReturnThese)
  } else {return(Subset)}
}
