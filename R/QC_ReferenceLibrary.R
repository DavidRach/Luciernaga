#' Queries the available reference library for available fluorophores (and their naming conventions)
#'
#' @param FluorNameContains A character string pattern to match, example "APC"
#' @param NumberDetectors The number of detectors the instrument contains.
#' For Cytek Instruments 5L = 64, 4L_UV = 54, 4L_YG = 48, 3L=38, 2L_VB=30,
#' 2L_BR=22, 1L=14
#' For BD S8=78, S6="48_S", A5="48_A"
#' For Sony ID7000 7L=184, 6L_DUV="182_DUV", 5L=147, 4L=112, 3L=86
#' For ThermoFisher BigFoot 7L_488-561=55, 7L_532-594="52_7L", 6L_445="52_6L", 6L_785=51 
#' @param returnPlots Whether to return signature plot as well. Default FALSE.
#' @param ListOverride Default FALSE
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#'
#' @return A dataframe column containing matching Fluorophores from your querry
#' @export
#'
#' @examples
#' QC_ReferenceLibrary(FluorNameContains = "FITC", NumberDetectors=64)
QC_ReferenceLibrary <- function(FluorNameContains, NumberDetectors,
                                returnPlots=FALSE, ListOverride=FALSE){
  ReferenceData <- InstrumentReferences(NumberDetectors=NumberDetectors)
  if (returnPlots == TRUE){ReferenceData1 <- ReferenceData}


  TheList <- ReferenceData %>% select(Fluorophore) %>% unique()
  rownames(TheList) <- NULL

  if (ListOverride == FALSE){
  Subset <- TheList %>% filter(str_detect(Fluorophore, FluorNameContains))
  } else {Subset <- TheList %>% dplyr::filter(
    Fluorophore %in% FluorNameContains)}

  if (returnPlots==TRUE){
    TheseFluorophores <- Subset %>% pull(Fluorophore)

    ThePlot <- SimilarFluorPlots(TheseFluorophores=TheseFluorophores,
                                 TheFluorophore=NULL, data=ReferenceData1)
    ReturnThese <- list(Subset, ThePlot)
    return(ReturnThese)
  } else {return(Subset)}
}
