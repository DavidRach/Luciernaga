
#' A wrapper for taking Reference Library Signatures derrived from SetupLog and generate
#' plots for each of them
#' 
#' @param x A
#' @param distinguish A
#' @param data A
#' @param columnname A
#' @param detectorcolumn A
#' @param valuecolumn A
#' @param Normalize A
#' @param TheFormat A
#'  
#' @importFrom dplyr filter pull
#' @importFrom purrr flatten map
#' 
#' @return A ggplot2 object
#' 
#' @noRd 
SetupLibraryPlotWrapper <- function(x, distinguish, data, columnname, detectorcolumn, 
  valuecolumn, Normalize, TheFormat){
  
  SubsetData <- data |> filter(.data[[columnname]] %in% x)
  Distinguisher <- SubsetData |> pull(distinguish) |> unique()

  ThePlots <- map(.x=Distinguisher, .f=InternalSetupLibraryWrap, distinguish=distinguish,
   data=SubsetData, columnname = columnname, detectorcolumn=detectorcolumn,
    valuecolumn=valuecolumn, Normalize=Normalize, TheFormat=TheFormat)

  ThePlots <- flatten(ThePlots)

  return(ThePlots)
}