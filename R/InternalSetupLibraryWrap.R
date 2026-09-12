
#' Internal for SetupLibraryPlotWrapper
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
#' @importFrom purrr map
#' 
#' @return A ggplot2 object
#' 
#' @noRd
InternalSetupLibraryWrap <- function(x, distinguish, data, columnname,
 detectorcolumn, valuecolumn, Normalize, TheFormat){

  SubsetData <- data |> filter(.data[[distinguish]] %in% x)
  TheTarget <- SubsetData |> pull(columnname) |> unique()

  InnerPlots <- map(.x=TheTarget, .f=QC_ViewSignature, data=SubsetData,
   columnname = columnname, detectorcolumn=detectorcolumn,
    valuecolumn=valuecolumn, Normalize=Normalize, TheFormat=TheFormat)

   return(InnerPlots)
}