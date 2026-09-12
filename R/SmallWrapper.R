#' Internal for QC_ReferenceLibrary
#' 
#' @importFrom dplyr filter
#' 
#' @noRd
SmallWrapper <- function(x, data, TheseFluorophores, unstained){
    Internal <- data |> filter(Instrument %in% x)
    ThePlot <- SimilarFluorPlots(TheseFluorophores=TheseFluorophores,
      TheFluorophore=NULL, data=Internal, plotlinecolor=plotlinecolor,
      legend=legend, plotname=plotname, unstained=unstained)
  }