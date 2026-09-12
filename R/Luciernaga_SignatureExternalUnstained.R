#' Wrapper for LuciernagaQC for large scale profiling. Matches unstained
#'  with single colors for respective dates, extracts signature and plots,
#'  returns data.frame and list of patchwork objects.
#' 
#' @param x The list of Dates
#' @param UnstainedList The list of file.paths for the unstained files
#' @param FluorophoreList The list of file.paths for the fluorophore files
#' @param Multiple Default FALSE, if expecting multiple single color
#'  controls per unstained set TRUE
#' @param GateTemplatePath File.path to the openCyto gating template .csv
#' @param removestrings Default is ".fcs", removes from name
#' @param AFOverlap File.path to LuciernagaQC AFOverlap .csv to
#'  handle exceptions. 
#' @param controlType Either "beads" or "cells" (selects respective
#'  external AF protocol)
#' @param subsets The desired openCyto gating population to extract
#'  signature from
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom progressr with_progress progressor
#' 
#' @export
#' 
#' @return A data.frame and a list of patchwork plots
#' 
#' @examples 
#' 
#' library(Luciernaga)
Luciernaga_SignatureExternalUnstained <- function(x, UnstainedList,
   FluorophoreList, Multiple=FALSE, GateTemplatePath,
   removestrings=".fcs", AFOverlap, controlType, subsets){
  
  with_progress({
    p <- progressor(along = x)  
    
    Iterated <- map(.x = x, .f = function(date) {

      result <- ControlWrapper(date, UnstainedList = UnstainedList,
        FluorophoreList = FluorophoreList, Multiple = Multiple,
        GateTemplatePath = GateTemplatePath, AFOverlap = AFOverlap,
        controlType=controlType, subsets=subsets)
      p()
      return(result)
    })
  })

  ThePlots <- lapply(Iterated, function(x) x$Plot)
  TheData <- lapply(Iterated, function(x) x$data)
  TheDataset <- bind_rows(TheData)
  TheFinalList <- list(data=TheDataset, plot=ThePlots)
  return(TheFinalList)
}


