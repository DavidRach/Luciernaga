
#' Internal for RetrieveAveragedSignature, iterates on per gate
#'  basis to grab the corresponding averaged values for the cells
#'  contained within
#' 
#' @param template The template data.frame
#' @param gs A gating set
#' @param inverse.transform Whether to inverse.transform, default
#'  is set to TRUE
#' @param excludeThese Removes columns whose name contain these
#'  values via str_detect, default set for "FSC|SSC|Time|-H" 
#' which is typical for Cytek Aurora raw .fcs files
#' 
#' @importFrom flowWorkspace gs_pop_get_data 
#' @importFrom flowCore exprs
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate relocate
#' 
#' @return A data.frame row containing the fluorophore name and
#'  the corresponding detector columns with the average signature
#'  value for the cells within that gate
#' 
#' @export
#' 
#' @examples A <- 2+2
#' 
GateExprsIterated <- function(gs, sample, gate,
   inverse.transform=TRUE, excludeThese="FSC|SSC|Time|-H"){

  PopulationInterest <- gs_pop_get_data(gs[[sample]],
     gate, inverse.transform=inverse.transform)
  TheDataValues <- exprs(PopulationInterest[[1]])

  TheDataValues <- data.frame(TheDataValues, check.names=FALSE)

  DetectorsOnly <- TheDataValues[, !stringr::str_detect(
    names(TheDataValues), excludeThese)]
  
  Signature <- AveragedSignature(DetectorsOnly, stats="median",
   normalize=FALSE)

  Signature <- Signature |> mutate(Fluorophore=gate) |> 
    relocate(Fluorophore, .before=1)

  return(Signature)
  }