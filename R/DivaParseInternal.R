#' Internal for TubeIterate, returns Gains and Metadata for the tube
#' 
#' @param x The iterated xml_node for fluorophore being parsed
#' @param y The name of the iterated fluorophore
#' @importFrom xml2 xml_children xml_name xml_text
#' @importFrom dplyr mutate
#' 
#' @return The iterated data.frame row of Fluorophore and Gain
#' @noRd
DivaParseInternal <- function(x, y){
  TheData <- xml_children(x)
  Fluorophore <- y
  Gain <- TheData[xml_name(TheData) == "voltage"]
  Gain <- xml_text(Gain)

  if (length(Gain) == 0){
      Data <- data.frame(Fluorophore)
      Data <- Data |> mutate(Gain="0")
  } else {
  Data <- data.frame(cbind(Fluorophore, Gain))
  }

  Data$Gain <- as.numeric(Data$Gain)
  
  return(Data)
  }