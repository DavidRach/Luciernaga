#' Internal for Chorizo, setup Display parameters for Description
#'
#' @param x The iterated rowname
#' @param data The Display data.frame to be filtered
#'
#' @importFrom dplyr pull
#'
#' @return The list of display parameters for the detector
#' @noRd
DisplayInternal <- function(x, data){
  Subset <- data[x,]

  TheX <- gsub("$", "", fixed=TRUE, x)

  Display <- paste0(TheX, "DISPLAY")
  DisplayVal <- Subset |> pull(Display)

  DisplayList <- list(Display=DisplayVal)
  names(DisplayList) <- Display
  return(DisplayList)
}