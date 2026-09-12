#' Imports Detector, GainBaseline and RCVCutoff for non-Cytek instruments
#'
#' @param x  A data.frame or file.path to the .csv with the Detector, GainBaseline and RCVCutoff
#'
#' @importFrom utils read.csv
#' 
#' @return Data for subsequent use
#' 
#' @noRd
NotCytekDailyQC <- function(x){
  if (!is.data.frame(x)){
      Data <- read.csv(x, check.names = FALSE)
  } else {Data <- x}

  TheseColumns <- c("Detector", "GainBaseline", "RCVCutoff")

  if (any(!TheseColumns %in% colnames(Data))){
      stop("CSV file should contain `Detector`, `GainBaseline` and `RCVCutoff` as column names")
  }

  Data$Detector <- as.character(Data$Detector)
  Data$GainBaseline <- as.numeric(Data$GainBaseline)
  Data$RCVCutoff <- as.numeric(Data$RCVCutoff)

  return(Data)
}