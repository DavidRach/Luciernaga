#' Internal for QCHistoryArchive
#' 
#' @param x Something
#' @param TheInstrument Something
#' @param TheSubset Something
#' 
#' @importFrom dplyr filter
#' 
#' @return A value of some form
#' 
#' @noRd
InternalColorDateFilter <- function(x, TheInstrument, TheSubset){
  DateTime <- as.Date(x)
  TheDateSummary <- TheSubset |> filter(Date %in% DateTime)
  InstrumentStatus <- ColorCodeStatus(x=TheInstrument, y=TheDateSummary)
  Snapshot <- cbind(DateTime, InstrumentStatus)
  Snapshot$DateTime <- as.Date(Snapshot$DateTime)
  return(Snapshot)
}