#' Dashboard Mapper, individual date history retrieval
#'
#' @param x The individual date
#' @param data The instrument data
#' @param Instrument The name of the instrument
#'
#' @importFrom dplyr filter
#'
#' @return Individual date QC summary
#' @noRd
DateMapper <- function(x, data, Instrument){
  DateTime <- x
  TheDateData <- data |> dplyr::filter(DateTime %in% x)
  TheDateSummary <- VisualQCSummary(x=TheDateData)
  InstrumentStatus <- ColorCodeStatus(x=Instrument, y=TheDateSummary)
  Snapshot <- cbind(DateTime, InstrumentStatus)
  Snapshot$DateTime <- as.Date(Snapshot$DateTime)
  return(Snapshot)
}