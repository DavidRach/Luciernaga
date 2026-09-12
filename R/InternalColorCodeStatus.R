#' Internal for QCHistoryArchive
#' 
#' @param x Something
#' @param TheInstrument Something
#' @param TheSubset Something
#' 
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' 
#' @return A value of some form
#' 
#' @noRd
InternalColorCodeStatus <- function(x, data){
  TheInstrument <- x
  TheSubset <- data |> filter(Instrument %in% TheInstrument)
  TheDates <- data |> pull(Date) |> unique()

  TheInstrumentHistory <- map(.x=TheDates, .f=InternalColorDateFilter,
     TheInstrument=TheInstrument, TheSubset=TheSubset) |> bind_rows()
  return(TheInstrumentHistory)
}