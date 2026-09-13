#' Internal for QCHistoryArchive
#' 
#' @param x Something
#' @param data Something
#' 
#' @importFrom dplyr filter pull bind_rows
#' @importFrom purrr map
#' 
#' @return A value of some form
#' 
#' @noRd
InternalColorCodeStatus <- function(x, data) {
  TheInstrument <- x
  TheSubset <- data |> filter(Instrument %in% TheInstrument)
  TheDates <- data |> pull(Date) |> unique()

  TheInstrumentHistory <- map(.x = TheDates, .f = InternalColorDateFilter,
                               TheInstrument = TheInstrument,
                               TheSubset = TheSubset) |>
    bind_rows()
  return(TheInstrumentHistory)
}