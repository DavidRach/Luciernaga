#' Dashboard Internal, wrapper for individual instrument history
#'
#' @param x The Instrument name
#' @param y The Instrument data
#' @param timewindow The number desired months
#'
#' @importFrom dplyr filter pull bind_rows
#' @importFrom purrr map
#'
#' @return Individual instrument QC history summary
#' @noRd
AcrossTime <- function(x, y, timewindow) {

  WindowOfInterest <- Sys.time() - months(timewindow)
  data <- y
  data <- data |> filter(DateTime >= WindowOfInterest)
  TheDates <- data |> pull(DateTime) |> unique()

  Instrument <- x

  # x <- TheDates[1]

  InstrumentHistory <- map(.x=TheDates, data=data,
     .f=DateMapper, Instrument=Instrument) |> bind_rows()

  return(InstrumentHistory)
}